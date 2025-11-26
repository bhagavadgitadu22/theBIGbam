//! Shared processing logic for both CLI and Python bindings.

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use walkdir::WalkDir;

use crate::bam_reader::{detect_sequencing_type, process_reads_for_contig};
use crate::compress::compress_signal;
use crate::db::{create_metadata_db, create_temp_sample_db, write_features_to_temp_db, write_presences_to_temp_db, merge_temp_db_into_main, finalize_db};
use crate::features::{calculate_assemblycheck, calculate_coverage, calculate_phagetermini};
use crate::genbank::parse_genbank;
use crate::types::{get_plot_type, ContigInfo, FeaturePoint, PresenceData};

/// Configuration for processing.
pub struct ProcessConfig {
    pub threads: usize,
    pub min_coverage: f64,
    pub step: usize,
    pub z_thresh: f64,
    pub deriv_thresh: f64,
    pub max_points: usize,
}

impl Default for ProcessConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            min_coverage: 50.0,
            step: 50,
            z_thresh: 3.0,
            deriv_thresh: 3.0,
            max_points: 10000,
        }
    }
}

/// Result of processing all samples.
pub struct ProcessResult {
    pub samples_processed: usize,
    pub samples_failed: usize,
    pub total_time_secs: f64,
}

/// Compress and add feature points to the output vector.
fn add_compressed_feature(
    values: &[f64],
    feature: &str,
    contig_name: &str,
    config: &ProcessConfig,
    output: &mut Vec<FeaturePoint>,
) {
    let plot_type = get_plot_type(feature);
    let (xs, ys) = compress_signal(values, plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
    output.extend(xs.into_iter().zip(ys).map(|(x, y)| FeaturePoint {
        contig_name: contig_name.to_string(),
        feature: feature.to_string(),
        position: x,
        value: y,
    }));
}

/// Process one BAM file (sample) and return all features.
pub fn process_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
) -> Result<(Vec<FeaturePoint>, Vec<PresenceData>, String)> {
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "")
        .to_string();

    let seq_type = detect_sequencing_type(bam_path)?;

    let mut bam = IndexedReader::from_path(bam_path)?;
    bam.set_threads(2)?;

    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();

    for contig in contigs {
        let tid = bam.header().tid(contig.name.as_bytes());
        let ref_length = if let Some(tid) = tid {
            let bam_length = bam.header().target_len(tid).unwrap_or(contig.length as u64) as usize;
            bam_length / 2
        } else {
            contig.length
        };

        let reads = process_reads_for_contig(&mut bam, &contig.name, ref_length, modules, seq_type)?;

        if reads.is_empty() {
            continue;
        }

        // Coverage check
        let mut covered = vec![false; ref_length];
        for read in &reads {
            let start = (read.ref_start.max(0) as usize).min(ref_length);
            let end = (read.ref_end as usize).min(ref_length);
            for i in start..end {
                covered[i] = true;
            }
        }
        let covered_bp: usize = covered.iter().filter(|&&x| x).count();
        let coverage_pct = (covered_bp as f64 / ref_length as f64) * 100.0;

        if coverage_pct < config.min_coverage {
            continue;
        }

        all_presences.push(PresenceData {
            contig_name: contig.name.clone(),
            coverage_pct: coverage_pct as f32,
        });

        // Coverage
        if modules.contains(&"coverage".to_string()) || modules.contains(&"assemblycheck".to_string()) {
            let coverage = calculate_coverage(&reads, ref_length);
            let values: Vec<f64> = coverage.iter().map(|&x| x as f64).collect();
            add_compressed_feature(&values, "coverage", &contig.name, config, &mut all_features);
        }

        // Phagetermini
        if modules.contains(&"phagetermini".to_string()) {
            let pt_features = calculate_phagetermini(&reads, ref_length, seq_type);

            for feature_name in ["coverage_reduced", "reads_starts", "reads_ends"] {
                if let Some(data) = pt_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    add_compressed_feature(&values, feature_name, &contig.name, config, &mut all_features);
                }
            }

            // Tau (derived: (starts + ends) / coverage_reduced)
            if let (Some(cov_red), Some(starts), Some(ends)) = (
                pt_features.get("coverage_reduced"),
                pt_features.get("reads_starts"),
                pt_features.get("reads_ends"),
            ) {
                let tau: Vec<f64> = cov_red.iter().zip(starts).zip(ends)
                    .map(|((&c, &s), &e)| if c > 0 { (s + e) as f64 / c as f64 } else { 0.0 })
                    .collect();
                add_compressed_feature(&tau, "tau", &contig.name, config, &mut all_features);
            }
        }

        // Assembly check
        if modules.contains(&"assemblycheck".to_string()) {
            let ac_features = calculate_assemblycheck(&reads, ref_length, seq_type);

            for feature_name in ["left_clippings", "right_clippings", "insertions", "deletions", "mismatches", "bad_orientations"] {
                if let Some(data) = ac_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    add_compressed_feature(&values, feature_name, &contig.name, config, &mut all_features);
                }
            }

            // Read lengths (derived: sum / count, for long reads)
            if let (Some(sum_rl), Some(count_rl)) = (ac_features.get("sum_read_lengths"), ac_features.get("count_read_lengths")) {
                let values: Vec<f64> = sum_rl.iter().zip(count_rl)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                add_compressed_feature(&values, "read_lengths", &contig.name, config, &mut all_features);
            }

            // Insert sizes (derived: sum / count, for short-paired reads)
            if let (Some(sum_is), Some(count_is)) = (ac_features.get("sum_insert_sizes"), ac_features.get("count_insert_sizes")) {
                let values: Vec<f64> = sum_is.iter().zip(count_is)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                add_compressed_feature(&values, "insert_sizes", &contig.name, config, &mut all_features);
            }
        }
    }

    Ok((all_features, all_presences, sample_name))
}

/// Run processing on all BAM files in a directory.
/// This is the main entry point shared by both CLI and Python bindings.
pub fn run_all_samples(
    genbank_path: &Path,
    bam_dir: &Path,
    output_db: &Path,
    modules: &[String],
    annotation_tool: &str,
    config: &ProcessConfig,
) -> Result<ProcessResult> {
    // Suppress htslib warnings
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok(); // Ignore if already built

    eprintln!("### Parsing GenBank file...");
    let (contigs, annotations) = parse_genbank(genbank_path, annotation_tool)?;
    eprintln!("Found {} contigs with {} annotations", contigs.len(), annotations.len());

    // Get BAM files
    let bam_files: Vec<PathBuf> = if bam_dir.is_dir() {
        WalkDir::new(bam_dir)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().map(|x| x == "bam").unwrap_or(false))
            .map(|e| e.path().to_path_buf())
            .collect()
    } else {
        vec![bam_dir.to_path_buf()]
    };

    if bam_files.is_empty() {
        anyhow::bail!("No BAM files found");
    }

    // Sort by file size (largest first) for load balancing
    let mut bam_files = bam_files;
    bam_files.sort_by_key(|p| std::cmp::Reverse(fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    eprintln!("Found {} BAM files\n", bam_files.len());

    // Create parent directory if needed and temp directory for per-sample DBs
    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent)?;
    }
    let temp_dir = output_db.with_extension("temp_dbs");
    fs::create_dir_all(&temp_dir)?;

    std::thread::sleep(std::time::Duration::from_millis(50));

    // Create main database with schema
    let db_path = output_db;
    create_metadata_db(db_path, &contigs, &annotations)?;

    eprintln!("### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!("Modules: {}\n", modules.join(", "));

    // Create progress bar
    let pb = ProgressBar::new(bam_files.len() as u64);
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let start_time = std::time::Instant::now();
    let completed_count = AtomicUsize::new(0);
    let total = bam_files.len();

    // Process samples in parallel - each writes to its own temp DB
    let results: Vec<Result<(String, PathBuf, f64)>> = bam_files
        .par_iter()
        .map(|bam_path| {
            let sample_start = std::time::Instant::now();

            let result = process_sample(bam_path, &contigs, modules, config);

            match result {
                Ok((features, presences, sample_name)) => {
                    // Create temp DB for this sample
                    let temp_db_path = temp_dir.join(format!("{}.db", sample_name));
                    let temp_conn = create_temp_sample_db(&temp_db_path)?;

                    // Write features and presences to temp DB
                    if !features.is_empty() {
                        write_features_to_temp_db(&temp_conn, &features)?;
                    }
                    if !presences.is_empty() {
                        write_presences_to_temp_db(&temp_conn, &sample_name, &presences)?;
                    }

                    // Close temp DB connection
                    drop(temp_conn);

                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    pb.set_message(format!("[{}/{}] {} ({:.2}s)", done, total, sample_name, sample_time));
                    pb.inc(1);
                    Ok((sample_name, temp_db_path, sample_time))
                }
                Err(e) => {
                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    pb.set_message(format!("[{}/{}] ERR ({:.2}s)", done, total, sample_time));
                    pb.inc(1);
                    Err(e)
                }
            }
        })
        .collect();

    pb.finish_with_message("Done");

    // Merge all temp DBs into main DB sequentially (avoids lock contention)
    eprintln!("Merging sample data into database...");
    let merge_pb = ProgressBar::new(results.len() as u64);
    merge_pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} merging...")
            .unwrap()
            .progress_chars("=>-"),
    );

    for result in &results {
        if let Ok((_, temp_db_path, _)) = result {
            merge_temp_db_into_main(&db_path, temp_db_path, &contigs)?;
            // Clean up temp DB
            let _ = fs::remove_file(temp_db_path);
        }
        merge_pb.inc(1);
    }
    merge_pb.finish_with_message("Done");

    // Clean up temp directory
    let _ = fs::remove_dir(&temp_dir);

    // Finalize database (checkpoint WAL)
    finalize_db(&db_path)?;

    let elapsed = start_time.elapsed();
    let successful = results.iter().filter(|r| r.is_ok()).count();
    let failed = results.len() - successful;

    // Print summary
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!("  Samples processed: {}/{}", successful, results.len());
    eprintln!("  Total time:        {:.2}s", elapsed.as_secs_f64());
    eprintln!();
    eprintln!("  Output: {:?}", output_db);
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    Ok(ProcessResult {
        samples_processed: successful,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
    })
}
