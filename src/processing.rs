//! Shared processing logic for both CLI and Python bindings.

use anyhow::Result;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
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
    /// If true, replicate Python's buggy derivative outlier detection for comparison
    pub python_compat: bool,
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
            python_compat: false,
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
    let (xs, ys) = compress_signal(values, plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points, config.python_compat);
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

    // Iterate over BAM header references (like Python does) instead of GenBank contigs
    // This ensures we only process contigs that are in the BAM file
    let header = bam.header().clone();
    let n_refs = header.target_count();

    for tid in 0..n_refs {
        let ref_name = std::str::from_utf8(header.tid2name(tid)).unwrap_or("");
        let bam_length = header.target_len(tid).unwrap_or(0) as usize;
        let ref_length = bam_length / 2;  // py: lengths = [l // 2 for l in bam.lengths]

        // Skip if this contig is not in our GenBank contigs list
        let contig = match contigs.iter().find(|c| c.name == ref_name) {
            Some(c) => c,
            None => continue,
        };

        let reads = process_reads_for_contig(&mut bam, ref_name, ref_length, modules, seq_type)?;

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

    // Create main database with schema (skip indexes in python_compat mode)
    let db_path = output_db;
    create_metadata_db(db_path, &contigs, &annotations, !config.python_compat)?;

    eprintln!("### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!("Modules: {}\n", modules.join(", "));

    // Create multi-progress with two bars
    let mp = MultiProgress::new();

    let process_pb = mp.add(ProgressBar::new(bam_files.len() as u64));
    process_pb.set_style(
        ProgressStyle::with_template("{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );
    process_pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let merge_pb = mp.add(ProgressBar::new(bam_files.len() as u64));
    merge_pb.set_style(
        ProgressStyle::with_template("{spinner:.yellow} Writing:     [{elapsed_precise}] [{bar:40.yellow/red}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );
    merge_pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let merge_pb = Arc::new(merge_pb);

    let start_time = std::time::Instant::now();
    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);
    let total = bam_files.len();

    // Create channel for background merging - samples are merged as they complete
    let (merge_tx, merge_rx) = mpsc::channel::<PathBuf>();

    // Clone data needed by merge thread
    let merge_db_path = db_path.to_path_buf();
    let merge_contigs = contigs.clone();
    let merge_temp_dir = temp_dir.clone();
    let merge_pb_clone = Arc::clone(&merge_pb);

    // Spawn background merge thread
    let merge_handle = thread::spawn(move || {
        let mut merged = 0usize;
        for temp_db_path in merge_rx {
            let sample_name = temp_db_path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();

            if merge_temp_db_into_main(&merge_db_path, &temp_db_path, &merge_contigs).is_err() {
                merge_pb_clone.set_message(format!("ERR: {}", sample_name));
            } else {
                merge_pb_clone.set_message(format!("{}", sample_name));
            }
            // Clean up temp DB
            let _ = fs::remove_file(&temp_db_path);
            merged += 1;
            merge_pb_clone.inc(1);
        }
        // Clean up temp directory
        let _ = fs::remove_dir(&merge_temp_dir);
        merge_pb_clone.finish_with_message("Done");
        merged
    });

    // Process samples in parallel - each writes to its own temp DB and sends to merge thread
    bam_files
        .par_iter()
        .for_each(|bam_path| {
            let sample_start = std::time::Instant::now();

            let result = process_sample(bam_path, &contigs, modules, config);

            match result {
                Ok((features, presences, sample_name)) => {
                    // Create temp DB for this sample
                    let temp_db_path = temp_dir.join(format!("{}.db", sample_name));
                    let temp_conn = match create_temp_sample_db(&temp_db_path) {
                        Ok(c) => c,
                        Err(_) => {
                            failed_count.fetch_add(1, Ordering::SeqCst);
                            return;
                        }
                    };

                    // Write features and presences to temp DB
                    if !features.is_empty() {
                        if write_features_to_temp_db(&temp_conn, &features).is_err() {
                            failed_count.fetch_add(1, Ordering::SeqCst);
                            return;
                        }
                    }
                    if !presences.is_empty() {
                        if write_presences_to_temp_db(&temp_conn, &sample_name, &presences).is_err() {
                            failed_count.fetch_add(1, Ordering::SeqCst);
                            return;
                        }
                    }

                    // Close temp DB connection
                    drop(temp_conn);

                    // Send to merge thread (non-blocking)
                    let _ = merge_tx.send(temp_db_path);

                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    process_pb.set_message(format!("[{}/{}] {} ({:.2}s)", done, total, sample_name, sample_time));
                    process_pb.inc(1);
                }
                Err(_) => {
                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    failed_count.fetch_add(1, Ordering::SeqCst);
                    process_pb.set_message(format!("[{}/{}] ERR ({:.2}s)", done, total, sample_time));
                    process_pb.inc(1);
                }
            }
        });

    process_pb.finish_with_message("Done");

    // Close channel to signal merge thread to finish
    drop(merge_tx);

    // Wait for merge thread to complete
    let _merged = merge_handle.join().expect("Merge thread panicked");

    // Finalize database (checkpoint WAL)
    finalize_db(db_path)?;

    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);
    let successful = total - failed;

    // Print summary
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!("  Samples processed: {}/{}", successful, total);
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
