//! Shared processing logic for both CLI and Python bindings.

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use walkdir::WalkDir;

use crate::bam_reader::{detect_sequencing_type, process_reads_for_contig};
use crate::compress::compress_signal;
use crate::db::{create_metadata_db, update_sample_presences};
use crate::features::{calculate_assemblycheck, calculate_coverage, calculate_phagetermini};
use crate::genbank::parse_genbank;
use crate::parquet::{write_features_parquet, write_presences_parquet};
use crate::types::{get_feature_config, ContigInfo, FeaturePoint, PresenceData};

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

/// Process one BAM file (sample) and return all features.
/// Python equivalent: `calculating_features_per_sample()` in calculating_data.py:651-673
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
            let cfg = get_feature_config("coverage");
            let (xs, ys) = compress_signal(&values, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
            for (x, y) in xs.into_iter().zip(ys) {
                all_features.push(FeaturePoint {
                    contig_name: contig.name.clone(),
                    feature: "coverage".to_string(),
                    position: x,
                    value: y,
                });
            }
        }

        // Phagetermini
        if modules.contains(&"phagetermini".to_string()) {
            let pt_features = calculate_phagetermini(&reads, ref_length, seq_type);

            for feature_name in ["coverage_reduced", "reads_starts", "reads_ends"] {
                if let Some(data) = pt_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let cfg = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push(FeaturePoint {
                            contig_name: contig.name.clone(),
                            feature: feature_name.to_string(),
                            position: x,
                            value: y,
                        });
                    }
                }
            }

            // Tau
            if let (Some(cov_red), Some(starts), Some(ends)) = (
                pt_features.get("coverage_reduced"),
                pt_features.get("reads_starts"),
                pt_features.get("reads_ends"),
            ) {
                let tau: Vec<f64> = cov_red
                    .iter()
                    .zip(starts)
                    .zip(ends)
                    .map(|((&c, &s), &e)| if c > 0 { (s + e) as f64 / c as f64 } else { 0.0 })
                    .collect();
                let cfg = get_feature_config("tau");
                let (xs, ys) = compress_signal(&tau, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "tau".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }
        }

        // Assembly check
        if modules.contains(&"assemblycheck".to_string()) {
            let ac_features = calculate_assemblycheck(&reads, ref_length, seq_type);

            for feature_name in ["left_clippings", "right_clippings", "insertions", "deletions", "mismatches", "bad_orientations"] {
                if let Some(data) = ac_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let cfg = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push(FeaturePoint {
                            contig_name: contig.name.clone(),
                            feature: feature_name.to_string(),
                            position: x,
                            value: y,
                        });
                    }
                }
            }

            // Read lengths (long reads)
            if let (Some(sum_rl), Some(count_rl)) = (
                ac_features.get("sum_read_lengths"),
                ac_features.get("count_read_lengths"),
            ) {
                let values: Vec<f64> = sum_rl
                    .iter()
                    .zip(count_rl)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let cfg = get_feature_config("read_lengths");
                let (xs, ys) = compress_signal(&values, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "read_lengths".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }

            // Insert sizes (short-paired)
            if let (Some(sum_is), Some(count_is)) = (
                ac_features.get("sum_insert_sizes"),
                ac_features.get("count_insert_sizes"),
            ) {
                let values: Vec<f64> = sum_is
                    .iter()
                    .zip(count_is)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let cfg = get_feature_config("insert_sizes");
                let (xs, ys) = compress_signal(&values, cfg.plot_type, config.step, config.z_thresh, config.deriv_thresh, config.max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "insert_sizes".to_string(),
                        position: x,
                        value: y,
                    });
                }
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
    output_dir: &Path,
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
    let bam_files: Vec<std::path::PathBuf> = if bam_dir.is_dir() {
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

    // Create output directories
    fs::create_dir_all(output_dir)?;
    fs::create_dir_all(output_dir.join("features"))?;
    fs::create_dir_all(output_dir.join("presences"))?;

    std::thread::sleep(std::time::Duration::from_millis(50));

    // Create metadata database
    let db_path = output_dir.join("metadata.db");
    create_metadata_db(&db_path, &contigs, &annotations)?;

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

    // Process samples in parallel
    let results: Vec<Result<(String, Vec<PresenceData>, f64)>> = bam_files
        .par_iter()
        .map(|bam_path| {
            let sample_start = std::time::Instant::now();

            let result = process_sample(bam_path, &contigs, modules, config);

            match result {
                Ok((features, presences, sample_name)) => {
                    if !features.is_empty() {
                        let features_path = output_dir.join("features").join(format!("{}.parquet", sample_name));
                        write_features_parquet(&features, &features_path)?;
                    }

                    if !presences.is_empty() {
                        let presences_path = output_dir.join("presences").join(format!("{}.parquet", sample_name));
                        write_presences_parquet(&presences, &presences_path)?;
                    }

                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    pb.set_message(format!("[{}/{}] {} ({:.2}s)", done, total, sample_name, sample_time));
                    pb.inc(1);
                    Ok((sample_name, presences, sample_time))
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

    // Update SQLite with sample/presences
    eprintln!("Updating metadata database...");
    update_sample_presences(&db_path, &results, &contigs)?;

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
    eprintln!("  Output: {:?}", output_dir);
    eprintln!("    - metadata.db");
    eprintln!("    - features/*.parquet ({} files)", successful);
    eprintln!("    - presences/*.parquet ({} files)", successful);
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    Ok(ProcessResult {
        samples_processed: successful,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
    })
}
