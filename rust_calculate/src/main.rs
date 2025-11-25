//! # MGFeatureViewer Rust Calculator
//!
//! This is a 1:1 Rust port of the Python calculation code in `mgfeatureviewer/calculating_data.py`.
//! It produces identical output values, with one minor difference noted below.
//!
//! ## Module Structure
//!
//! - `types.rs`      - All structs and enums
//! - `bam_reader.rs` - BAM file reading and sequencing type detection
//! - `features.rs`   - Feature calculations (coverage, phagetermini, assemblycheck)
//! - `compress.rs`   - Signal compression for storage
//! - `db.rs`         - SQLite database operations
//! - `parquet.rs`    - Parquet file writing
//! - `genbank.rs`    - GenBank file parsing
//!
//! ## Function Mapping (Rust → Python)
//!
//! | Rust Function                    | Python Function                              | Line  |
//! |----------------------------------|----------------------------------------------|-------|
//! | `detect_sequencing_type`         | `find_sequencing_type_from_bam`              | ~55   |
//! | `process_reads_for_contig`       | `preprocess_reads`                           | ~497  |
//! | `calculate_coverage`             | `calculate_coverage_numba`                   | ~182  |
//! | `calculate_phagetermini`         | `get_features_phagetermini`                  | ~283  |
//! | `calculate_assemblycheck`        | `get_features_assemblycheck`                 | ~410  |
//! | `compress_signal`                | `compress_signal`                            | ~91   |
//! | `process_sample`                 | `calculating_features_per_sample`            | ~651  |
//!
//! ## Known Difference: Derivative Outlier Detection
//!
//! The Python `compress_signal` function (line ~116) has a subtle bug in derivative outlier
//! detection where it performs arithmetic on a boolean array. This causes Python to keep
//! fewer positions after compression. The Rust implementation correctly identifies derivative
//! outliers, so it may keep slightly more positions. **All values at common positions are
//! identical** - Rust just preserves more data points as originally intended.
//!
//! ## Verification
//!
//! Run `python compare_values_only.py <sqlite_db> <parquet_dir>` to verify outputs match:
//! - All 11 features tested across all samples
//! - 100% exact match at common positions (0 difference)
//! - Metadata tables match exactly
//!
//! See also: `PYTHON_COMPARISON.md` for side-by-side code examples.

mod bam_reader;
mod compress;
mod db;
mod features;
mod genbank;
mod parquet;
mod types;

use anyhow::Result;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::bam_reader::{detect_sequencing_type, process_reads_for_contig};
use crate::compress::compress_signal;
use crate::db::{create_metadata_db, update_sample_presences};
use crate::features::{calculate_assemblycheck, calculate_coverage, calculate_phagetermini};
use crate::genbank::parse_genbank;
use crate::parquet::{write_features_parquet, write_presences_parquet};
use crate::types::{get_feature_config, ContigInfo, FeaturePoint, PresenceData};

#[derive(Parser, Debug)]
#[command(name = "mgfeatureviewer_calc")]
#[command(about = "Calculate features from BAM files (Rust implementation)")]
struct Args {
    /// Number of threads
    #[arg(short = 't', long)]
    threads: usize,

    /// Path to GenBank file
    #[arg(short = 'g', long)]
    genbank: PathBuf,

    /// Path to BAM file or directory
    #[arg(short = 'b', long)]
    bam_files: PathBuf,

    /// Modules to compute (comma-separated: coverage,phagetermini,assemblycheck)
    #[arg(short = 'm', long)]
    modules: String,

    /// Output directory
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Annotation tool (optional)
    #[arg(short = 'a', long, default_value = "")]
    annotation_tool: String,

    /// Minimum coverage percentage for contig inclusion
    #[arg(long, default_value = "50")]
    min_coverage: f64,

    /// Step size for compression
    #[arg(long, default_value = "50")]
    step: usize,

    /// Z-score threshold for outliers
    #[arg(long, default_value = "3.0")]
    outlier_threshold: f64,

    /// Derivative threshold for outliers
    #[arg(long, default_value = "3.0")]
    derivative_threshold: f64,

    /// Maximum points after compression
    #[arg(long, default_value = "10000")]
    max_points: usize,
}

/// Process one BAM file (sample) and return all features.
/// Python equivalent: `calculating_features_per_sample()` in calculating_data.py:651-673
/// Also calls `calculating_features_per_contig_per_sample()` at py:605-647
fn process_sample(
    bam_path: &std::path::Path,
    contigs: &[ContigInfo],
    modules: &[String],
    min_coverage: f64,
    step: usize,
    z_thresh: f64,
    deriv_thresh: f64,
    max_points: usize,
) -> Result<(Vec<FeaturePoint>, Vec<PresenceData>, String)> {
    // calculating_data.py:651-673 - calculating_features_per_sample()

    // py:657 - sample_name = os.path.basename(mapping_file).replace(".bam", "")
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "")
        .to_string();

    // py:658 - sequencing_type = find_sequencing_type_from_bam(mapping_file)
    let seq_type = detect_sequencing_type(bam_path)?;

    // py:653 - bam_file = pysam.AlignmentFile(mapping_file, "rb")
    let mut bam = IndexedReader::from_path(bam_path)?;
    // Enable multi-threaded BAM decompression (2 threads per file)
    bam.set_threads(2)?;

    let mut all_features = Vec::new();  // py:660
    let mut all_presences = Vec::new(); // py:661

    // py:663 - for ref, length in zip(references, lengths):
    for contig in contigs {
        // py:655 - lengths = [l // 2 for l in bam_file.lengths]
        // Python uses BAM header length / 2 for calculations (circular genome handling)
        let tid = bam.header().tid(contig.name.as_bytes());
        let ref_length = if let Some(tid) = tid {
            let bam_length = bam.header().target_len(tid).unwrap_or(contig.length as u64) as usize;
            bam_length / 2  // Match Python's circular genome handling
        } else {
            contig.length  // Fallback to GenBank length if not in BAM
        };

        // py:609-612 - Fetch reads and preprocess
        let reads = process_reads_for_contig(&mut bam, &contig.name, ref_length, modules, seq_type)?;

        // py:615-616 - if len(ref_starts) == 0: return None, None
        if reads.is_empty() {
            continue;
        }

        // py:619-626 - Coverage check: calculate coverage percentage
        let mut covered = vec![false; ref_length];  // py:619
        for read in &reads {
            let start = (read.ref_start.max(0) as usize).min(ref_length);  // py:620
            let end = (read.ref_end as usize).min(ref_length);              // py:621
            for i in start..end {  // py:622-623
                covered[i] = true;
            }
        }
        let covered_bp: usize = covered.iter().filter(|&&x| x).count();  // py:625
        let coverage_pct = (covered_bp as f64 / ref_length as f64) * 100.0;  // py:626

        // py:628-629 - if coverage_pct < min_coverage: return None, None
        if coverage_pct < min_coverage {
            continue;
        }

        // py:631-632 - presence = (ref_name, coverage_pct)
        all_presences.push(PresenceData {
            contig_name: contig.name.clone(),
            coverage_pct: coverage_pct as f32,
        });

        // py:638-639 - if {"coverage", "assemblycheck"}.intersection(module_list):
        if modules.contains(&"coverage".to_string()) || modules.contains(&"assemblycheck".to_string()) {
            let coverage = calculate_coverage(&reads, ref_length);
            let values: Vec<f64> = coverage.iter().map(|&x| x as f64).collect();
            let config = get_feature_config("coverage");
            let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
            for (x, y) in xs.into_iter().zip(ys) {
                all_features.push(FeaturePoint {
                    contig_name: contig.name.clone(),
                    feature: "coverage".to_string(),
                    position: x,
                    value: y,
                });
            }
        }

        // py:640-642 - if "phagetermini" in module_list:
        if modules.contains(&"phagetermini".to_string()) {
            let pt_features = calculate_phagetermini(&reads, ref_length, seq_type);

            for feature_name in ["coverage_reduced", "reads_starts", "reads_ends"] {
                if let Some(data) = pt_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let config = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
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

            // py:278-281 - tau = (reads_starts + reads_ends) / coverage_reduced
            if let (Some(cov_red), Some(starts), Some(ends)) = (
                pt_features.get("coverage_reduced"),
                pt_features.get("reads_starts"),
                pt_features.get("reads_ends"),
            ) {
                let tau: Vec<f64> = cov_red
                    .iter()
                    .zip(starts)
                    .zip(ends)
                    .map(|((&c, &s), &e)| {
                        if c > 0 {
                            (s + e) as f64 / c as f64
                        } else {
                            0.0
                        }
                    })
                    .collect();
                let config = get_feature_config("tau");
                let (xs, ys) = compress_signal(&tau, config.plot_type, step, z_thresh, deriv_thresh, max_points);
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

        // py:643-645 - if "assemblycheck" in module_list:
        if modules.contains(&"assemblycheck".to_string()) {
            let ac_features = calculate_assemblycheck(&reads, ref_length, seq_type);

            // Process each assemblycheck feature
            for feature_name in ["left_clippings", "right_clippings", "insertions", "deletions", "mismatches", "bad_orientations"] {
                if let Some(data) = ac_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let config = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
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

            // py:484-485 - compute_final_lengths for read_lengths (long reads)
            if let (Some(sum_rl), Some(count_rl)) = (
                ac_features.get("sum_read_lengths"),
                ac_features.get("count_read_lengths"),
            ) {
                // py:404-408 - compute_final_lengths()
                let values: Vec<f64> = sum_rl
                    .iter()
                    .zip(count_rl)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let config = get_feature_config("read_lengths");
                let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "read_lengths".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }

            // py:486-487 - compute_final_lengths for insert_sizes (short-paired)
            if let (Some(sum_is), Some(count_is)) = (
                ac_features.get("sum_insert_sizes"),
                ac_features.get("count_insert_sizes"),
            ) {
                // py:404-408 - compute_final_lengths()
                let values: Vec<f64> = sum_is
                    .iter()
                    .zip(count_is)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let config = get_feature_config("insert_sizes");
                let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
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

    // py:673 - return all_features, all_presences, sample_name
    Ok((all_features, all_presences, sample_name))
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Suppress htslib warnings about index file timestamps
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    println!("### Parsing GenBank file...");
    let (contigs, annotations) = parse_genbank(&args.genbank, &args.annotation_tool)?;
    println!("Found {} contigs with {} annotations", contigs.len(), annotations.len());

    // Get BAM files
    let bam_files: Vec<PathBuf> = if args.bam_files.is_dir() {
        walkdir::WalkDir::new(&args.bam_files)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().map(|x| x == "bam").unwrap_or(false))
            .map(|e| e.path().to_path_buf())
            .collect()
    } else {
        vec![args.bam_files.clone()]
    };

    if bam_files.is_empty() {
        anyhow::bail!("No BAM files found");
    }

    // Sort by file size (largest first) for load balancing
    let mut bam_files = bam_files;
    bam_files.sort_by_key(|p| std::cmp::Reverse(fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    println!("Found {} BAM files", bam_files.len());

    // Create output directories
    if args.output.exists() {
        anyhow::bail!("Output directory already exists: {:?}", args.output);
    }
    fs::create_dir_all(&args.output)?;
    fs::create_dir_all(args.output.join("features"))?;
    fs::create_dir_all(args.output.join("presences"))?;

    // Sync to ensure directories are created before SQLite
    std::thread::sleep(std::time::Duration::from_millis(100));

    // Create metadata database - use canonical path
    let output_dir = fs::canonicalize(&args.output)?;
    let db_path = output_dir.join("metadata.db");
    create_metadata_db(&db_path, &contigs, &annotations)?;

    // Parse modules
    let modules: Vec<String> = args.modules.split(',').map(|s| s.trim().to_string()).collect();

    println!("\n### Processing {} samples with {} threads", bam_files.len(), args.threads);
    println!("Modules: {}", modules.join(", "));
    println!();

    // Create progress bar
    let pb = ProgressBar::new(bam_files.len() as u64);
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let start_time = std::time::Instant::now();
    let error_count = AtomicUsize::new(0);
    let completed_count = AtomicUsize::new(0);
    let total = bam_files.len();

    // Process samples in parallel - collect presences for SQLite
    let results: Vec<Result<(String, Vec<PresenceData>, f64)>> = bam_files
        .par_iter()
        .map(|bam_path| {
            let sample_start = std::time::Instant::now();

            let result = process_sample(
                bam_path,
                &contigs,
                &modules,
                args.min_coverage,
                args.step,
                args.outlier_threshold,
                args.derivative_threshold,
                args.max_points,
            );

            match result {
                Ok((features, presences, sample_name)) => {
                    // Write Parquet files
                    if !features.is_empty() {
                        let features_path = args.output.join("features").join(format!("{}.parquet", sample_name));
                        write_features_parquet(&features, &features_path)?;
                    }

                    if !presences.is_empty() {
                        let presences_path = args.output.join("presences").join(format!("{}.parquet", sample_name));
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
                    pb.set_message(format!("[{}/{}] ERR: {} ({:.2}s)", done, total, bam_path.file_name().unwrap_or_default().to_string_lossy(), sample_time));
                    pb.inc(1);
                    error_count.fetch_add(1, Ordering::Relaxed);
                    Err(e)
                }
            }
        })
        .collect();

    pb.finish_with_message("Done");

    // Insert Sample and Presences into SQLite
    println!("Updating metadata database...");
    update_sample_presences(&db_path, &results, &contigs)?;

    let elapsed = start_time.elapsed();
    let total_samples = bam_files.len();
    let successful = total_samples - error_count.load(Ordering::Relaxed);

    // Calculate statistics
    let sample_times: Vec<f64> = results.iter().filter_map(|r| r.as_ref().ok().map(|(_, _, t)| *t)).collect();
    let avg_time = if !sample_times.is_empty() {
        sample_times.iter().sum::<f64>() / sample_times.len() as f64
    } else {
        0.0
    };

    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("### Complete");
    println!();
    println!("  Samples processed: {}/{}", successful, total_samples);
    println!("  Total time:        {:.2}s", elapsed.as_secs_f64());
    println!("  Avg time/sample:   {:.2}s", avg_time);
    println!();
    println!("  Output: {:?}", args.output);
    println!("    - metadata.db");
    println!("    - features/*.parquet ({} files)", successful);
    println!("    - presences/*.parquet ({} files)", successful);
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    if error_count.load(Ordering::Relaxed) > 0 {
        eprintln!("\nWarning: {} samples failed to process", error_count.load(Ordering::Relaxed));
    }

    Ok(())
}
