//! Shared processing logic for both CLI and Python bindings.
//!
//! # Overview
//!
//! This module orchestrates the parallel processing of BAM files:
//! 1. Parse GenBank file to get contig metadata
//! 2. Discover BAM files in input directory
//! 3. Process each BAM file in parallel using rayon
//! 4. Write results to SQLite database
//!
//! The processing is optimized for single-pass streaming of BAM files,
//! avoiding the need to load all reads into memory.
//!
//! # Python Equivalent
//!
//! This module corresponds to the main processing functions in `calculating_data.py`:
//! ```python
//! # Python: calculating_all_features_parallel() - parallel sample processing
//! def calculating_all_features_parallel(list_modules, bam_files, db_path, ...):
//!     with Pool(processes=n_sample_cores) as pool:
//!         pool.starmap(_process_single_sample, args_list)
//!
//! # Python: _process_single_sample() - single sample processing
//! def _process_single_sample(list_modules, bam_file, db_path, ...):
//!     # Create temp DB, process all contigs, merge into main DB
//!     calculating_features_per_sample(...)
//!     merge_temp_db_into_main(db_main, str(temp_db))
//! ```
//!
//! # Architecture Differences
//!
//! ## Python
//! - Uses multiprocessing.Pool for parallelism
//! - Each process opens its own BAM file and database connections
//! - Temp DBs are merged sequentially after all processing
//!
//! ## Rust
//! - Uses rayon for thread-based parallelism (lower overhead than processes)
//! - Dedicated merge thread receives temp DBs via channel and merges incrementally
//! - Progress bars show both processing and merge progress
//!
//! # Rust Concepts
//!
//! ## Rayon
//! `bam_files.par_iter().for_each(...)` parallelizes iteration over BAM files.
//! Rayon automatically distributes work across threads.
//!
//! ## Channels (`mpsc`)
//! `mpsc::channel()` creates a multi-producer, single-consumer channel.
//! Processing threads send temp DB paths to a dedicated merge thread.
//!
//! ## Atomics
//! `AtomicUsize` provides thread-safe counters without mutex overhead.
//! Used for tracking completed/failed samples.

use anyhow::{Context, Result};
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

use crate::bam_reader::{detect_sequencing_type, process_contig_streaming};
use crate::compress::compress_signal;
use crate::db::{
    create_metadata_db, create_temp_sample_db, finalize_db, merge_temp_db_into_main,
    write_features_to_temp_db, write_presences_to_temp_db,
};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::genbank::parse_genbank;
use crate::types::{
    get_plot_type, ContigInfo, FeaturePoint, PresenceData, ASSEMBLYCHECK_FEATURES,
    PHAGETERMINI_FEATURES,
};

/// Configuration for processing.
#[derive(Clone)]
pub struct ProcessConfig {
    pub threads: usize,
    pub min_coverage: f64,
    pub step: usize,
    pub z_thresh: f64,
    pub deriv_thresh: f64,
    pub max_points: usize,
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
#[inline]
fn add_compressed_feature(
    values: &[f64],
    feature: &str,
    contig_name: &str,
    config: &ProcessConfig,
    output: &mut Vec<FeaturePoint>,
) {
    let plot_type = get_plot_type(feature);
    let (xs, ys) = compress_signal(
        values,
        plot_type,
        config.step,
        config.z_thresh,
        config.deriv_thresh,
        config.max_points,
        config.python_compat,
    );

    output.extend(xs.into_iter().zip(ys).map(|(x, y)| FeaturePoint {
        contig_name: contig_name.to_string(),
        feature: feature.to_string(),
        position: x,
        value: y,
    }));
}

/// Add features from FeatureArrays to output (optimized path).
fn add_features_from_arrays(
    arrays: &FeatureArrays,
    contig_name: &str,
    config: &ProcessConfig,
    flags: ModuleFlags,
    seq_type: crate::types::SequencingType,
    output: &mut Vec<FeaturePoint>,
) {
    // Coverage
    if flags.coverage {
        let values: Vec<f64> = arrays.coverage.iter().map(|&x| x as f64).collect();
        add_compressed_feature(&values, "coverage", contig_name, config, output);
    }

    // Phagetermini features
    if flags.phagetermini {
        for &feature_name in PHAGETERMINI_FEATURES {
            let data = match feature_name {
                "coverage_reduced" => &arrays.coverage_reduced,
                "reads_starts" => &arrays.reads_starts,
                "reads_ends" => &arrays.reads_ends,
                _ => continue,
            };
            let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
            add_compressed_feature(&values, feature_name, contig_name, config, output);
        }

        // Tau (derived)
        let tau: Vec<f64> = arrays
            .coverage_reduced
            .iter()
            .zip(&arrays.reads_starts)
            .zip(&arrays.reads_ends)
            .map(|((&c, &s), &e)| {
                if c > 0 {
                    (s + e) as f64 / c as f64
                } else {
                    0.0
                }
            })
            .collect();
        add_compressed_feature(&tau, "tau", contig_name, config, output);
    }

    // Assemblycheck features
    if flags.assemblycheck {
        for &feature_name in ASSEMBLYCHECK_FEATURES {
            let data = match feature_name {
                "left_clippings" => &arrays.left_clippings,
                "right_clippings" => &arrays.right_clippings,
                "insertions" => &arrays.insertions,
                "deletions" => &arrays.deletions,
                "mismatches" => &arrays.mismatches,
                "bad_orientations" => &arrays.bad_orientations,
                _ => continue,
            };
            let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
            add_compressed_feature(&values, feature_name, contig_name, config, output);
        }

        // Read lengths (derived, for long reads)
        if seq_type.is_long() {
            let values: Vec<f64> = arrays
                .sum_read_lengths
                .iter()
                .zip(&arrays.count_read_lengths)
                .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                .collect();
            add_compressed_feature(&values, "read_lengths", contig_name, config, output);
        }

        // Insert sizes (derived, for short-paired)
        if seq_type.is_short_paired() {
            let values: Vec<f64> = arrays
                .sum_insert_sizes
                .iter()
                .zip(&arrays.count_insert_sizes)
                .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                .collect();
            add_compressed_feature(&values, "insert_sizes", contig_name, config, output);
        }
    }
}

/// Process one BAM file using streaming (single-pass, optimized).
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
        .replace("_with_MD", "");

    let seq_type = detect_sequencing_type(bam_path)?;
    let flags = ModuleFlags::from_modules(modules);

    let mut bam = IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;
    // 4 decompression threads: benchmarked optimal for parallel BAM processing (7% faster than 2)
    bam.set_threads(4.min(config.threads.max(2)))?;

    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();

    let header = bam.header().clone();
    let n_refs = header.target_count();

    for tid in 0..n_refs {
        let ref_name = std::str::from_utf8(header.tid2name(tid)).unwrap_or("");
        let bam_length = header.target_len(tid).unwrap_or(0) as usize;
        let ref_length = bam_length / 2;

        // Skip if contig not in GenBank list
        let contig = match contigs.iter().find(|c| c.name == ref_name) {
            Some(c) => c,
            None => continue,
        };

        // Process contig using streaming - single pass over reads
        let arrays = match process_contig_streaming(&mut bam, ref_name, ref_length, seq_type, flags)? {
            Some(a) => a,
            None => continue,
        };

        // Check coverage percentage
        let coverage_pct = arrays.coverage_percentage();
        if coverage_pct < config.min_coverage {
            continue;
        }

        all_presences.push(PresenceData {
            contig_name: contig.name.clone(),
            coverage_pct: coverage_pct as f32,
        });

        // Add features from arrays
        add_features_from_arrays(&arrays, &contig.name, config, flags, seq_type, &mut all_features);
    }

    Ok((all_features, all_presences, sample_name))
}

/// Run processing on all BAM files in a directory.
pub fn run_all_samples(
    genbank_path: &Path,
    bam_dir: &Path,
    output_db: &Path,
    modules: &[String],
    annotation_tool: &str,
    config: &ProcessConfig,
) -> Result<ProcessResult> {
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    eprintln!("### Parsing GenBank file...");
    let (contigs, annotations) = parse_genbank(genbank_path, annotation_tool)?;
    eprintln!(
        "Found {} contigs with {} annotations",
        contigs.len(),
        annotations.len()
    );

    let bam_files = discover_bam_files(bam_dir)?;
    eprintln!("Found {} BAM files\n", bam_files.len());

    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent).context("Failed to create output directory")?;
    }
    let temp_dir = output_db.with_extension("temp_dbs");
    fs::create_dir_all(&temp_dir).context("Failed to create temp directory")?;

    std::thread::sleep(std::time::Duration::from_millis(50));

    create_metadata_db(output_db, &contigs, &annotations, !config.python_compat)?;

    eprintln!(
        "### Processing {} samples with {} threads",
        bam_files.len(),
        config.threads
    );
    eprintln!("Modules: {}\n", modules.join(", "));

    let result = process_samples_parallel(&bam_files, &contigs, modules, config, output_db, &temp_dir)?;

    finalize_db(output_db)?;

    print_summary(&result, output_db);

    Ok(result)
}

fn discover_bam_files(bam_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut bam_files: Vec<PathBuf> = if bam_dir.is_dir() {
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
        anyhow::bail!("No BAM files found in {}", bam_dir.display());
    }

    bam_files.sort_by_key(|p| std::cmp::Reverse(fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    Ok(bam_files)
}

fn process_samples_parallel(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    db_path: &Path,
    temp_dir: &Path,
) -> Result<ProcessResult> {
    let mp = MultiProgress::new();
    let total = bam_files.len();

    let process_pb = mp.add(ProgressBar::new(total as u64));
    process_pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
        )
        .unwrap()
        .progress_chars("=>-"),
    );
    process_pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let merge_pb = mp.add(ProgressBar::new(total as u64));
    merge_pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.yellow} Writing:     [{elapsed_precise}] [{bar:40.yellow/red}] {pos}/{len} {msg}",
        )
        .unwrap()
        .progress_chars("=>-"),
    );
    merge_pb.enable_steady_tick(std::time::Duration::from_millis(100));
    let merge_pb = Arc::new(merge_pb);

    let start_time = std::time::Instant::now();
    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);

    let (merge_tx, merge_rx) = mpsc::channel::<PathBuf>();

    let merge_handle = spawn_merge_thread(
        merge_rx,
        db_path.to_path_buf(),
        contigs.to_vec(),
        temp_dir.to_path_buf(),
        Arc::clone(&merge_pb),
    );

    bam_files.par_iter().for_each(|bam_path| {
        process_single_sample(
            bam_path,
            contigs,
            modules,
            config,
            temp_dir,
            &merge_tx,
            &completed_count,
            &failed_count,
            &process_pb,
            total,
        );
    });

    process_pb.finish_with_message("Done");
    drop(merge_tx);

    merge_handle.join().expect("Merge thread panicked");

    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);

    Ok(ProcessResult {
        samples_processed: total - failed,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
    })
}

#[allow(clippy::too_many_arguments)]
fn process_single_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    temp_dir: &Path,
    merge_tx: &mpsc::Sender<PathBuf>,
    completed_count: &AtomicUsize,
    failed_count: &AtomicUsize,
    process_pb: &ProgressBar,
    total: usize,
) {
    let sample_start = std::time::Instant::now();
    let result = process_sample(bam_path, contigs, modules, config);

    match result {
        Ok((features, presences, sample_name)) => {
            let temp_db_path = temp_dir.join(format!("{}.db", sample_name));

            if let Ok(temp_conn) = create_temp_sample_db(&temp_db_path) {
                let write_ok = (!features.is_empty()
                    && write_features_to_temp_db(&temp_conn, &features).is_err())
                    || (!presences.is_empty()
                        && write_presences_to_temp_db(&temp_conn, &sample_name, &presences).is_err());

                if write_ok {
                    failed_count.fetch_add(1, Ordering::SeqCst);
                    return;
                }

                drop(temp_conn);
                let _ = merge_tx.send(temp_db_path);
            } else {
                failed_count.fetch_add(1, Ordering::SeqCst);
                return;
            }

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
}

fn spawn_merge_thread(
    merge_rx: mpsc::Receiver<PathBuf>,
    merge_db_path: PathBuf,
    merge_contigs: Vec<ContigInfo>,
    merge_temp_dir: PathBuf,
    merge_pb: Arc<ProgressBar>,
) -> thread::JoinHandle<usize> {
    thread::spawn(move || {
        let mut merged = 0usize;

        for temp_db_path in merge_rx {
            let sample_name = temp_db_path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();

            if merge_temp_db_into_main(&merge_db_path, &temp_db_path, &merge_contigs).is_err() {
                merge_pb.set_message(format!("ERR: {}", sample_name));
            } else {
                merge_pb.set_message(sample_name);
            }

            let _ = fs::remove_file(&temp_db_path);
            merged += 1;
            merge_pb.inc(1);
        }

        let _ = fs::remove_dir(&merge_temp_dir);
        merge_pb.finish_with_message("Done");
        merged
    })
}

fn print_summary(result: &ProcessResult, output_db: &Path) {
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!(
        "  Samples processed: {}/{}",
        result.samples_processed,
        result.samples_processed + result.samples_failed
    );
    eprintln!("  Total time:        {:.2}s", result.total_time_secs);
    eprintln!();
    eprintln!("  Output: {:?}", output_db);
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
}
