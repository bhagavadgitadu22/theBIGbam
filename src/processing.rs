//! Parallel BAM file processing and feature calculation.
//!
//! Orchestrates parallel processing of BAM files:
//! 1. Each BAM file processed by separate thread (rayon)
//! 2. Features written to temporary SQLite database per sample
//! 3. Temp databases merged sequentially into main database
//!
//! BAM processing (95% of time) is fully parallelized.
//! Database merging (5% of time) runs sequentially after processing completes.

use anyhow::{Context, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::sync::Mutex;

use crate::bam_reader::{detect_sequencing_type, process_contig_streaming};
use crate::compress::compress_signal_with_reference;
use crate::db::{
    create_metadata_db, create_temp_sample_db, finalize_db, merge_temp_db_into_main,
    write_features_to_temp_db, write_presences_to_temp_db,
};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::genbank::parse_genbank;
use crate::types::{
    get_plot_type, ContigInfo, FeaturePoint, PlotType, PresenceData, ASSEMBLYCHECK_FEATURES
};

/// Configuration for processing.
#[derive(Clone)]
pub struct ProcessConfig {
    pub threads: usize,
    pub min_coverage: f64,
    /// Relative tolerance for RLE compression (e.g., 0.1 = 10% change threshold)
    pub curve_ratio: f64,
    pub bar_ratio: f64,
    /// Circular genome flag: if true, assembly was doubled during mapping (enables modulo logic)
    pub circular: bool,
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
    add_compressed_feature_with_reference(values, None, feature, contig_name, config, output);
}

/// Compress and add feature points with optional coverage reference.
///
/// For bar plots, uses coverage as reference for context-aware compression.
#[inline]
fn add_compressed_feature_with_reference(
    values: &[f64],
    reference: Option<&[f64]>,
    feature: &str,
    contig_name: &str,
    config: &ProcessConfig,
    output: &mut Vec<FeaturePoint>,
) {
    let plot_type = get_plot_type(feature);
    let runs = compress_signal_with_reference(values, reference, plot_type, config.curve_ratio, config.bar_ratio);

    output.extend(runs.into_iter().map(|run| FeaturePoint {
        contig_name: contig_name.to_string(),
        feature: feature.to_string(),
        start_pos: run.start_pos,
        end_pos: run.end_pos,
        value: run.value,
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
    // Coverage (always compress self-referentially)
    let primary_reads_f64: Vec<f64> = arrays.primary_reads.iter().map(|&x| x as f64).collect();
    if flags.coverage {
        add_compressed_feature(&primary_reads_f64, "primary_reads", contig_name, config, output);
        
        // Secondary reads (self-referential curve)
        // When circular=true, subtract primary coverage to remove artifact secondary alignments from doubled assembly
        let secondary_reads_f64: Vec<f64> = if config.circular {
            arrays.secondary_reads.iter()
                .zip(&arrays.primary_reads)
                .map(|(&sec, &cov)| if sec > cov { (sec - cov) as f64 } else { 0.0 })
                .collect()
        } else {
            arrays.secondary_reads.iter().map(|&x| x as f64).collect()
        };
        add_compressed_feature(&secondary_reads_f64, "secondary_reads", contig_name, config, output);
        
        // Supplementary reads (self-referential curve)
        let supplementary_reads_f64: Vec<f64> = arrays.supplementary_reads.iter().map(|&x| x as f64).collect();
        add_compressed_feature(&supplementary_reads_f64, "supplementary_reads", contig_name, config, output);
    }

    // Coverage_reduced for phagetermini (use as reference for phagetermini bar features)
    let coverage_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
    
    // Phagetermini features
    if flags.phagetermini {
        // coverage_reduced (self-referential curve)
        add_compressed_feature(&coverage_reduced_f64, "coverage_reduced", contig_name, config, output);
        
        // reads_starts and reads_ends (bars, use coverage_reduced as reference)
        let reads_starts: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_starts_runs = compress_signal_with_reference(&reads_starts, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        
        let reads_ends: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let reads_ends_runs = compress_signal_with_reference(&reads_ends, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);

        // Add reads_starts and reads_ends to output
        output.extend(reads_starts_runs.iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "reads_starts".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
        }));
        
        output.extend(reads_ends_runs.iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "reads_ends".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
        }));

        // Tau: calculate tau values directly for positions where reads_starts or reads_ends were saved
        // Merge the runs from both features and calculate tau for those positions
        for run in reads_starts_runs.iter().chain(reads_ends_runs.iter()) {
            // Calculate average tau value for this run
            let mut tau_sum = 0.0;
            let mut count = 0;
            
            for pos in run.start_pos..=run.end_pos {
                let idx = (pos - 1) as usize; // Convert from 1-indexed to 0-indexed
                if idx < arrays.coverage_reduced.len() {
                    let c = arrays.coverage_reduced[idx];
                    let s = arrays.reads_starts[idx];
                    let e = arrays.reads_ends[idx];
                    if c > 0 {
                        tau_sum += (s + e) as f64 / c as f64;
                        count += 1;
                    }
                }
            }
            
            if count > 0 {
                let tau_value = tau_sum / count as f64;
                output.push(FeaturePoint {
                    contig_name: contig_name.to_string(),
                    feature: "tau".to_string(),
                    start_pos: run.start_pos,
                    end_pos: run.end_pos,
                    value: tau_value as f32,
                });
            }
        }
    }

    // Assemblycheck features (all bars, use main coverage as reference)
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
            add_compressed_feature_with_reference(&values, Some(&primary_reads_f64), feature_name, contig_name, config, output);
        }

        // Read lengths (curve for long reads, self-referential)
        if seq_type.is_long() {
            let values: Vec<f64> = arrays
                .sum_read_lengths
                .iter()
                .zip(&arrays.count_read_lengths)
                .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                .collect();
            add_compressed_feature(&values, "read_lengths", contig_name, config, output);
        }

        // Insert sizes (curve for paired reads, self-referential)  
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
/// Parallelizes at the contig level for samples with many contigs.
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

    // Get number of reference sequences using temporary reader
    let temp_bam = IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;
    let n_refs = temp_bam.header().target_count();
    drop(temp_bam);

    // Process contigs in parallel - each gets its own BAM reader
    // This is critical for samples with many contigs (e.g., 50,000 contigs with 1 sample)
    let results: Vec<_> = (0..n_refs)
        .into_par_iter()
        .filter_map(|tid| {
            // Each thread gets its own BAM reader
            let mut bam = IndexedReader::from_path(bam_path).ok()?;
            // Use 1 decompression thread per reader to avoid I/O contention
            // Total decompression threads = config.threads (one per parallel contig)
            bam.set_threads(1).ok()?;

            // Extract header info before mutable borrow
            let ref_name = std::str::from_utf8(bam.header().tid2name(tid)).ok()?.to_string();
            let bam_length = bam.header().target_len(tid).unwrap_or(0) as usize;
            let ref_length = if config.circular { bam_length / 2 } else { bam_length };

            // Skip if contig not in GenBank list
            if contigs.iter().find(|c| c.name == ref_name).is_none() {
                return None;
            }

            // Process contig using streaming - single pass over reads
            // Early coverage check happens inside process_contig_streaming
            let (arrays, coverage_pct) = match process_contig_streaming(&mut bam, &ref_name, ref_length, seq_type, flags, config.circular, config.min_coverage) {
                Ok(Some(result)) => result,
                Ok(None) => return None,
                Err(e) => {
                    eprintln!("Error processing contig {} in {}: {}", ref_name, bam_path.display(), e);
                    return None;
                }
            };

            // Coverage already checked and returned from process_contig_streaming
            let presence = PresenceData {
                contig_name: ref_name.clone(),
                coverage_pct: coverage_pct as f32,
            };

            // Calculate features for this contig
            let mut features = Vec::new();
            add_features_from_arrays(&arrays, &ref_name, config, flags, seq_type, &mut features);

            Some((features, presence))
        })
        .collect();

    // Merge results from all contigs
    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();
    
    for (features, presence) in results {
        all_features.extend(features);
        all_presences.push(presence);
    }

    Ok((all_features, all_presences, sample_name))
}

/// Extract contig information from BAM file headers.
/// Deduplicates contigs across all BAM files.
fn extract_contigs_from_bams(bam_files: &[PathBuf], circular: bool) -> Result<Vec<ContigInfo>> {
    use std::collections::HashMap;
    
    let mut contig_map: HashMap<String, usize> = HashMap::new();
    
    // Scan all BAM files to collect unique contigs
    for bam_path in bam_files {
        let bam = IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path.display()))?;
        let header = bam.header();
        
        for tid in 0..header.target_count() {
            let ref_name = std::str::from_utf8(header.tid2name(tid))
                .context("Invalid UTF-8 in reference name")?
                .to_string();
            let bam_length = header.target_len(tid).unwrap_or(0) as usize;
            
            // For circular genomes, BAM length is doubled, so actual length is half
            let actual_length = if circular { bam_length / 2 } else { bam_length };
            
            // Use the contig if not seen, or verify length matches
            contig_map.entry(ref_name.clone())
                .and_modify(|existing_len| {
                    if *existing_len != actual_length {
                        eprintln!("Warning: Contig '{}' has different lengths across BAM files ({} vs {})", 
                                  ref_name, *existing_len, actual_length);
                    }
                })
                .or_insert(actual_length);
        }
    }
    
    // Convert to sorted vector of ContigInfo
    let mut contigs: Vec<ContigInfo> = contig_map
        .into_iter()
        .map(|(name, length)| ContigInfo {
            name,
            length,
            annotation_tool: String::new(),
        })
        .collect();
    
    // Sort by name for consistent ordering
    contigs.sort_by(|a, b| a.name.cmp(&b.name));
    
    Ok(contigs)
}

/// Run processing on all BAM files.
pub fn run_all_samples(
    genbank_path: &Path,
    bam_files: &[PathBuf],
    output_db: &Path,
    modules: &[String],
    annotation_tool: &str,
    config: &ProcessConfig,
    create_indexes: bool,
) -> Result<ProcessResult> {
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    eprintln!("### Parsing GenBank file...");
    let (contigs, annotations) = if genbank_path.as_os_str().is_empty() {
        eprintln!("No GenBank file provided - extracting contigs from BAM headers");
        let contigs = extract_contigs_from_bams(bam_files, config.circular)?;
        eprintln!("Found {} contigs from BAM files", contigs.len());
        (contigs, Vec::new())
    } else {
        parse_genbank(genbank_path, annotation_tool)?
    };
    eprintln!(
        "Found {} contigs with {} annotations",
        contigs.len(),
        annotations.len()
    );

    eprintln!("Found {} BAM files\n", bam_files.len());

    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent).context("Failed to create output directory")?;
    }
    let temp_dir = output_db.with_extension("temp_dbs");
    fs::create_dir_all(&temp_dir).context("Failed to create temp directory")?;

    std::thread::sleep(std::time::Duration::from_millis(50));

    create_metadata_db(output_db, &contigs, &annotations, create_indexes)?;

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

    let start_time = std::time::Instant::now();
    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);

    // Collect temp DB paths from parallel processing
    let temp_db_paths = Arc::new(Mutex::new(Vec::new()));

    bam_files.par_iter().for_each(|bam_path| {
        process_single_sample(
            bam_path,
            contigs,
            modules,
            config,
            temp_dir,
            Arc::clone(&temp_db_paths),
            &completed_count,
            &failed_count,
            &process_pb,
            total,
        );
    });

    process_pb.finish_with_message("Done");

    // Merge all temp DBs sequentially (fast operation, not a bottleneck)
    let temp_paths = temp_db_paths.lock().unwrap();
    for temp_db_path in temp_paths.iter() {
        let sample_name = temp_db_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        if merge_temp_db_into_main(db_path, temp_db_path, contigs).is_err() {
            merge_pb.set_message(format!("ERR: {}", sample_name));
        } else {
            merge_pb.set_message(sample_name);
        }

        let _ = fs::remove_file(temp_db_path);
        merge_pb.inc(1);
    }
    drop(temp_paths);

    let _ = fs::remove_dir(temp_dir);
    merge_pb.finish_with_message("Done");

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
    temp_db_paths: Arc<Mutex<Vec<PathBuf>>>,
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

            let temp_conn = match create_temp_sample_db(&temp_db_path) {
                Ok(conn) => conn,
                Err(e) => {
                    eprintln!("\nError creating temp DB for {}: {}", sample_name, e);
                    failed_count.fetch_add(1, Ordering::SeqCst);
                    return;
                }
            };

            if !features.is_empty() {
                if let Err(e) = write_features_to_temp_db(&temp_conn, &features) {
                    eprintln!("\nError writing features for {}: {}", sample_name, e);
                    failed_count.fetch_add(1, Ordering::SeqCst);
                    return;
                }
            }

            if !presences.is_empty() {
                if let Err(e) = write_presences_to_temp_db(&temp_conn, &sample_name, &presences) {
                    eprintln!("\nError writing presences for {}: {}", sample_name, e);
                    failed_count.fetch_add(1, Ordering::SeqCst);
                    return;
                }
            }

            drop(temp_conn);
            
            // Add temp DB path to collection for later merging
            temp_db_paths.lock().unwrap().push(temp_db_path);

            let sample_time = sample_start.elapsed().as_secs_f64();
            let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
            process_pb.set_message(format!("[{}/{}] {} ({:.2}s)", done, total, sample_name, sample_time));
            process_pb.inc(1);
        }
        Err(e) => {
            eprintln!("\nError processing {}: {}", bam_path.display(), e);
            let sample_time = sample_start.elapsed().as_secs_f64();
            let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
            failed_count.fetch_add(1, Ordering::SeqCst);
            process_pb.set_message(format!("[{}/{}] ERR ({:.2}s)", done, total, sample_time));
            process_pb.inc(1);
        }
    }
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
