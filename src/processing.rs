//! Parallel BAM file processing and feature calculation.
//!
//! Orchestrates parallel processing of BAM files:
//! 1. Each BAM file processed by separate thread (rayon)
//! 2. Results collected in memory
//! 3. Written to DuckDB sequentially (DuckDB is single-writer)
//!
//! BAM processing (95% of time) is fully parallelized.
//! Database writing (5% of time) runs sequentially after processing completes.

use anyhow::{Context, Result};
use atty::Stream;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::sync_channel;
use std::thread;

use crate::bam_reader::{detect_sequencing_type, process_contig_streaming};
use crate::compress::{
    compress_signal_with_reference, Run,
    add_compressed_feature, add_compressed_feature_with_reference,
    add_compressed_feature_with_stats, merge_identical_runs, 
};
use crate::db::{DbWriter, CompletenessData, DuplicationData};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::genbank::parse_genbank;
use crate::types::{
    ContigInfo, FeaturePoint, PackagingData, PlotType, PresenceData, SequencingType
};

use crate::processing_phage_packaging::{
    apply_dtr_merging, classify_packaging, filter_peaks_by_clippings,
    is_valid_terminal_duplication, merge_peaks, DtrRegion, PhageTerminiConfig,
};
use crate::processing_completeness::compute_completeness;

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
    /// Sequencing type: if None, auto-detect per sample; if Some, use for all samples
    pub sequencing_type: Option<SequencingType>,
    /// Phage termini detection configuration
    pub phagetermini_config: PhageTerminiConfig,
}

impl ProcessConfig {
    /// Parse sequencing type string into SequencingType enum.
    /// If empty or invalid, falls back to auto-detection from BAM file.
    pub fn parse_sequencing_type(seq_type_str: &str) -> Option<SequencingType> {
        match seq_type_str.to_lowercase().as_str() {
            "long" => Some(SequencingType::Long),
            "paired-short" => Some(SequencingType::ShortPaired),
            "single-short" => Some(SequencingType::ShortSingle),
            _ => None,
        }
    }
}

/// Result of processing all samples.
pub struct ProcessResult {
    pub samples_processed: usize,
    pub samples_failed: usize,
    pub total_time_secs: f64,
    pub processing_time_secs: f64,
    pub writing_time_secs: f64,
}
/// Add features from FeatureArrays to output (optimized path).
/// Returns a tuple of:
/// - Optional packaging data for phagetermini module
/// - Optional completeness data for assemblycheck module
fn add_features_from_arrays(
    arrays: &mut FeatureArrays,
    contig_name: &str,
    contig_length: usize,
    config: &ProcessConfig,
    seq_type: SequencingType,
    flags: ModuleFlags,
    duplications: &[DuplicationData],
    output: &mut Vec<FeaturePoint>,
) -> (Option<PackagingData>, Option<CompletenessData>) {
    let pt_config = config.phagetermini_config;

    // Build DTR regions from duplications for this contig (used for both merging and classification)
    // Filter criteria:
    // - Same contig
    // - ≥90% identity
    // - Valid terminal duplication (one region at start, other at end of contig)
    let dtr_regions: Vec<DtrRegion> = if flags.phagetermini && !duplications.is_empty() {
        duplications
            .iter()
            .filter(|d| {
                d.contig_name == contig_name
                    && d.pident >= 90.0
                    && is_valid_terminal_duplication(d, contig_length, pt_config.max_distance_duplication)
            })
            .map(|d| {
                // Determine if direct (DTR) or inverted (ITR)
                let is_direct = (d.position1 < d.position2 && d.position1prime < d.position2prime)
                    || (d.position1 > d.position2 && d.position1prime > d.position2prime);

                // Determine first and second regions (first = lower start position)
                if d.position1 < d.position1prime {
                    DtrRegion {
                        first_start: d.position1.min(d.position2),
                        first_end: d.position1.max(d.position2),
                        second_start: d.position1prime.min(d.position2prime),
                        second_end: d.position1prime.max(d.position2prime),
                        is_direct,
                    }
                } else {
                    DtrRegion {
                        first_start: d.position1prime.min(d.position2prime),
                        first_end: d.position1prime.max(d.position2prime),
                        second_start: d.position1.min(d.position2),
                        second_end: d.position1.max(d.position2),
                        is_direct,
                    }
                }
            })
            .collect()
    } else {
        Vec::new()
    };

    // Coverage (always compress self-referentially)
    let primary_reads_f64: Vec<f64> = arrays.primary_reads.iter().map(|&x| x as f64).collect();
    if flags.coverage {
        // Save strand-specific coverage (VIEW will compute total primary_reads)
        let primary_reads_plus_f64: Vec<f64> = arrays.primary_reads_plus_only.iter().map(|&x| x as f64).collect();
        let primary_reads_minus_f64: Vec<f64> = arrays.primary_reads_minus_only.iter().map(|&x| x as f64).collect();
        add_compressed_feature(&primary_reads_plus_f64, "primary_reads_plus_only", contig_name, config, output);
        add_compressed_feature(&primary_reads_minus_f64, "primary_reads_minus_only", contig_name, config, output);
        
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

        // MAPQ - average mapping quality per position (sum_mapq / primary_reads)
        let mapq_f64: Vec<f64> = arrays.sum_mapq.iter()
            .zip(&arrays.primary_reads)
            .map(|(&sum, &count)| if count > 0 { sum as f64 / count as f64 } else { 0.0 })
            .collect();
        add_compressed_feature(&mapq_f64, "mapq", contig_name, config, output);
    }

    // Assemblycheck features
    let mut left_clip_runs: Vec<Run> = Vec::new();
    let mut right_clip_runs: Vec<Run> = Vec::new();
    let mut insertion_runs: Vec<Run> = Vec::new();
    let mut deletion_runs: Vec<Run> = Vec::new();
    let mut mismatch_runs: Vec<Run> = Vec::new();
    if flags.assemblycheck || flags.phagetermini {
        // Clippings and insertions with statistics
        let left_clip_counts: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| v.len() as f64).collect();
        let left_clip_means: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
        }).collect();
        let left_clip_medians: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();
        let left_clip_stds: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
            if v.len() <= 1 { 0.0 } else {
                let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                var.sqrt()
            }
        }).collect();
        left_clip_runs = add_compressed_feature_with_stats(&left_clip_counts, &left_clip_means, &left_clip_medians, &left_clip_stds,
            Some(&primary_reads_f64), "left_clippings", contig_name, config, output);

        let right_clip_counts: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| v.len() as f64).collect();
        let right_clip_means: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
        }).collect();
        let right_clip_medians: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();
        let right_clip_stds: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
            if v.len() <= 1 { 0.0 } else {
                let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                var.sqrt()
            }
        }).collect();
        right_clip_runs = add_compressed_feature_with_stats(&right_clip_counts, &right_clip_means, &right_clip_medians, &right_clip_stds,
            Some(&primary_reads_f64), "right_clippings", contig_name, config, output);
    }
            
    if flags.assemblycheck {
        let insertion_counts: Vec<f64> = arrays.insertion_lengths.iter().map(|v| v.len() as f64).collect();
        let insertion_means: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
        }).collect();
        let insertion_medians: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();
        let insertion_stds: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.len() <= 1 { 0.0 } else {
                let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                var.sqrt()
            }
        }).collect();
        insertion_runs = add_compressed_feature_with_stats(&insertion_counts, &insertion_means, &insertion_medians, &insertion_stds,
            Some(&primary_reads_f64), "insertions", contig_name, config, output);

        // Other assembly check features (no statistics)
        // Apply two-stage compression: first with coverage reference, then merge identical runs
        let deletions_f64: Vec<f64> = arrays.deletions.iter().map(|&x| x as f64).collect();
        deletion_runs = compress_signal_with_reference(&deletions_f64, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        let deletions_runs_merged = merge_identical_runs(deletion_runs.clone());
        output.extend(deletions_runs_merged.into_iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "deletions".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
            mean: None,
            median: None,
            std: None,
        }));

        let mismatches_f64: Vec<f64> = arrays.mismatches.iter().map(|&x| x as f64).collect();
        mismatch_runs = add_compressed_feature_with_reference(&mismatches_f64, Some(&primary_reads_f64), "mismatches", contig_name, config, output);
        
        if seq_type.is_short_paired() {
            let non_inward_f64: Vec<f64> = arrays.non_inward_pairs.iter().map(|&x| x as f64).collect();
            let non_inward_runs = compress_signal_with_reference(&non_inward_f64, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
            let non_inward_runs_merged = merge_identical_runs(non_inward_runs);
            output.extend(non_inward_runs_merged.into_iter().map(|run| FeaturePoint {
                contig_name: contig_name.to_string(),
                feature: "non_inward_pairs".to_string(),
                start_pos: run.start_pos,
                end_pos: run.end_pos,
                value: run.value,
                mean: None,
                median: None,
                std: None,
            }));
            
            let mate_unmapped_f64: Vec<f64> = arrays.mate_not_mapped.iter().map(|&x| x as f64).collect();
            let mate_unmapped_runs = compress_signal_with_reference(&mate_unmapped_f64, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
            let mate_unmapped_runs_merged = merge_identical_runs(mate_unmapped_runs);
            output.extend(mate_unmapped_runs_merged.into_iter().map(|run| FeaturePoint {
                contig_name: contig_name.to_string(),
                feature: "mate_not_mapped".to_string(),
                start_pos: run.start_pos,
                end_pos: run.end_pos,
                value: run.value,
                mean: None,
                median: None,
                std: None,
            }));
            
            let mate_other_contig_f64: Vec<f64> = arrays.mate_on_another_contig.iter().map(|&x| x as f64).collect();
            let mate_other_contig_runs = compress_signal_with_reference(&mate_other_contig_f64, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
            let mate_other_contig_runs_merged = merge_identical_runs(mate_other_contig_runs);
            output.extend(mate_other_contig_runs_merged.into_iter().map(|run| FeaturePoint {
                contig_name: contig_name.to_string(),
                feature: "mate_on_another_contig".to_string(),
                start_pos: run.start_pos,
                end_pos: run.end_pos,
                value: run.value,
                mean: None,
                median: None,
                std: None,
            }));
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

    // Phagetermini features
    let packaging_result = if flags.phagetermini {
        // === STEP 1: Save ORIGINAL data to database (before any DTR merging) ===

        // coverage_reduced (self-referential curve)
        let coverage_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
        add_compressed_feature(&coverage_reduced_f64, "coverage_reduced", contig_name, config, output);

        // reads_starts and reads_ends (bars) - save original data
        let reads_starts_original: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_starts_runs = compress_signal_with_reference(&reads_starts_original, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);

        let reads_ends_original: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let reads_ends_runs = compress_signal_with_reference(&reads_ends_original, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);

        // Add reads_starts and reads_ends to output (original data)
        output.extend(reads_starts_runs.iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "reads_starts".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
            mean: None,
            median: None,
            std: None,
        }));

        output.extend(reads_ends_runs.iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "reads_ends".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
            mean: None,
            median: None,
            std: None,
        }));

        // Tau: calculate tau values using original data
        for run in reads_starts_runs.iter().chain(reads_ends_runs.iter()) {
            let mut tau_sum = 0.0;
            let mut count = 0;

            for pos in run.start_pos..=run.end_pos {
                let idx = (pos - 1) as usize;
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
                    mean: None,
                    median: None,
                    std: None,
                });
            }
        }

        // === STEP 2: Apply DTR merging to arrays ===
        if !dtr_regions.is_empty() {
            apply_dtr_merging(arrays, &dtr_regions);
        }

        // === STEP 3: Recompute runs from MERGED data ===
        let coverage_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
        let reads_starts_merged: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_ends_merged: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let reads_starts_runs_merged = compress_signal_with_reference(&reads_starts_merged, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        let reads_ends_runs_merged = compress_signal_with_reference(&reads_ends_merged, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        
        let left_clip_counts: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| v.len() as f64).collect();
        let right_clip_counts: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| v.len() as f64).collect();
        let left_clip_runs = compress_signal_with_reference(&left_clip_counts, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        let right_clip_runs = compress_signal_with_reference(&right_clip_counts, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);

        // === STEP 4: Filter and merge peaks ===
        let (reads_starts_filtered, reads_ends_filtered) = filter_peaks_by_clippings(
            &reads_starts_runs_merged,
            &reads_ends_runs_merged,
            &left_clip_runs,
            &right_clip_runs,
            pt_config.min_events,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular,
            &dtr_regions,
        );

        let start_peaks = merge_peaks(
            &reads_starts_filtered,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular,
            &dtr_regions,
        );
        let end_peaks = merge_peaks(
            &reads_ends_filtered,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular,
            &dtr_regions,
        );

        let (mechanism, left_termini, right_termini, duplication) = classify_packaging(
            &start_peaks,
            &end_peaks,
            contig_length,
            config.circular,
            &dtr_regions,
        );

        // Only return PackagingData if there's a detected mechanism (not "No_packaging")
        if mechanism != "No_packaging" {
            Some(PackagingData {
                contig_name: contig_name.to_string(),
                mechanism,
                left_termini,
                right_termini,
                duplication,
            })
        } else {
            None
        }
    } else {
        None
    };

    // Compute completeness statistics when assemblycheck is enabled
    let completeness_result = if flags.assemblycheck {
        let completeness = compute_completeness(
            &arrays.left_clipping_lengths,
            &arrays.right_clipping_lengths,
            &arrays.insertion_lengths,
            &arrays.deletion_lengths,
            &arrays.primary_reads,
            &left_clip_runs,
            &right_clip_runs,
            &insertion_runs,
            &deletion_runs,
            &mismatch_runs,
            contig_name,
            contig_length,
        );
        if completeness.has_data() {
            Some(completeness)
        } else {
            None
        }
    } else {
        None
    };

    (packaging_result, completeness_result)
}

/// Process one BAM file using streaming (single-pass, optimized).
/// Parallelizes at the contig level for samples with many contigs.
pub fn process_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    duplications: &[DuplicationData],
) -> Result<(Vec<FeaturePoint>, Vec<PresenceData>, Vec<PackagingData>, Vec<CompletenessData>, String)> {
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "");

    let flags = ModuleFlags::from_modules(modules);

    // Get number of reference sequences using temporary reader
    let temp_bam = IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;
    let n_refs = temp_bam.header().target_count();
    drop(temp_bam);

    // Determine sequencing type: use provided value or auto-detect from this BAM file
    let seq_type = match config.sequencing_type {
        Some(st) => st,
        None => detect_sequencing_type(bam_path)
            .with_context(|| format!("Failed to detect sequencing type from {}", bam_path.display()))?,
    };

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
            let (mut arrays, coverage_pct) = match process_contig_streaming(&mut bam, &ref_name, ref_length, seq_type, flags, config.circular, config.min_coverage) {
                Ok(Some(result)) => result,
                Ok(None) => return None,
                Err(e) => {
                    eprintln!("Error processing contig {} in {}: {}", ref_name, bam_path.display(), e);
                    return None;
                }
            };

            // Coverage already checked and returned from process_contig_streaming

            // Calculate features for this contig
            let mut features = Vec::new();
            let (packaging_info, completeness_info) = add_features_from_arrays(&mut arrays, &ref_name, ref_length, config, seq_type, flags, duplications, &mut features);

            // Calculate Coverage-Weighted Total Variation (CWTV)
            // CWTV = 1/(n-1) * Σ|cov(i+1) - cov(i)| for i=1 to n (circular wrap-around)
            let coverage_variation = if arrays.primary_reads.len() > 1 {
                let n = arrays.primary_reads.len();
                // Sum of consecutive differences (n-1 terms)
                let total_variation: f64 = arrays.primary_reads
                    .windows(2)
                    .map(|w| (w[1] as f64 - w[0] as f64).abs())
                    .sum();
                (total_variation / (n - 1) as f64) as f32
            } else {
                0.0
            };

            let presence = PresenceData {
                contig_name: ref_name.clone(),
                coverage_pct: coverage_pct as f32,
                coverage_variation,
            };

            Some((features, presence, packaging_info, completeness_info))
        })
        .collect();

    // Merge results from all contigs
    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();
    let mut all_packaging = Vec::new();
    let mut all_completeness = Vec::new();

    for (features, presence, packaging, completeness) in results {
        all_features.extend(features);
        all_presences.push(presence);
        if let Some(pkg) = packaging {
            all_packaging.push(pkg);
        }
        if let Some(comp) = completeness {
            all_completeness.push(comp);
        }
    }

    Ok((all_features, all_presences, all_packaging, all_completeness, sample_name))
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
    _create_indexes: bool, // Ignored - DuckDB uses zone maps instead of indexes
    max_samples_in_memory: usize,
    autoblast_file: &Path,
) -> Result<ProcessResult> {
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    eprintln!("\n### Parsing input files...");
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

    eprintln!("Found {} BAM files", bam_files.len());

    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent).context("Failed to create output directory")?;
    }

    std::thread::sleep(std::time::Duration::from_millis(50));

    // Create database with schema and initial data
    let db_writer = DbWriter::create(output_db, &contigs, &annotations)?;

    // Parse and write duplications from autoblast file (if provided)
    let duplications = if !autoblast_file.as_os_str().is_empty() && autoblast_file.exists() {
        eprintln!("\n### Parsing autoblast results...");
        let dups = crate::db::parse_autoblast_file(autoblast_file)?;
        if !dups.is_empty() {
            db_writer.write_duplications(&dups)?;
        } else {
            eprintln!("No duplications found in autoblast file");
        }
        dups
    } else {
        Vec::new()
    };

    eprintln!("\n### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!("Modules: {}\n", modules.join(", "));

    let result = process_samples_parallel(&bam_files, &contigs, modules, config, db_writer, max_samples_in_memory, &duplications)?;

    print_summary(&result, output_db);

    Ok(result)
}

/// Holds processed sample data ready for database writing.
struct SampleResult {
    sample_name: String,
    features: Vec<FeaturePoint>,
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    completeness: Vec<CompletenessData>,
}

fn process_samples_parallel(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    db_writer: DbWriter,
    max_samples_in_memory: usize,
    duplications: &[DuplicationData],
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let is_tty = atty::is(Stream::Stderr);
    let start_time = std::time::Instant::now();

    // Sequential mode for single thread: process one → write one → repeat
    if config.threads == 1 {
        return process_samples_sequential(bam_files, contigs, modules, config, db_writer, is_tty, duplications);
    }

    // Parallel mode: producer-consumer with bounded channel
    let mp = MultiProgress::new();

    // For TTY: use interactive progress bars
    // For non-TTY (SLURM logs): use hidden bars and print log messages
    let process_pb = if is_tty {
        let pb = mp.add(ProgressBar::new(total as u64));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        mp.add(ProgressBar::hidden())
    };

    let write_pb = if is_tty {
        let pb = mp.add(ProgressBar::new(total as u64));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.yellow} Writing:     [{elapsed_precise}] [{bar:40.yellow/red}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        mp.add(ProgressBar::hidden())
    };

    // Bounded channel limits memory to ~max_samples_in_memory samples
    let (tx, rx) = sync_channel::<SampleResult>(max_samples_in_memory);

    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);

    // Clone values needed by writer thread
    let write_pb_clone = write_pb.clone();
    let is_tty_writer = is_tty;
    let total_writer = total;

    // Spawn dedicated writer thread
    let writer_handle = thread::spawn(move || -> Result<(usize, std::time::Duration)> {
        let write_start = std::time::Instant::now();
        let mut written_count = 0usize;

        // Receive and write samples as they arrive
        for result in rx {
            written_count += 1;

            // Insert sample
            if let Err(e) = db_writer.insert_sample(&result.sample_name) {
                eprintln!("\nError inserting sample {}: {}", result.sample_name, e);
                let msg = format!("ERR: {}", result.sample_name);
                write_pb_clone.set_message(msg.clone());
                write_pb_clone.inc(1);
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                }
                continue;
            }

            // Write all data for this sample
            if let Err(e) = db_writer.write_sample_data(
                &result.sample_name,
                &result.presences,
                &result.packaging,
                &result.completeness,
                &result.features,
            ) {
                eprintln!("\nError writing data for {}: {}", result.sample_name, e);
                let msg = format!("ERR: {}", result.sample_name);
                write_pb_clone.set_message(msg.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                }
            } else {
                write_pb_clone.set_message(result.sample_name.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, result.sample_name);
                }
            }
            write_pb_clone.inc(1);
        }

        // Finalize database
        db_writer.finalize()?;
        write_pb_clone.finish_with_message("Done");
        if !is_tty_writer {
            eprintln!("Writing:    Done ({:.2}s)", write_start.elapsed().as_secs_f64());
        }

        Ok((written_count, write_start.elapsed()))
    });

    // Process samples in parallel, sending to channel immediately
    // send() blocks if channel is full (backpressure)
    bam_files.par_iter().for_each(|bam_path| {
        let sample_start = std::time::Instant::now();

        match process_sample(bam_path, contigs, modules, config, duplications) {
            Ok((features, presences, packaging, completeness, sample_name)) => {
                let sample_time = sample_start.elapsed().as_secs_f64();
                completed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("{} ({:.2}s)", sample_name, sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("Processing: {}", msg);
                }

                // Send to writer thread (blocks if channel full)
                let _ = tx.send(SampleResult {
                    sample_name,
                    features,
                    presences,
                    packaging,
                    completeness,
                });
            }
            Err(e) => {
                eprintln!("\nError processing {}: {}", bam_path.display(), e);
                let sample_time = sample_start.elapsed().as_secs_f64();
                completed_count.fetch_add(1, Ordering::SeqCst);
                failed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("ERR ({:.2}s)", sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("Processing: {}", msg);
                }
            }
        }
    });

    // Drop sender to signal writer thread that processing is complete
    drop(tx);

    process_pb.finish_with_message("Done");
    let processing_time = start_time.elapsed();
    if !is_tty {
        eprintln!("Processing: Done ({:.2}s)", processing_time.as_secs_f64());
    }

    // Wait for writer thread to finish
    let (written_count, writing_time) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);

    Ok(ProcessResult {
        samples_processed: written_count,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
        processing_time_secs: processing_time.as_secs_f64(),
        writing_time_secs: writing_time.as_secs_f64(),
    })
}

/// Sequential processing for single-threaded mode.
/// Process one sample → write it → repeat. No channel overhead.
fn process_samples_sequential(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    db_writer: DbWriter,
    is_tty: bool,
    duplications: &[DuplicationData],
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let start_time = std::time::Instant::now();

    let pb = if is_tty {
        let pb = ProgressBar::new(total as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        ProgressBar::hidden()
    };

    let mut processed = 0usize;
    let mut failed = 0usize;
    let mut processing_time_total = std::time::Duration::ZERO;
    let mut writing_time_total = std::time::Duration::ZERO;

    for bam_path in bam_files {
        let sample_start = std::time::Instant::now();

        match process_sample(bam_path, contigs, modules, config, duplications) {
            Ok((features, presences, packaging, completeness, sample_name)) => {
                let process_elapsed = sample_start.elapsed();
                processing_time_total += process_elapsed;

                let write_start = std::time::Instant::now();

                // Insert sample
                if let Err(e) = db_writer.insert_sample(&sample_name) {
                    eprintln!("\nError inserting sample {}: {}", sample_name, e);
                    failed += 1;
                    pb.inc(1);
                    continue;
                }

                // Write sample data
                if let Err(e) = db_writer.write_sample_data(
                    &sample_name,
                    &presences,
                    &packaging,
                    &completeness,
                    &features,
                ) {
                    eprintln!("\nError writing data for {}: {}", sample_name, e);
                    failed += 1;
                } else {
                    processed += 1;
                }

                writing_time_total += write_start.elapsed();

                let total_sample_time = sample_start.elapsed().as_secs_f64();
                let msg = format!("{} ({:.2}s)", sample_name, total_sample_time);
                pb.set_message(msg.clone());
                pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
            Err(e) => {
                eprintln!("\nError processing {}: {}", bam_path.display(), e);
                failed += 1;
                processing_time_total += sample_start.elapsed();
                let msg = "ERR".to_string();
                pb.set_message(msg.clone());
                pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
        }
    }

    // Finalize database
    db_writer.finalize()?;
    pb.finish_with_message("Done");
    if !is_tty {
        eprintln!("Done ({:.2}s)", start_time.elapsed().as_secs_f64());
    }

    Ok(ProcessResult {
        samples_processed: processed,
        samples_failed: failed,
        total_time_secs: start_time.elapsed().as_secs_f64(),
        processing_time_secs: processing_time_total.as_secs_f64(),
        writing_time_secs: writing_time_total.as_secs_f64(),
    })
}

fn print_summary(result: &ProcessResult, output_db: &Path) {
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!("  Samples processed: {}/{}", result.samples_processed, result.samples_processed + result.samples_failed);
    eprintln!("  Total time:        {:.2}s", result.total_time_secs);
    eprintln!();
    eprintln!("  Output: {:?}", output_db);
}
