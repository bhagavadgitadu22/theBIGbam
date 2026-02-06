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

use crate::bam_reader::{detect_sequencing_type, get_total_read_count, process_contig_streaming};
use crate::compress::{
    compress_signal_with_reference, Run,
    add_compressed_feature, add_compressed_feature_with_reference,
    add_compressed_feature_with_stats, merge_identical_runs, 
};
use crate::db::{DbWriter, CompletenessData, GCContentData, RepeatsData};
use crate::gc_content::{compute_gc_content, compute_gc_skew, GCParams};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::parser::parse_annotations;
use crate::types::{
    ContigInfo, FeaturePoint, PackagingData, PlotType, PresenceData, SequencingType
};

use crate::processing_phage_packaging::{
    classify_packaging_areas, filter_and_merge_to_areas_with_diagnostics,
    is_valid_terminal_repeat, DtrRegion, PhageTerminiConfig,
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
    /// RLE compression tolerance for contig-level features like GC content (default 0.1%)
    pub contig_variation_percentage: f64,
    /// Circular genome flag: if true, assembly was doubled during mapping (enables modulo logic)
    pub circular: bool,
    /// Sequencing type: if None, auto-detect per sample; if Some, use for all samples
    pub sequencing_type: Option<SequencingType>,
    /// Phage termini detection configuration
    pub phagetermini_config: PhageTerminiConfig,
    /// GC content and GC skew parameters
    pub gc_params: GCParams,
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

/// Check if a BAM file is missing MD tags by sampling first few reads.
/// Returns true if MD tags are missing, false otherwise.
fn check_missing_md_tags(bam: &mut IndexedReader) -> bool {
    let header = bam.header().clone();
    let target_names = header.target_names();
    let first_tid = target_names.first();

    if first_tid.is_none() {
        return false;
    }

    let contig_name = match std::str::from_utf8(first_tid.unwrap()) {
        Ok(s) => s,
        Err(_) => return false,
    };

    if bam.fetch(contig_name).is_err() {
        return false;
    }

    let mut checked = 0;
    let mut missing_md = 0;

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped() {
            continue;
        }

        checked += 1;
        if record.aux(b"MD").is_err() {
            missing_md += 1;
        }

        if checked >= 10 {
            break;
        }
    }

    checked > 0 && missing_md == checked
}

/// Check for circularity mismatch between BAM and GenBank/FASTA.
/// Returns Some(warning_message) if mismatch detected, None otherwise.
fn check_circularity_mismatch(
    bam: &IndexedReader,
    contigs: &[ContigInfo],
    circular: bool,
    sample_name: &str,
) -> Option<String> {
    let header = bam.header();

    // Check first contig only to avoid spam
    if header.target_count() == 0 {
        return None;
    }

    let ref_name = match std::str::from_utf8(header.tid2name(0)) {
        Ok(s) => s,
        Err(_) => return None,
    };
    let bam_length = header.target_len(0).unwrap_or(0) as usize;

    // Find matching contig in GenBank/FASTA
    let contig_info = contigs.iter().find(|c| c.name == ref_name)?;
    let genbank_length = contig_info.length;

    if circular {
        // --circular used: expect bam_length = 2 * genbank_length
        if bam_length == genbank_length {
            return Some(format!(
                "Warning: in sample '{}', BAM references are the same size as the FASTA contigs, \
                 you likely did a normal mapping and should not use the --circular flag",
                sample_name
            ));
        }
    } else {
        // --circular NOT used: expect bam_length = genbank_length
        if bam_length == genbank_length * 2 {
            return Some(format!(
                "Warning: in sample '{}', BAM references are twice the size of FASTA contigs, \
                 you likely did a circular mapping and should use the --circular flag",
                sample_name
            ));
        }
    }

    None
}

/// Run validation checks on all BAM files before processing.
/// Prints warnings for MD tags and circularity mismatches.
fn validate_all_samples(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    modules: &[String],
    circular: bool,
) {
    let flags = ModuleFlags::from_modules(modules);
    let needs_md = flags.needs_md();

    for bam_path in bam_files {
        let sample_name = bam_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .replace("_with_MD", "");

        // Open BAM for checks
        let bam = match IndexedReader::from_path(bam_path) {
            Ok(b) => b,
            Err(_) => continue, // Skip, will fail later during processing
        };

        // Circularity check
        if let Some(warning) = check_circularity_mismatch(&bam, contigs, circular, &sample_name) {
            eprintln!("{}", warning);
        }

        // MD tag check
        if needs_md {
            let mut check_bam = match IndexedReader::from_path(bam_path) {
                Ok(b) => b,
                Err(_) => continue,
            };
            if check_missing_md_tags(&mut check_bam) {
                eprintln!(
                    "Warning: in sample '{}', BAM file is missing MD tags. \
                     MD tags are required for: Misalignment, Phage termini. \
                     Consider using 'samtools calmd' to add them.",
                    sample_name
                );
            }
        }
    }
}

/// Add features from FeatureArrays to output (optimized path).
/// Returns a tuple of:
/// - Optional packaging data for phagetermini module (includes all termini with filtering metadata)
/// - Optional completeness data for assemblycheck module
fn add_features_from_arrays(
    arrays: &mut FeatureArrays,
    contig_name: &str,
    contig_length: usize,
    config: &ProcessConfig,
    seq_type: SequencingType,
    flags: ModuleFlags,
    repeats: &[RepeatsData],
    output: &mut Vec<FeaturePoint>,
) -> (Option<PackagingData>, Option<CompletenessData>) {
    let pt_config = config.phagetermini_config;

    // Build DTR regions from repeats for this contig (used for both merging and classification)
    // Filter criteria:
    // - Same contig
    // - ≥90% identity
    // - Valid terminal repeat (one region at start, other at end of contig)
    let dtr_regions: Vec<DtrRegion> = if flags.phagetermini && !repeats.is_empty() {
        repeats
            .iter()
            .filter(|d| {
                d.contig_name == contig_name
                    && d.pident >= pt_config.min_identity_dtr as f64
                    && is_valid_terminal_repeat(d, contig_length, pt_config.max_distance_duplication)
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
    if flags.mapping_metrics || flags.phagetermini {
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
            
    if flags.mapping_metrics {
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

        // Other mapping metrics features (no statistics)
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
    }

    // Paired-reads module
    if flags.paired_read_metrics && seq_type.is_short_paired() {
        let non_inward_f64: Vec<f64> = arrays.non_inward_pairs.iter().map(|&x| x as f64).collect();
        let non_inward_runs = compress_signal_with_reference(&non_inward_f64, Some(&primary_reads_f64), PlotType::Curve, config.curve_ratio, config.bar_ratio);
        output.extend(non_inward_runs.into_iter().map(|run| FeaturePoint {
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
        let mate_unmapped_runs = compress_signal_with_reference(&mate_unmapped_f64, Some(&primary_reads_f64), PlotType::Curve, config.curve_ratio, config.bar_ratio);
        output.extend(mate_unmapped_runs.into_iter().map(|run| FeaturePoint {
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
        let mate_other_contig_runs = compress_signal_with_reference(&mate_other_contig_f64, Some(&primary_reads_f64), PlotType::Curve, config.curve_ratio, config.bar_ratio);
        output.extend(mate_other_contig_runs.into_iter().map(|run| FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: "mate_on_another_contig".to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
            mean: None,
            median: None,
            std: None,
        }));

        // Insert sizes (curve for paired reads, self-referential)
        let values: Vec<f64> = arrays
            .sum_insert_sizes
            .iter()
            .zip(&arrays.count_insert_sizes)
            .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
            .collect();
        add_compressed_feature(&values, "insert_sizes", contig_name, config, output);
    }

    // Long-reads module
    if flags.long_read_metrics && seq_type.is_long() {
        let values: Vec<f64> = arrays
            .sum_read_lengths
            .iter()
            .zip(&arrays.count_read_lengths)
            .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
            .collect();
        add_compressed_feature(&values, "read_lengths", contig_name, config, output);
    }

    // Phagetermini features
    let packaging_result = if flags.phagetermini {
        // === STEP 1: Save ORIGINAL data to database (before any DTR merging) ===
        // coverage_reduced (self-referential curve)
        let coverage_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
        add_compressed_feature(&coverage_reduced_f64, "coverage_reduced", contig_name, config, output);

        // reads_starts and reads_ends (bars) - save original data
        let reads_starts_original: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_starts_runs = compress_signal_with_reference(
            &reads_starts_original, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio
        );

        let reads_ends_original: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let reads_ends_runs = compress_signal_with_reference(
            &reads_ends_original, Some(&coverage_reduced_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio
        );

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

        // === STEP 2: Calculate aligned fraction and global metrics ===
        let aligned_count = arrays.coverage_reduced.iter().filter(|&&x| x > 0).count();
        let aligned_fraction = if arrays.coverage_reduced.is_empty() {
            0.0
        } else {
            (aligned_count as f64 / arrays.coverage_reduced.len() as f64) * 100.0
        };

        if aligned_fraction >= pt_config.min_aligned_fraction as f64 {
            // Calculate global metrics for statistical filtering
            let total_coverage_reduced: u64 = arrays.coverage_reduced.iter().sum();
            let total_primary_reads: u64 = arrays.primary_reads.iter().sum();

            // Compute clipped_ratio = (total_primary_reads - total_coverage_reduced) / total_coverage_reduced
            // This represents the fraction of reads that are soft-clipped at any position on average.
            // For long reads, coverage_reduced can contain 2× the reads because each
            // read is split in two (start half and end half). Normalize by dividing by 2.
            let clipped_ratio = if total_coverage_reduced > 0 {
                let effective_coverage_reduced = if seq_type.is_long() {
                    total_coverage_reduced as f64 / 2.0
                } else {
                    total_coverage_reduced as f64
                };
                ((total_primary_reads as f64 - effective_coverage_reduced) / effective_coverage_reduced).max(0.0)
            } else {
                1.0
            };

            // Get clipping counts as u64 arrays
            let left_clip_counts: Vec<u64> = arrays.left_clipping_lengths.iter().map(|v| v.len() as u64).collect();
            let right_clip_counts: Vec<u64> = arrays.right_clipping_lengths.iter().map(|v| v.len() as u64).collect();

            // === STEP 3-5: Filter, statistical test, merge into areas, global filter ===
            // Returns ALL areas (both kept and discarded) with filtering metadata
            let (start_areas, end_areas) = filter_and_merge_to_areas_with_diagnostics(
                &arrays.reads_starts,
                &arrays.reads_ends,
                &left_clip_counts,
                &right_clip_counts,
                &arrays.coverage_reduced,
                &arrays.primary_reads,
                clipped_ratio,
                &pt_config,
                contig_length,
                config.circular,
                &dtr_regions,
            );

            // === STEP 6: Classify packaging based on peak-areas ===
            let (mechanism, left_termini, right_termini, duplication, repeat_length) = classify_packaging_areas(
                &start_areas,
                &end_areas,
                contig_length,
                config.circular,
                &dtr_regions,
            );

            // Only return PackagingData if there's a detected mechanism (not "No_packaging")
            // All termini (both kept and discarded) are included with filtering metadata
            if mechanism != "No_packaging" {
                Some(PackagingData {
                    contig_name: contig_name.to_string(),
                    mechanism,
                    left_termini,
                    right_termini,
                    duplication,
                    repeat_length,
                })
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    // Compute completeness statistics when mapping_metrics is enabled
    let completeness_result = if flags.mapping_metrics {
        // Use actual primary read count (not sum of coverage array which gives coverage × positions)
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
            arrays.circularising_reads_count,
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
/// Returns (features, presences, packaging, completeness, sample_name, seq_type, total_reads, mapped_reads)
pub fn process_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    repeats: &[RepeatsData],
) -> Result<(Vec<FeaturePoint>, Vec<PresenceData>, Vec<PackagingData>, Vec<CompletenessData>, String, SequencingType, u64, u64)> {
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

    // Get total read count from BAM index (fast, no iteration)
    let total_reads = get_total_read_count(bam_path)?;

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
            let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming(&mut bam, &ref_name, ref_length, seq_type, flags, config.circular, config.min_coverage) {
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
            let (packaging_info, completeness_info) = add_features_from_arrays(&mut arrays, &ref_name, ref_length, config, seq_type, flags, repeats, &mut features);

            // Calculate mean and median coverage depth
            let coverage_mean = arrays.coverage_mean() as f32;
            let coverage_median = arrays.coverage_median() as f32;

            // Calculate coverage variation using Fano factor style normalization
            // Formula: 1000000 × (1/(n-1)) × Σ(cov(i+1) - cov(i))² / mean_cov
            // Divides by mean (not mean²) to remove spurious correlation with coverage depth
            // Multiplied by 1000000 for 6 decimal digits precision when stored as integer
            let coverage_variation = if arrays.primary_reads.len() > 1 && coverage_mean > 0.0 {
                let n = arrays.primary_reads.len();
                let mean_cov = coverage_mean as f64;
                let sum_squared_diff: f64 = arrays.primary_reads
                    .windows(2)
                    .map(|w| {
                        let diff = w[1] as f64 - w[0] as f64;
                        diff * diff
                    })
                    .sum::<f64>() / (n - 1) as f64;
                let root_diff = sum_squared_diff.sqrt() / mean_cov;
                (root_diff * 1000000.0) as f32
            } else {
                0.0
            };

            // Coverage SD: Coefficient of Variation (CV) = std_dev / mean
            // Normalized to remove correlation with coverage depth
            // Scaled by 1,000,000 for integer storage (same as coverage_variation)
            let coverage_sd = if arrays.primary_reads.len() > 0 && coverage_mean > 0.0 {
                let mean_cov = coverage_mean as f64;
                let n = arrays.primary_reads.len() as f64;
                let variance: f64 = arrays.primary_reads
                    .iter()
                    .map(|&x| {
                        let diff = x as f64 - mean_cov;
                        diff * diff
                    })
                    .sum::<f64>() / n;
                let cv = variance.sqrt() / mean_cov;  // Coefficient of variation
                (cv * 1_000_000.0) as f32
            } else {
                0.0
            };

            let presence = PresenceData {
                contig_name: ref_name.clone(),
                coverage_pct: coverage_pct as f32,
                coverage_variation,
                coverage_sd,
                coverage_mean,
                coverage_median,
            };

            Some((features, presence, packaging_info, completeness_info, primary_count))
        })
        .collect();

    // Merge results from all contigs
    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();
    let mut all_packaging = Vec::new();
    let mut all_completeness = Vec::new();
    let mut mapped_reads: u64 = 0;

    for (features, presence, packaging, completeness, primary_count) in results {
        all_features.extend(features);
        all_presences.push(presence);
        mapped_reads += primary_count;
        if let Some(pkg) = packaging {
            all_packaging.push(pkg);
        }
        if let Some(comp) = completeness {
            all_completeness.push(comp);
        }
    }

    Ok((all_features, all_presences, all_packaging, all_completeness, sample_name, seq_type, total_reads, mapped_reads))
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
            sequence: None, // BAM headers don't contain sequence data
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
    config: &ProcessConfig,
    _create_indexes: bool, // Ignored - DuckDB uses zone maps instead of indexes
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
        parse_annotations(genbank_path)?
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

    // Parse and write repeats from autoblast file (if provided)
    let repeats = if !autoblast_file.as_os_str().is_empty() && autoblast_file.exists() {
        eprintln!("\n### Parsing autoblast results...");
        let reps = crate::db::parse_autoblast_file(autoblast_file)?;
        if !reps.is_empty() {
            db_writer.write_repeats(&reps)?;
        } else {
            eprintln!("No repeats found in autoblast file");
        }
        reps
    } else {
        Vec::new()
    };

    // Compute and write GC content and GC skew from sequence data (if available)
    // Uses configurable window sizes (default 500bp for GC content, 1000bp for GC skew)
    // and configurable RLE compression (default 0.1%)
    let gc_data: Vec<GCContentData> = contigs
        .iter()
        .filter_map(|contig| {
            contig.sequence.as_ref().map(|seq| {
                let (runs, stats) = compute_gc_content(seq, config.gc_params.gc_content_window_size, config.contig_variation_percentage);
                let (skew_runs, skew_stats) = compute_gc_skew(seq, config.gc_params.gc_skew_window_size, config.contig_variation_percentage);
                GCContentData {
                    contig_name: contig.name.clone(),
                    runs,
                    skew_runs,
                    stats,
                    skew_stats,
                }
            })
        })
        .collect();

    if !gc_data.is_empty() {
        eprintln!("\n### Computing GC content for {} contigs...", gc_data.len());
        db_writer.write_gc_content(&gc_data)?;
        db_writer.write_gc_skew(&gc_data)?;
        db_writer.update_contig_gc_stats(&gc_data)?;
    }

    // If no BAM files provided, skip sample processing (genbank-only mode)
    if bam_files.is_empty() {
        eprintln!("\n### No BAM files provided - skipping sample processing");
        eprintln!("Database populated with {} contigs and {} annotations", contigs.len(), annotations.len());
        
        return Ok(ProcessResult {
            samples_processed: 0,
            samples_failed: 0,
            total_time_secs: 0.0,
            processing_time_secs: 0.0,
            writing_time_secs: 0.0,
        });
    }

    eprintln!("\n### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!("Modules: {}\n", modules.join(", "));

    // Pre-validation: check all samples for warnings before processing starts
    validate_all_samples(bam_files, &contigs, modules, config.circular);

    let result = process_samples_parallel(&bam_files, &contigs, modules, config, db_writer, &repeats)?;

    print_summary(&result, output_db);

    Ok(result)
}

/// Holds processed sample data ready for database writing.
struct SampleResult {
    sample_name: String,
    sequencing_type: SequencingType,
    total_reads: u64,
    mapped_reads: u64,
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
    repeats: &[RepeatsData],
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let is_tty = atty::is(Stream::Stderr);
    let start_time = std::time::Instant::now();

    // Smart parallelization based on dataset characteristics:
    // - Many contigs (metagenomic ≥500): limit concurrency to prevent OOM
    // - Few contigs (phage/small genomes <500): full parallelism
    //
    // Sequential mode only for single-threaded operation
    if config.threads == 1 {
        return process_samples_sequential(bam_files, contigs, modules, config, db_writer, is_tty, repeats);
    }

    // Adaptive channel size based on memory pressure from contig count
    // - Metagenomic (≥500 contigs): channel_size=1 → process one sample at a time (prevents OOM)
    // - Normal (<500 contigs): channel_size=threads → full parallelism
    let channel_size = if contigs.len() >= 500 {
        1  // Large samples → process one at a time to prevent OOM
    } else {
        config.threads  // Normal samples → full parallelism
    };

    // Sample-parallel mode: producer-consumer with bounded channel
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

    // Bounded channel limits memory based on dataset characteristics
    let (tx, rx) = sync_channel::<SampleResult>(channel_size);

    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);

    // Clone values needed by writer thread
    let write_pb_clone = write_pb.clone();
    let is_tty_writer = is_tty;
    let total_writer = total;
    let circular = config.circular;

    // Spawn dedicated writer thread
    let writer_handle = thread::spawn(move || -> Result<(usize, std::time::Duration)> {
        let write_start = std::time::Instant::now();
        let mut written_count = 0usize;

        // Receive and write samples as they arrive
        for result in rx {
            written_count += 1;

            // Insert sample
            if let Err(e) = db_writer.insert_sample(&result.sample_name, result.sequencing_type.as_str(), result.total_reads, result.mapped_reads) {
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
                circular,
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
        write_pb_clone.finish();

        Ok((written_count, write_start.elapsed()))
    });

    // Process samples in parallel, sending to channel immediately
    // send() blocks if channel is full (backpressure)
    bam_files.par_iter().for_each(|bam_path| {
        let sample_start = std::time::Instant::now();

        match process_sample(bam_path, contigs, modules, config, repeats) {
            Ok((features, presences, packaging, completeness, sample_name, sequencing_type, total_reads, mapped_reads)) => {
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
                    sequencing_type,
                    total_reads,
                    mapped_reads,
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

    process_pb.finish();
    let processing_time = start_time.elapsed();

    // Wait for writer thread to finish
    let (written_count, writing_time) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    // Clear MultiProgress to avoid duplicate bar display
    mp.clear().ok();

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
    repeats: &[RepeatsData],
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

        match process_sample(bam_path, contigs, modules, config, repeats) {
            Ok((features, presences, packaging, completeness, sample_name, sequencing_type, total_reads, mapped_reads)) => {
                let process_elapsed = sample_start.elapsed();
                processing_time_total += process_elapsed;

                let write_start = std::time::Instant::now();

                // Insert sample
                if let Err(e) = db_writer.insert_sample(&sample_name, sequencing_type.as_str(), total_reads, mapped_reads) {
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
                    config.circular,
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
