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

use crate::bam_reader::process_contig_streaming;
use crate::compress::{
    compress_signal_with_reference, Run,
};
use crate::db::{
    create_metadata_db, create_temp_sample_db, finalize_db, merge_temp_db_into_main,
    write_features_to_temp_db, write_presences_to_temp_db,
};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::genbank::parse_genbank;
use crate::types::{
    get_plot_type, ContigInfo, FeaturePoint, PlotType, PresenceData, SequencingType
};

/// Merge consecutive runs with identical values (0% tolerance RLE).
/// 
/// Applied to features that should be constant along entire reads (deletions, mate flags).
/// Only merges runs that are BOTH adjacent (end+1 == start) AND have the same value.
#[inline]
fn merge_identical_runs(runs: Vec<Run>) -> Vec<Run> {
    if runs.is_empty() {
        return runs;
    }
    
    let mut merged = Vec::new();
    let mut current = runs[0].clone();
    
    for run in runs.into_iter().skip(1) {
        // Only merge if runs are adjacent AND have the same value
        if current.end_pos + 1 == run.start_pos && (run.value - current.value).abs() < f32::EPSILON {
            // Adjacent runs with same value: extend the current run
            current.end_pos = run.end_pos;
        } else {
            // Either not adjacent or different value: push current and start new
            merged.push(current);
            current = run;
        }
    }
    merged.push(current);
    
    merged
}

/// Peak for phage termini detection.
#[derive(Clone, Debug)]
struct Peak {
    position: i32,
    value: f64,
}

/// Calculate circular distance between two positions.
fn circular_distance(pos1: i32, pos2: i32, genome_length: usize) -> i32 {
    let gl = genome_length as i32;
    let dist_forward = if pos2 >= pos1 {
        pos2 - pos1
    } else {
        (gl - pos1) + pos2
    };
    let dist_backward = if pos1 >= pos2 {
        pos1 - pos2
    } else {
        (gl - pos2) + pos1
    };
    dist_forward.min(dist_backward)
}

/// Filters out peaks that have clipping events within max_distance.
fn has_nearby_clip(
    peak_pos: i32,
    clip_positions: &[u32],
    min_distance: i32,
    genome_length: usize,
    circular: bool,
) -> bool {
    for &clip in clip_positions {
        let clip_pos = clip as i32;

        let dist = if circular {
            circular_distance(peak_pos, clip_pos, genome_length)
        } else {
            (peak_pos - clip_pos).abs()
        };

        if dist <= min_distance {
            return true;
        }
    }
    false
}

fn filter_peaks_by_spc_and_clippings(
    reads_starts: &[Run],
    reads_ends: &[Run],
    min_spc: f64,
    min_distance: i32,
    genome_length: usize,
    circular: bool,
    left_clip_pos: &Option<Vec<u32>>,
    right_clip_pos: &Option<Vec<u32>>,
) -> (Vec<Run>, Vec<Run>) {
    // flatten clip positions once (cheaper & simpler)
    let left_clips: Vec<u32> = left_clip_pos.iter().flatten().copied().collect();
    let right_clips: Vec<u32> = right_clip_pos.iter().flatten().copied().collect();

    // filter reads_starts using right clipping
    let filtered_starts: Vec<Run> = reads_starts
        .iter()
        .filter(|run| {
            run.value as f64 >= min_spc &&
            !has_nearby_clip(
                run.start_pos,
                &right_clips,
                min_distance,
                genome_length,
                circular,
            )
        })
        .cloned()
        .collect();

    // filter reads_ends using left clipping
    let filtered_ends: Vec<Run> = reads_ends
        .iter()
        .filter(|run| {
            run.value as f64 >= min_spc &&
            !has_nearby_clip(
                run.start_pos,
                &left_clips,
                min_distance,
                genome_length,
                circular,
            )
        })
        .cloned()
        .collect();

    (filtered_starts, filtered_ends)
}

/// Merge nearby peaks within max_distance, keeping the highest value peak.
fn merge_peaks(
    runs: &[Run],
    max_distance: i32,
    genome_length: usize,
    circular: bool,
) -> Vec<Peak> {
    if runs.is_empty() {
        return Vec::new();
    }

    // Convert runs to peaks
    let mut peaks: Vec<Peak> = runs
        .iter()
        .map(|r| Peak {
            position: r.start_pos,
            value: r.value as f64,
        })
        .collect();

    // Sort by genomic position
    peaks.sort_by_key(|p| p.position);

    let mut merged: Vec<Peak> = Vec::new();
    let mut current = peaks[0].clone();

    for peak in peaks.into_iter().skip(1) {
        let dist = if circular {
            circular_distance(current.position, peak.position, genome_length)
        } else {
            peak.position - current.position
        };

        if dist <= max_distance {
            // keep strongest
            if peak.value > current.value {
                current = peak;
            }
        } else {
            merged.push(current);
            current = peak;
        }
    }
    merged.push(current);

    // Circular wrap-around merge (first & last)
    if circular && merged.len() > 1 {
        let first = &merged[0];
        let last = &merged[merged.len() - 1];

        let wrap_dist = circular_distance(last.position, first.position, genome_length);

        if wrap_dist <= max_distance {
            let keeper = if first.value >= last.value {
                first.clone()
            } else {
                last.clone()
            };

            merged.remove(merged.len() - 1);
            merged[0] = keeper;
        }
    }

    merged
}

/// Classify phage packaging mechanism based on peak configuration.
fn classify_packaging(
    start_peaks: &[Peak],
    end_peaks: &[Peak],
    genome_length: usize,
    circular: bool,
) -> (String, Option<i32>, Option<i32>) {
    match (start_peaks.len(), end_peaks.len()) {
        (0, 0) => ("No_packaging".to_string(), None, None),
        
        (1, 0) => {
            // Only start peak - PAC
            ("PAC".to_string(), Some(start_peaks[0].position), None)
        }
        
        (0, 1) => {
            // Only end peak - PAC
            ("PAC".to_string(), None, Some(end_peaks[0].position))
        }
        
        (1, 1) => {
            // Two peaks - classify based on distance and order
            let start_pos = start_peaks[0].position;
            let end_pos = end_peaks[0].position;
            
            // Calculate actual genomic distance (shortest path in circular)
            let distance = if circular {
                circular_distance(start_pos, end_pos, genome_length)
            } else {
                (start_pos - end_pos).abs()
            };
            
            let genome_10pct = (genome_length as f64 * 0.1) as i32;
            
            // Determine order: which comes first going forward (clockwise in circular)
            // For linear: simply check positions
            // For circular: check which direction is shorter to go from start to end
            let end_before_start = if circular {
                // Calculate forward distance from start to end
                let gl = genome_length as i32;
                let forward_dist = if end_pos >= start_pos {
                    end_pos - start_pos
                } else {
                    (gl - start_pos) + end_pos
                };
                // If forward distance equals circular_distance, end is ahead
                // If they're different, end is behind (going backward is shorter)
                forward_dist != distance
            } else {
                end_pos < start_pos
            };
            
            let (mechanism, left, right) = if distance < 2 {
                // Very close - cohesive ends
                ("COS".to_string(), Some(start_pos.min(end_pos)), Some(start_pos.max(end_pos)))
            } else if distance <= 20 {
                // Close cohesive ends with directionality
                if end_before_start {
                    ("COS_3'".to_string(), Some(end_pos), Some(start_pos))
                } else {
                    ("COS_5'".to_string(), Some(start_pos), Some(end_pos))
                }
            } else if distance <= 1000 {
                // Short direct terminal repeats
                if end_before_start {
                    ("DTR_short_3'".to_string(), Some(end_pos), Some(start_pos))
                } else {
                    ("DTR_short_5'".to_string(), Some(start_pos), Some(end_pos))
                }
            } else if distance <= genome_10pct {
                // Long direct terminal repeats
                if end_before_start {
                    ("DTR_long_3'".to_string(), Some(end_pos), Some(start_pos))
                } else {
                    ("DTR_long_5'".to_string(), Some(start_pos), Some(end_pos))
                }
            } else {
                // Very distant - outlier
                if end_before_start {
                    ("DTR_outlier_3'".to_string(), Some(end_pos), Some(start_pos))
                } else {
                    ("DTR_outlier_5'".to_string(), Some(start_pos), Some(end_pos))
                }
            };
            
            (mechanism, left, right)
        }
        
        _ => {
            // Multiple peaks or other configurations
            ("Unknown_packaging".to_string(), None, None)
        }
    }
}

/// Configuration for phage termini detection.
#[derive(Clone, Copy)]
pub struct PhageTerminiConfig {
    /// Minimum SPC (Significant Peak Count) to consider a peak significant
    pub min_spc: f64,
    /// Maximum distance (bp) between peaks to merge them
    pub max_distance_peaks: i32,
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_spc: 10.0,
            max_distance_peaks: 20,
        }
    }
}

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
    /// Sequencing type: determines which features to calculate
    pub sequencing_type: SequencingType,
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
) -> Vec<Run> {
    let plot_type = get_plot_type(feature);
    let runs = compress_signal_with_reference(values, reference, plot_type, config.curve_ratio, config.bar_ratio);

    output.extend(runs.iter().map(|run| FeaturePoint {
        contig_name: contig_name.to_string(),
        feature: feature.to_string(),
        start_pos: run.start_pos,
        end_pos: run.end_pos,
        value: run.value,
        mean: None,
        median: None,
        std: None,
    }));

    runs
}

/// Compress and add feature points with statistics (for clippings/insertions).
///
/// Includes mean, median, and standard deviation from the length vectors.
#[inline]
fn add_compressed_feature_with_stats(
    counts: &[f64],
    means: &[f64],
    medians: &[f64],
    stds: &[f64],
    reference: Option<&[f64]>,
    feature: &str,
    contig_name: &str,
    config: &ProcessConfig,
    output: &mut Vec<FeaturePoint>,
) -> Vec<Run> {
    let plot_type = get_plot_type(feature);
    let runs = compress_signal_with_reference(counts, reference, plot_type, config.curve_ratio, config.bar_ratio);

    output.extend(runs.iter().map(|run| {
        // For the run's position range, compute average statistics
        let start_idx = (run.start_pos - 1) as usize;
        let end_idx = run.end_pos as usize;
        
        let (mean_val, median_val, std_val) = if start_idx < means.len() {
            let range_mean: f64 = means[start_idx..end_idx.min(means.len())].iter().sum::<f64>() 
                / (end_idx - start_idx) as f64;
            let range_median: f64 = medians[start_idx..end_idx.min(medians.len())].iter().sum::<f64>() 
                / (end_idx - start_idx) as f64;
            let range_std: f64 = stds[start_idx..end_idx.min(stds.len())].iter().sum::<f64>() 
                / (end_idx - start_idx) as f64;
            (Some(range_mean as f32), Some(range_median as f32), Some(range_std as f32))
        } else {
            (None, None, None)
        };

        FeaturePoint {
            contig_name: contig_name.to_string(),
            feature: feature.to_string(),
            start_pos: run.start_pos,
            end_pos: run.end_pos,
            value: run.value,
            mean: mean_val,
            median: median_val,
            std: std_val,
        }
    }));

    runs
}

/// Add features from FeatureArrays to output (optimized path).
/// Returns optional packaging classification for phagetermini module.
fn add_features_from_arrays(
    arrays: &FeatureArrays,
    contig_name: &str,
    contig_length: usize,
    config: &ProcessConfig,
    flags: ModuleFlags,
    output: &mut Vec<FeaturePoint>,
) -> Option<(String, Option<i32>, Option<i32>)> {
    let seq_type = config.sequencing_type;
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
    }

    // Assemblycheck features
    let mut left_clip_pos: Option<Vec<u32>> = None;
    let mut right_clip_pos: Option<Vec<u32>> = None;
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
        let left_clip_run = add_compressed_feature_with_stats(&left_clip_counts, &left_clip_means, &left_clip_medians, &left_clip_stds,
            Some(&primary_reads_f64), "left_clippings", contig_name, config, output);
        left_clip_pos = Some(left_clip_run.iter().map(|r| r.start_pos as u32).collect());

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
        let right_clip_run = add_compressed_feature_with_stats(&right_clip_counts, &right_clip_means, &right_clip_medians, &right_clip_stds,
            Some(&primary_reads_f64), "right_clippings", contig_name, config, output);
        right_clip_pos = Some(right_clip_run.iter().map(|r| r.start_pos as u32).collect());
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
        add_compressed_feature_with_stats(&insertion_counts, &insertion_means, &insertion_medians, &insertion_stds,
            Some(&primary_reads_f64), "insertions", contig_name, config, output);

        // Other assembly check features (no statistics)
        // Apply two-stage compression: first with coverage reference, then merge identical runs
        let deletions_f64: Vec<f64> = arrays.deletions.iter().map(|&x| x as f64).collect();
        let deletions_runs = compress_signal_with_reference(&deletions_f64, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        let deletions_runs_merged = merge_identical_runs(deletions_runs);
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
        add_compressed_feature_with_reference(&mismatches_f64, Some(&primary_reads_f64), "mismatches", contig_name, config, output);
        
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
    // Coverage_reduced for phagetermini (use as reference for phagetermini bar features)
    let coverage_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
    let packaging_result = if flags.phagetermini {
        // coverage_reduced (self-referential curve)
        add_compressed_feature(&coverage_reduced_f64, "coverage_reduced", contig_name, config, output);
        
        // reads_starts and reads_ends (bars)
        // use primary_reads as reference instead of coverage_reduced to give more weight to mapping errors revealing assembly errors that can give wrong termini
        let reads_starts: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_starts_runs = compress_signal_with_reference(&reads_starts, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);
        
        let reads_ends: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let reads_ends_runs = compress_signal_with_reference(&reads_ends, Some(&primary_reads_f64), PlotType::Bars, config.curve_ratio, config.bar_ratio);

        // Add reads_starts and reads_ends to output
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

        // Tau: calculate tau values directly for positions where reads_starts or reads_ends were saved
        // Merge the runs from both features and calculate tau for those positions
        for run in reads_starts_runs.iter().chain(reads_ends_runs.iter()) {
            // Calculate average tau value for this run
            let mut tau_sum = 0.0;
            let mut count = 0;
            
            for pos in run.start_pos..=run.end_pos {
                let idx = (pos - 1) as usize; // Convert from 1-indexed to 0-indexed
                if idx < primary_reads_f64.len() {
                    let c = primary_reads_f64[idx];
                    let s = arrays.reads_starts[idx];
                    let e = arrays.reads_ends[idx];
                    if c > 0.0 {
                        tau_sum += (s + e) as f64 / c;
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

        // Classify packaging mechanism
        let pt_config = config.phagetermini_config;

        // Filter runs by minimum SPC threshold AND no clipping events within max_distance
        let (reads_starts_filtered, reads_ends_filtered) = filter_peaks_by_spc_and_clippings(
            &reads_starts_runs,
            &reads_ends_runs,
            pt_config.min_spc,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular,
            &left_clip_pos,
            &right_clip_pos,
        );

        let start_peaks = merge_peaks(
            &reads_starts_filtered,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular
        );
        let end_peaks = merge_peaks(
            &reads_ends_filtered,
            pt_config.max_distance_peaks,
            contig_length,
            config.circular
        );

        let (mechanism, left_terminus, right_terminus) = classify_packaging(
            &start_peaks,
            &end_peaks,
            contig_length,
            config.circular,
        );

        Some((mechanism, left_terminus, right_terminus))
    } else {
        None
    };

    packaging_result
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
            let (arrays, coverage_pct) = match process_contig_streaming(&mut bam, &ref_name, ref_length, config.sequencing_type, flags, config.circular, config.min_coverage) {
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
            let packaging_info = add_features_from_arrays(&arrays, &ref_name, ref_length, config, flags, &mut features);

            let presence = PresenceData {
                contig_name: ref_name.clone(),
                coverage_pct: coverage_pct as f32,
                phage_packaging_mechanism: packaging_info.as_ref().map(|(m, _, _)| m.clone()),
                phage_left_terminus: packaging_info.as_ref().and_then(|(_, l, _)| *l),
                phage_right_terminus: packaging_info.as_ref().and_then(|(_, _, r)| *r),
            };

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
    let processing_time = start_time.elapsed();

    // Merge all temp DBs sequentially (fast operation, not a bottleneck)
    let merge_start = std::time::Instant::now();
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

    let writing_time = merge_start.elapsed();
    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);

    Ok(ProcessResult {
        samples_processed: total - failed,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
        processing_time_secs: processing_time.as_secs_f64(),
        writing_time_secs: writing_time.as_secs_f64(),
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
