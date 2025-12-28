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

use crate::bam_reader::process_contig_streaming;
use crate::compress::{
    compress_signal_with_reference, Run,
};
use crate::db::{DbWriter, CompletenessData, DuplicationData};
use crate::features::{FeatureArrays, ModuleFlags};
use crate::genbank::parse_genbank;
use crate::types::{
    get_plot_type, ContigInfo, FeaturePoint, PackagingData, PlotType, PresenceData, SequencingType
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

/// DTR region info for distance calculations.
/// first_start/first_end: positions of the first (kept) duplicated region
/// second_start/second_end: positions of the second (zeroed) duplicated region
#[derive(Clone, Debug)]
struct DtrRegion {
    first_start: i32,
    first_end: i32,
    second_start: i32,
    second_end: i32,
    /// true = DTR (direct), false = ITR (inverted)
    is_direct: bool,
}

/// Calculate minimum distance between two positions, considering DTR equivalence.
///
/// If a position is near the second DTR region, it's treated as being near
/// the equivalent position in the first DTR region.
fn dtr_aware_distance(
    pos1: i32,
    pos2: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
    max_distance: i32,
) -> i32 {
    // Start with regular distance
    let regular_dist = if circular {
        circular_distance(pos1, pos2, genome_length)
    } else {
        (pos1 - pos2).abs()
    };

    // If no DTR regions or already close enough, return regular distance
    if dtr_regions.is_empty() || regular_dist <= max_distance {
        return regular_dist;
    }

    // Check if either position is near a second DTR region
    // If so, calculate distance using the equivalent first region position
    let mut min_dist = regular_dist;

    for dtr in dtr_regions {
        // Check if pos1 is near second region - calculate virtual position
        let pos1_virtual = get_virtual_position(pos1, dtr, max_distance);
        // Check if pos2 is near second region - calculate virtual position
        let pos2_virtual = get_virtual_position(pos2, dtr, max_distance);

        // Calculate distance using virtual positions
        let virtual_dist = if circular {
            circular_distance(pos1_virtual, pos2_virtual, genome_length)
        } else {
            (pos1_virtual - pos2_virtual).abs()
        };

        min_dist = min_dist.min(virtual_dist);
    }

    min_dist
}

/// Get the virtual position for a point near a DTR region.
/// If the position is near the second DTR region, return the equivalent
/// position near the first DTR region. Otherwise, return the original position.
fn get_virtual_position(pos: i32, dtr: &DtrRegion, max_distance: i32) -> i32 {
    // Distance from pos to the start of second region
    let dist_to_second_start = (pos - dtr.second_start).abs();
    // Distance from pos to the end of second region
    let dist_to_second_end = (pos - dtr.second_end).abs();

    if dist_to_second_start <= max_distance {
        // Position is near second region's start
        // Map to equivalent position near first region's start
        let offset = pos - dtr.second_start;
        dtr.first_start + offset
    } else if dist_to_second_end <= max_distance {
        // Position is near second region's end
        // Map to equivalent position near first region's end
        let offset = pos - dtr.second_end;
        dtr.first_end + offset
    } else if pos >= dtr.second_start && pos <= dtr.second_end {
        // Position is inside second region
        let offset = pos - dtr.second_start;
        dtr.first_start + offset
    } else {
        // Position is not near second region, keep as is
        pos
    }
}

/// Filters out peaks that have clipping events within max_distance.
fn has_nearby_clip(
    peak_pos: i32,
    clip_positions: &[u32],
    min_distance: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> bool {
    for &clip in clip_positions {
        let clip_pos = clip as i32;

        let dist = dtr_aware_distance(
            peak_pos,
            clip_pos,
            genome_length,
            circular,
            dtr_regions,
            min_distance,
        );

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
    dtr_regions: &[DtrRegion],
) -> (Vec<Run>, Vec<Run>) {
    // flatten clip positions once (cheaper & simpler)
    let left_clips: Vec<u32> = left_clip_pos.iter().flatten().copied().collect();
    let right_clips: Vec<u32> = right_clip_pos.iter().flatten().copied().collect();

    // filter reads_starts using left clipping only
    let filtered_starts: Vec<Run> = reads_starts
        .iter()
        .filter(|run| {
            run.value as f64 >= min_spc &&
            !has_nearby_clip(
                run.start_pos,
                &left_clips,
                min_distance,
                genome_length,
                circular,
                dtr_regions,
            )
        })
        .cloned()
        .collect();

    // filter reads_ends using right clipping only
    let filtered_ends: Vec<Run> = reads_ends
        .iter()
        .filter(|run| {
            run.value as f64 >= min_spc &&
            !has_nearby_clip(
                run.start_pos,
                &right_clips,
                min_distance,
                genome_length,
                circular,
                dtr_regions,
            )
        })
        .cloned()
        .collect();

    (filtered_starts, filtered_ends)
}

/// Merge nearby peaks within max_distance, keeping the highest value peak.
/// Considers DTR-equivalent positions when calculating distances.
fn merge_peaks(
    runs: &[Run],
    max_distance: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
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
        let dist = dtr_aware_distance(
            current.position,
            peak.position,
            genome_length,
            circular,
            dtr_regions,
            max_distance,
        );

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

    // Circular wrap-around merge (first & last) - also DTR-aware
    if circular && merged.len() > 1 {
        let first = &merged[0];
        let last = &merged[merged.len() - 1];

        let wrap_dist = dtr_aware_distance(
            last.position,
            first.position,
            genome_length,
            circular,
            dtr_regions,
            max_distance,
        );

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

/// Check if a duplication is valid for phage termini analysis.
/// A valid duplication must have:
/// - One region starting within max_distance from position 1 (beginning)
/// - The other region ending within max_distance from contig_length (end)
fn is_valid_terminal_duplication(dup: &DuplicationData, contig_length: usize, max_distance: i32) -> bool {
    let first_start = dup.position1.min(dup.position2);
    let first_end = dup.position1.max(dup.position2);
    let second_start = dup.position1prime.min(dup.position2prime);
    let second_end = dup.position1prime.max(dup.position2prime);

    let contig_end = contig_length as i32;

    // Check if first region is at beginning and second region is at end
    let first_at_start = first_start <= max_distance;
    let second_at_end = second_end >= (contig_end - max_distance);

    // Check the reverse: first at end, second at beginning
    let first_at_end = first_end >= (contig_end - max_distance);
    let second_at_start = second_start <= max_distance;

    (first_at_start && second_at_end) || (first_at_end && second_at_start)
}

/// Translate a terminus position from second region to first region if applicable.
/// This ensures termini are always reported in the first (canonical) region.
/// Returns the translated position (or original if not in second region).
///
/// Mapping (regions are same length for true duplications):
/// - DTR: second_start → first_start, second_end → first_end
/// - ITR: second_start → first_end, second_end → first_start
fn translate_to_first_region(pos: i32, dtr_regions: &[DtrRegion]) -> i32 {
    for dtr in dtr_regions {
        // Check if pos is in second region
        if pos >= dtr.second_start && pos <= dtr.second_end {
            return dtr.first_start + pos - dtr.second_start;
        }
    }
    // Not in any second region, return as-is
    pos
}

/// Expand a terminus position to include equivalent DTR positions.
/// First translates to first region if in second region, then adds the equivalent position.
fn expand_terminus_with_dtr(pos: i32, dtr_regions: &[DtrRegion]) -> Vec<i32> {
    // First, translate to first region if in second region
    let primary_pos = translate_to_first_region(pos, dtr_regions);
    let mut positions = vec![primary_pos];

    for dtr in dtr_regions {
        // Check if primary_pos is in first region - add equivalent in second region
        if primary_pos >= dtr.first_start && primary_pos <= dtr.first_end {
            let offset = primary_pos - dtr.first_start;
            let equivalent = dtr.second_start + offset;
            if !positions.contains(&equivalent) {
                positions.push(equivalent);
            }
        }
    }

    positions.sort();
    positions
}

/// Check if a position is within any duplication region.
/// Returns Some(is_direct) if in a region, None otherwise.
fn position_in_duplication(pos: i32, dtr_regions: &[DtrRegion]) -> Option<bool> {
    for dtr in dtr_regions {
        if (pos >= dtr.first_start && pos <= dtr.first_end) || (pos >= dtr.second_start && pos <= dtr.second_end)
        {
            return Some(dtr.is_direct);
        }
    }
    None
}

/// Determine duplication status for all termini.
/// Returns Some(true) if ALL termini are in DTR regions,
/// Some(false) if ALL termini are in ITR regions,
/// None if no termini, mixed, or any terminus outside duplications.
fn determine_duplication_status(
    left_termini: &[i32],
    right_termini: &[i32],
    dtr_regions: &[DtrRegion],
) -> Option<bool> {
    // Collect all terminus positions (use just the first/original position, not expanded)
    let all_termini: Vec<i32> = left_termini.iter().chain(right_termini.iter()).copied().collect();

    if all_termini.is_empty() || dtr_regions.is_empty() {
        return None;
    }

    let mut all_direct = true;
    let mut all_inverted = true;

    for &pos in &all_termini {
        // Translate to first region before checking (in case position is in second region)
        let translated_pos = translate_to_first_region(pos, dtr_regions);
        match position_in_duplication(translated_pos, dtr_regions) {
            Some(true) => all_inverted = false,  // Found DTR
            Some(false) => all_direct = false,   // Found ITR
            None => return None,                  // Position not in any duplication
        }
    }

    if all_direct && !all_inverted {
        Some(true)  // All DTR
    } else if all_inverted && !all_direct {
        Some(false) // All ITR
    } else {
        None // Mixed or none
    }
}

/// Classify phage packaging mechanism based on peak configuration.
/// Returns (mechanism, left_termini, right_termini, duplication_status).
fn classify_packaging(
    start_peaks: &[Peak],
    end_peaks: &[Peak],
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (String, Vec<i32>, Vec<i32>, Option<bool>) {
    match (start_peaks.len(), end_peaks.len()) {
        (0, 0) => ("No_packaging".to_string(), vec![], vec![], None),

        (1, 0) => {
            // Only start peak - PAC
            let left = expand_terminus_with_dtr(start_peaks[0].position, dtr_regions);
            let dup_status = determine_duplication_status(&[start_peaks[0].position], &[], dtr_regions);
            ("PAC".to_string(), left, vec![], dup_status)
        }

        (0, 1) => {
            // Only end peak - PAC
            let right = expand_terminus_with_dtr(end_peaks[0].position, dtr_regions);
            let dup_status = determine_duplication_status(&[], &[end_peaks[0].position], dtr_regions);
            ("PAC".to_string(), vec![], right, dup_status)
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

            let (mechanism, left_pos, right_pos) = if distance < 2 {
                // Very close - cohesive ends
                ("COS".to_string(), start_pos.min(end_pos), start_pos.max(end_pos))
            } else if distance <= 20 {
                // Close cohesive ends with directionality
                if end_before_start {
                    ("COS_3'".to_string(), end_pos, start_pos)
                } else {
                    ("COS_5'".to_string(), start_pos, end_pos)
                }
            } else if distance <= 1000 {
                // Short direct terminal repeats
                if end_before_start {
                    ("DTR_short_3'".to_string(), end_pos, start_pos)
                } else {
                    ("DTR_short_5'".to_string(), start_pos, end_pos)
                }
            } else if distance <= genome_10pct {
                // Long direct terminal repeats
                if end_before_start {
                    ("DTR_long_3'".to_string(), end_pos, start_pos)
                } else {
                    ("DTR_long_5'".to_string(), start_pos, end_pos)
                }
            } else {
                // Very distant - outlier
                if end_before_start {
                    ("DTR_outlier_3'".to_string(), end_pos, start_pos)
                } else {
                    ("DTR_outlier_5'".to_string(), start_pos, end_pos)
                }
            };

            let left = expand_terminus_with_dtr(left_pos, dtr_regions);
            let right = expand_terminus_with_dtr(right_pos, dtr_regions);
            let dup_status = determine_duplication_status(&[left_pos], &[right_pos], dtr_regions);
            (mechanism, left, right, dup_status)
        }

        _ => {
            // Multiple peaks or other configurations - collect all original peak positions
            let left_original: Vec<i32> = start_peaks.iter().map(|p| p.position).collect();
            let right_original: Vec<i32> = end_peaks.iter().map(|p| p.position).collect();

            // Expand to include DTR equivalents
            let mut left_termini: Vec<i32> = start_peaks
                .iter()
                .flat_map(|p| expand_terminus_with_dtr(p.position, dtr_regions))
                .collect();
            let mut right_termini: Vec<i32> = end_peaks
                .iter()
                .flat_map(|p| expand_terminus_with_dtr(p.position, dtr_regions))
                .collect();

            // Remove duplicates and sort
            left_termini.sort();
            left_termini.dedup();
            right_termini.sort();
            right_termini.dedup();

            let dup_status = determine_duplication_status(&left_original, &right_original, dtr_regions);
            ("Unknown_packaging".to_string(), left_termini, right_termini, dup_status)
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
    /// Maximum distance (bp) from reference ends for duplication regions to be considered valid.
    /// Duplications must have one region starting within this distance from position 1
    /// and the other region ending within this distance from the contig end.
    pub max_distance_duplication: i32,
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_spc: 10.0,
            max_distance_peaks: 20,
            max_distance_duplication: 100,
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

/// Compute completeness statistics from clipping, insertion, mismatch, and deletion data.
///
/// For the left side: finds the first left-clipping event where clipped_length > distance_to_start.
///
/// For the right side: finds the first right-clipping event where clipped_length > distance_to_end.
///
/// score_completeness: sum of weighted missing basepairs (in reads but not reference):
/// - clippings: (count/coverage) * median_length
/// - insertions: (count/coverage) * median_length
/// - mismatches: (count/coverage) * 1
///
/// score_contamination: sum of weighted contamination basepairs (in reference but not reads):
/// - mismatches: (count/coverage) * 1
/// - deletions: (count/coverage) * median_length
/// - paired clippings (right followed by left): average_prevalence * distance
///
/// Returns None for a side if the first clipping event doesn't satisfy the condition.
fn compute_completeness(
    left_clipping_lengths: &[Vec<u32>],
    right_clipping_lengths: &[Vec<u32>],
    insertion_lengths: &[Vec<u32>],
    deletion_lengths: &[Vec<u32>],
    primary_reads: &[u64],
    left_clip_runs: &[Run],
    right_clip_runs: &[Run],
    insertion_runs: &[Run],
    deletion_runs: &[Run],
    mismatch_runs: &[Run],
    contig_name: &str,
    contig_length: usize,
) -> CompletenessData {
    // Helper function to compute median from a vector of lengths
    fn compute_median(lengths: &[u32]) -> i32 {
        if lengths.is_empty() {
            return 0;
        }
        let mut sorted = lengths.to_vec();
        sorted.sort_unstable();
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            ((sorted[mid - 1] + sorted[mid]) / 2) as i32
        } else {
            sorted[mid] as i32
        }
    }

    // Find the FIRST left-clipping event (leftmost position), then check if clipped_length > distance_to_start
    // Only valid if the first left-clipping satisfies the condition
    let left_result = left_clip_runs
        .iter()
        .min_by_key(|run| run.start_pos) // Find the first (leftmost) left-clipping
        .and_then(|run| {
            let pos = (run.start_pos - 1) as usize; // Convert to 0-indexed
            let distance_to_start = run.start_pos - 1; // Distance from position 1
            let min_missing = if pos < left_clipping_lengths.len() {
                compute_median(&left_clipping_lengths[pos])
            } else {
                0
            };
            // Only valid if clipped length > distance to start
            if min_missing > distance_to_start {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                Some((prevalence, distance_to_start, min_missing))
            } else {
                None
            }
        });

    // Find the FIRST right-clipping event (rightmost position), then check if clipped_length > distance_to_end
    // Only valid if the first right-clipping satisfies the condition
    let right_result = right_clip_runs
        .iter()
        .max_by_key(|run| run.start_pos) // Find the first (rightmost) right-clipping
        .and_then(|run| {
            let pos = (run.start_pos - 1) as usize; // Convert to 0-indexed
            let distance_to_end = contig_length as i32 - run.start_pos; // Distance to end
            let min_missing = if pos < right_clipping_lengths.len() {
                compute_median(&right_clipping_lengths[pos])
            } else {
                0
            };
            // Only valid if clipped length > distance to end
            if min_missing > distance_to_end {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                Some((prevalence, distance_to_end, min_missing))
            } else {
                None
            }
        });

    // Compute individual score components

    // Total mismatches: Σ(count/coverage) - each mismatch contributes 1 bp
    let mut total_mismatches = 0.0;
    for run in mismatch_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if coverage > 0.0 {
            total_mismatches += count / coverage;
        }
    }

    // Total deletions: Σ(count/coverage * median_length)
    let mut total_deletions = 0.0;
    for run in deletion_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if coverage > 0.0 && pos < deletion_lengths.len() {
            let length = compute_median(&deletion_lengths[pos]) as f64;
            total_deletions += (count / coverage) * length;
        }
    }

    // Total insertions: Σ(count/coverage * median_length)
    let mut total_insertions = 0.0;
    for run in insertion_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if coverage > 0.0 && pos < insertion_lengths.len() {
            let length = compute_median(&insertion_lengths[pos]) as f64;
            total_insertions += (count / coverage) * length;
        }
    }

    // Total reads clipped: Σ(count/coverage * median_length) for left + right clippings
    let mut total_reads_clipped = 0.0;
    for run in left_clip_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if coverage > 0.0 && pos < left_clipping_lengths.len() {
            let length = compute_median(&left_clipping_lengths[pos]) as f64;
            total_reads_clipped += (count / coverage) * length;
        }
    }
    for run in right_clip_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if coverage > 0.0 && pos < right_clipping_lengths.len() {
            let length = compute_median(&right_clipping_lengths[pos]) as f64;
            total_reads_clipped += (count / coverage) * length;
        }
    }

    // Total reference clipped: Σ(avg_prevalence * distance) for paired clips (right followed by left)
    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum ClipType { Right, Left }

    let mut all_clips: Vec<(i32, ClipType, f64)> = Vec::new(); // (position, type, prevalence)

    for run in right_clip_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        let prevalence = if coverage > 0.0 { count / coverage } else { 0.0 };
        all_clips.push((run.start_pos, ClipType::Right, prevalence));
    }

    for run in left_clip_runs {
        let pos = (run.start_pos - 1) as usize;
        let count = run.value as f64;
        let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        let prevalence = if coverage > 0.0 { count / coverage } else { 0.0 };
        all_clips.push((run.start_pos, ClipType::Left, prevalence));
    }

    // Sort by position
    all_clips.sort_by_key(|(pos, _, _)| *pos);

    // Find consecutive pairs: right-clipping followed immediately by left-clipping
    let mut total_reference_clipped = 0.0;
    for i in 0..all_clips.len().saturating_sub(1) {
        let (pos_right, type_right, prev_right) = all_clips[i];
        let (pos_left, type_left, prev_left) = all_clips[i + 1];

        if type_right == ClipType::Right && type_left == ClipType::Left {
            let distance = pos_left - pos_right;
            if distance > 0 {
                let avg_prevalence = (prev_right + prev_left) / 2.0;
                total_reference_clipped += avg_prevalence * distance as f64;
            }
        }
    }

    CompletenessData {
        contig_name: contig_name.to_string(),
        prevalence_left: left_result.map(|(p, _, _)| p * 100.0), // Store as percentage
        distance_left: left_result.map(|(_, d, _)| d),
        min_missing_left: left_result.map(|(_, _, m)| m),
        prevalence_right: right_result.map(|(p, _, _)| p * 100.0), // Store as percentage
        distance_right: right_result.map(|(_, d, _)| d),
        min_missing_right: right_result.map(|(_, _, m)| m),
        total_mismatches: if total_mismatches > 0.0 { Some(total_mismatches) } else { None },
        total_deletions: if total_deletions > 0.0 { Some(total_deletions) } else { None },
        total_insertions: if total_insertions > 0.0 { Some(total_insertions) } else { None },
        total_reads_clipped: if total_reads_clipped > 0.0 { Some(total_reads_clipped) } else { None },
        total_reference_clipped: if total_reference_clipped > 0.0 { Some(total_reference_clipped) } else { None },
    }
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

/// Apply terminal repeat merging for phagetermini analysis.
///
/// For duplications with ≥90% identity (both DTR and ITR):
/// - Merge signals from second region into first region (position i → position i)
/// - Zero out second region
///
/// This ensures termini detection isn't confused by duplicated regions.
fn apply_dtr_merging(
    arrays: &mut FeatureArrays,
    dtr_regions: &[DtrRegion],
) {
    if dtr_regions.is_empty() {
        return;
    }

    // Save debug CSV with position data for ALL positions (before merging)
    if let Ok(mut file) = std::fs::File::create("/tmp/dtr_debug.csv") {
        use std::io::Write;
        let _ = writeln!(file, "position,coverage,coverage_reduced,left_clipping,right_clipping,read_start,read_end,tau");
        let array_len = arrays.primary_reads.len();
        for idx in 0..array_len {
            let pos = idx + 1; // 1-indexed position
            let coverage = arrays.primary_reads[idx] as f64;
            let coverage_reduced = arrays.coverage_reduced.get(idx).copied().unwrap_or(0);
            let left_clipping = arrays.left_clipping_lengths.get(idx).map(|v| v.len()).unwrap_or(0);
            let right_clipping = arrays.right_clipping_lengths.get(idx).map(|v| v.len()).unwrap_or(0);
            let read_start = arrays.reads_starts.get(idx).copied().unwrap_or(0);
            let read_end = arrays.reads_ends.get(idx).copied().unwrap_or(0);
            let tau = if coverage > 0.0 { (read_start + read_end) as f64 / coverage } else { 0.0 };
            let _ = writeln!(file, "{},{},{},{},{},{},{},{:.6}", pos, coverage as u64, coverage_reduced, left_clipping, right_clipping, read_start, read_end, tau);
        }
        eprintln!("DEBUG: Wrote {} positions to /tmp/dtr_debug.csv", array_len);
    }

    // Process each duplication pair using pre-built DtrRegion
    for dtr in dtr_regions {
        // Convert from 1-indexed (DtrRegion) to 0-indexed (array indices)
        let first_start = (dtr.first_start - 1) as usize;
        let first_end = (dtr.first_end - 1) as usize;
        let second_start = (dtr.second_start - 1) as usize;
        let second_end = (dtr.second_end - 1) as usize;

        let region_len = first_end.saturating_sub(first_start) + 1;
        let second_region_len = second_end.saturating_sub(second_start) + 1;

        // Skip if regions are different lengths or out of bounds
        if region_len != second_region_len || second_end >= arrays.coverage_reduced.len() {
            continue;
        }

        // Merge signals from second region to first region
        for i in 0..region_len {
            let second_idx = second_start + i;
            let first_idx = first_start + i;

            // Merge coverage_reduced
            if first_idx < arrays.coverage_reduced.len() && second_idx < arrays.coverage_reduced.len() {
                arrays.coverage_reduced[first_idx] += arrays.coverage_reduced[second_idx];
                arrays.coverage_reduced[second_idx] = 0;
            }

            // Merge reads_starts
            if first_idx < arrays.reads_starts.len() && second_idx < arrays.reads_starts.len() {
                arrays.reads_starts[first_idx] += arrays.reads_starts[second_idx];
                arrays.reads_starts[second_idx] = 0;
            }

            // Merge reads_ends
            if first_idx < arrays.reads_ends.len() && second_idx < arrays.reads_ends.len() {
                arrays.reads_ends[first_idx] += arrays.reads_ends[second_idx];
                arrays.reads_ends[second_idx] = 0;
            }

            // Left clippings: keep only if present at both
            if first_idx < arrays.left_clipping_lengths.len() && second_idx < arrays.left_clipping_lengths.len() {
                let has_first = !arrays.left_clipping_lengths[first_idx].is_empty();
                let has_second = !arrays.left_clipping_lengths[second_idx].is_empty();
                if has_first && has_second {
                    let second_clips: Vec<u32> = arrays.left_clipping_lengths[second_idx].drain(..).collect();
                    arrays.left_clipping_lengths[first_idx].extend(second_clips);
                } else {
                    arrays.left_clipping_lengths[first_idx].clear();
                    arrays.left_clipping_lengths[second_idx].clear();
                }
            }

            // Right clippings: keep only if present at both
            if first_idx < arrays.right_clipping_lengths.len() && second_idx < arrays.right_clipping_lengths.len() {
                let has_first = !arrays.right_clipping_lengths[first_idx].is_empty();
                let has_second = !arrays.right_clipping_lengths[second_idx].is_empty();
                if has_first && has_second {
                    let second_clips: Vec<u32> = arrays.right_clipping_lengths[second_idx].drain(..).collect();
                    arrays.right_clipping_lengths[first_idx].extend(second_clips);
                } else {
                    arrays.right_clipping_lengths[first_idx].clear();
                    arrays.right_clipping_lengths[second_idx].clear();
                }
            }
        }
    }
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

    // Apply DTR merging if phagetermini is enabled and there are valid terminal duplications
    if flags.phagetermini && !dtr_regions.is_empty() {
        apply_dtr_merging(arrays, &dtr_regions);
    }

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

        // MAPQ - average mapping quality per position (sum_mapq / primary_reads)
        let mapq_f64: Vec<f64> = arrays.sum_mapq.iter()
            .zip(&arrays.primary_reads)
            .map(|(&sum, &count)| if count > 0 { sum as f64 / count as f64 } else { 0.0 })
            .collect();
        add_compressed_feature(&mapq_f64, "mapq", contig_name, config, output);
    }

    // Assemblycheck features
    let mut left_clip_pos: Option<Vec<u32>> = None;
    let mut right_clip_pos: Option<Vec<u32>> = None;
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
        left_clip_pos = Some(left_clip_runs.iter().map(|r| r.start_pos as u32).collect());

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
        right_clip_pos = Some(right_clip_runs.iter().map(|r| r.start_pos as u32).collect());
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
            let (mut arrays, coverage_pct) = match process_contig_streaming(&mut bam, &ref_name, ref_length, config.sequencing_type, flags, config.circular, config.min_coverage) {
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
            let (packaging_info, completeness_info) = add_features_from_arrays(&mut arrays, &ref_name, ref_length, config, flags, duplications, &mut features);

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
