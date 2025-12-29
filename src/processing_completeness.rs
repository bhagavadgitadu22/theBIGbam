//! Completeness statistics computation for assembly quality assessment.

use crate::compress::Run;
use crate::db::CompletenessData;

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
pub fn compute_completeness(
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