//! Completeness statistics computation for assembly quality assessment.

use crate::compress::Run;
use crate::db::CompletenessData;

/// Compute completeness statistics from clipping, insertion, mismatch, and deletion data.
///
/// For the left side: finds all significant left-clipping events where clipped_length > distance_to_start,
/// then picks the one with the highest prevalence (clipping_count / coverage).
///
/// For the right side: same logic with right-clipping events where clipped_length > distance_to_end.
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
/// Returns None for a side if no clipping event satisfies the condition.
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
    circularising_reads_count: u64,
    circularising_inserts_count: u64,
    circularising_insert_sizes: &[i32],
    all_proper_insert_sizes: &[i32],
    contig_end_unmapped_mates: u64,
    contig_end_mates_mapped_on_another_contig: u64,
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

    // Find all left-clipping events where median clipped length > distance to start,
    // then pick the one with the highest prevalence (clipping_count / coverage).
    let left_result = left_clip_runs
        .iter()
        .filter_map(|run| {
            let pos = (run.start_pos - 1) as usize;
            let distance_to_start = run.start_pos - 1;
            let min_missing = if pos < left_clipping_lengths.len() {
                compute_median(&left_clipping_lengths[pos])
            } else {
                0
            };
            if min_missing > distance_to_start {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                Some((prevalence, distance_to_start, min_missing))
            } else {
                None
            }
        })
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Find all right-clipping events where median clipped length > distance to end,
    // then pick the one with the highest prevalence (clipping_count / coverage).
    let right_result = right_clip_runs
        .iter()
        .filter_map(|run| {
            let pos = (run.start_pos - 1) as usize;
            let distance_to_end = contig_length as i32 - run.start_pos;
            let min_missing = if pos < right_clipping_lengths.len() {
                compute_median(&right_clipping_lengths[pos])
            } else {
                0
            };
            if min_missing > distance_to_end {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                Some((prevalence, distance_to_end, min_missing))
            } else {
                None
            }
        })
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

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

    // Total reference clipped: Σ(min_prevalence * distance) for paired clips (right followed by left)
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
                let min_prevalence = prev_right.min(prev_left);
                total_reference_clipped += min_prevalence * distance as f64;
            }
        }
    }

    // Compute circularising reads percentage based on mean coverage at junction
    let circularising_reads = if circularising_reads_count > 0 {
        Some(circularising_reads_count)
    } else {
        None
    };

    let circularising_reads_percentage = if circularising_reads_count > 0 && !primary_reads.is_empty() {
        let first = primary_reads[0] as f64;
        let last = primary_reads[primary_reads.len() - 1] as f64;
        let mean_junction_coverage = (first + last) / 2.0;
        if mean_junction_coverage > 0.0 {
            Some(((circularising_reads_count as f64 / mean_junction_coverage) * 100.0).round() as i32)
        } else {
            None
        }
    } else {
        None
    };

    // Compute circularising inserts percentage based on mean coverage at junction
    let circularising_inserts_percentage = if circularising_inserts_count > 0 && !primary_reads.is_empty() {
        let first = primary_reads[0] as f64;
        let last = primary_reads[primary_reads.len() - 1] as f64;
        let mean_junction_coverage = (first + last) / 2.0;
        if mean_junction_coverage > 0.0 {
            Some(((circularising_inserts_count as f64 / mean_junction_coverage) * 100.0).round() as i32)
        } else {
            None
        }
    } else {
        None
    };

    // Compute contig end unmapped mates percentage based on mean coverage at junction
    let contig_end_unmapped_mates_pct = if contig_end_unmapped_mates > 0 && !primary_reads.is_empty() {
        let first = primary_reads[0] as f64;
        let last = primary_reads[primary_reads.len() - 1] as f64;
        let mean_junction_coverage = (first + last) / 2.0;
        if mean_junction_coverage > 0.0 {
            Some(((contig_end_unmapped_mates as f64 / mean_junction_coverage) * 100.0).round() as i32)
        } else {
            None
        }
    } else {
        None
    };

    // Compute contig end mates mapped on another contig percentage based on mean coverage at junction
    let contig_end_mates_another_contig_pct = if contig_end_mates_mapped_on_another_contig > 0 && !primary_reads.is_empty() {
        let first = primary_reads[0] as f64;
        let last = primary_reads[primary_reads.len() - 1] as f64;
        let mean_junction_coverage = (first + last) / 2.0;
        if mean_junction_coverage > 0.0 {
            Some(((contig_end_mates_mapped_on_another_contig as f64 / mean_junction_coverage) * 100.0).round() as i32)
        } else {
            None
        }
    } else {
        None
    };

    // Compute circularising inserts metrics
    let circularising_inserts = if circularising_inserts_count > 0 {
        Some(circularising_inserts_count)
    } else {
        None
    };

    let (mean_extra_insert, median_extra_insert) = if !circularising_insert_sizes.is_empty()
        && !all_proper_insert_sizes.is_empty()
    {
        // Mean extra
        let circ_mean = circularising_insert_sizes.iter().map(|&x| x as f64).sum::<f64>()
            / circularising_insert_sizes.len() as f64;
        let all_mean = all_proper_insert_sizes.iter().map(|&x| x as f64).sum::<f64>()
            / all_proper_insert_sizes.len() as f64;

        // Median extra (sort copies)
        fn median_i32(slice: &[i32]) -> f64 {
            let mut sorted = slice.to_vec();
            sorted.sort_unstable();
            let mid = sorted.len() / 2;
            if sorted.len() % 2 == 0 {
                (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
            } else {
                sorted[mid] as f64
            }
        }
        let circ_median = median_i32(circularising_insert_sizes);
        let all_median = median_i32(all_proper_insert_sizes);

        (
            Some((circ_mean - all_mean).round() as i32),
            Some((circ_median - all_median).round() as i32),
        )
    } else {
        (None, None)
    };

    let unmapped_mates_ends = if contig_end_unmapped_mates > 0 {
        Some(contig_end_unmapped_mates)
    } else {
        None
    };
    let mates_another_contig_ends = if contig_end_mates_mapped_on_another_contig > 0 {
        Some(contig_end_mates_mapped_on_another_contig)
    } else {
        None
    };

    CompletenessData {
        contig_name: contig_name.to_string(),
        prevalence_left: left_result.map(|(p, _, _)| p * 100.0), // Store as percentage
        left_contamination_length: left_result.map(|(_, d, _)| d),
        left_missing_length: left_result.map(|(_, _, m)| m),
        prevalence_right: right_result.map(|(p, _, _)| p * 100.0), // Store as percentage
        right_contamination_length: right_result.map(|(_, d, _)| d),
        right_missing_length: right_result.map(|(_, _, m)| m),
        total_mismatches: if total_mismatches > 0.0 { Some(total_mismatches) } else { None },
        total_deletions: if total_deletions > 0.0 { Some(total_deletions) } else { None },
        total_insertions: if total_insertions > 0.0 { Some(total_insertions) } else { None },
        total_reads_clipped: if total_reads_clipped > 0.0 { Some(total_reads_clipped) } else { None },
        total_reference_clipped: if total_reference_clipped > 0.0 { Some(total_reference_clipped) } else { None },
        circularising_reads,
        circularising_reads_percentage,
        circularising_inserts,
        mean_extra_insert_length: mean_extra_insert,
        median_extra_insert_length: median_extra_insert,
        contig_end_unmapped_mates: unmapped_mates_ends,
        contig_end_mates_mapped_on_another_contig: mates_another_contig_ends,
        circularising_inserts_percentage,
        contig_end_unmapped_mates_percentage: contig_end_unmapped_mates_pct,
        contig_end_mates_mapped_on_another_contig_percentage: contig_end_mates_another_contig_pct,
    }
}