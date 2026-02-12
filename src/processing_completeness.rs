//! Metric computation for assembly quality assessment.
//!
//! Computes Misassembly (≥50% prevalence), Microdiversity (≥10% prevalence),
//! Side misassembly, and Topology data from per-position arrays.

use crate::compress::Run;
use crate::db::{MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData};

/// Compute median from a slice of u32 values. Returns 0 if empty.
fn compute_median(lengths: &[u32]) -> i64 {
    if lengths.is_empty() {
        return 0;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as i64 + sorted[mid] as i64) / 2
    } else {
        sorted[mid] as i64
    }
}

/// Compute median from an i32 slice. Returns 0.0 if empty.
fn median_i32(slice: &[i32]) -> f64 {
    if slice.is_empty() {
        return 0.0;
    }
    let mut sorted = slice.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[mid] as f64
    }
}

/// Compute all metrics from raw per-position arrays in a single pass.
///
/// Returns (MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData).
pub fn compute_all_metrics(
    left_clipping_lengths: &[Vec<u32>],
    right_clipping_lengths: &[Vec<u32>],
    insertion_lengths: &[Vec<u32>],
    deletion_lengths: &[Vec<u32>],
    primary_reads: &[u64],
    mismatches: &[u64],
    deletions: &[u64],
    contig_name: &str,
    contig_length: usize,
    // Circularity inputs
    circularising_reads_count: u64,
    circularising_inserts_count: u64,
    circularising_insert_sizes: &[i32],
    all_proper_insert_sizes: &[i32],
    contig_end_mates_mapped_on_another_contig: u64,
    // Side misassembly: needs clip runs for left/right side detection
    left_clip_runs: &[Run],
    right_clip_runs: &[Run],
) -> (MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData) {
    // === Initialize accumulators ===
    let mut mis_mm: i64 = 0;
    let mut mis_del: i64 = 0;
    let mut mis_ins: i64 = 0;
    let mut mis_clip: i64 = 0;
    let mut collapse_bp: i64 = 0;
    let mut expansion_bp: i64 = 0;

    let mut micro_mm: i64 = 0;
    let mut micro_del: i64 = 0;
    let mut micro_ins: i64 = 0;
    let mut micro_clip: i64 = 0;
    let mut micro_ref_bp: i64 = 0;
    let mut micro_reads_bp: i64 = 0;

    // === Core loop: iterate all positions once ===
    for pos in 0..contig_length {
        let cov = primary_reads.get(pos).copied().unwrap_or(0) as f64;
        if cov == 0.0 {
            continue;
        }

        // Mismatch prevalence
        let mm_count = mismatches.get(pos).copied().unwrap_or(0) as f64;
        let mm_prev = mm_count / cov;

        // Deletion prevalence
        let del_count = deletions.get(pos).copied().unwrap_or(0) as f64;
        let del_prev = del_count / cov;
        let del_median = if pos < deletion_lengths.len() {
            compute_median(&deletion_lengths[pos])
        } else {
            0
        };

        // Insertion prevalence
        let ins_count = if pos < insertion_lengths.len() {
            insertion_lengths[pos].len() as f64
        } else {
            0.0
        };
        let ins_prev = ins_count / cov;
        let ins_median = if pos < insertion_lengths.len() {
            compute_median(&insertion_lengths[pos])
        } else {
            0
        };

        // Clipping prevalence (left + right combined)
        let lclip = if pos < left_clipping_lengths.len() {
            left_clipping_lengths[pos].len() as f64
        } else {
            0.0
        };
        let rclip = if pos < right_clipping_lengths.len() {
            right_clipping_lengths[pos].len() as f64
        } else {
            0.0
        };
        let clip_count = lclip + rclip;
        let clip_prev = clip_count / cov;

        // Combined median for clips
        let clip_median = if (lclip + rclip) > 0.0 {
            let mut combined: Vec<u32> = Vec::new();
            if pos < left_clipping_lengths.len() {
                combined.extend_from_slice(&left_clipping_lengths[pos]);
            }
            if pos < right_clipping_lengths.len() {
                combined.extend_from_slice(&right_clipping_lengths[pos]);
            }
            compute_median(&combined)
        } else {
            0
        };

        // === Misassembly (≥50%) ===
        if mm_prev >= 0.5 {
            mis_mm += 1;
            collapse_bp += 1; // mismatch = 1bp
            expansion_bp += 1;
        }
        if del_prev >= 0.5 {
            mis_del += 1;
            expansion_bp += del_median;
        }
        if ins_prev >= 0.5 {
            mis_ins += 1;
            collapse_bp += ins_median;
        }
        if clip_prev >= 0.5 {
            mis_clip += 1;
            collapse_bp += clip_median;
        }

        // === Microdiversity (≥10%) ===
        if mm_prev >= 0.1 {
            micro_mm += 1;
            micro_reads_bp += 1; // mismatch = 1bp in reads
            micro_ref_bp += 1; // mismatch = 1bp in reference
        }
        if del_prev >= 0.1 {
            micro_del += 1;
            micro_ref_bp += del_median;
        }
        if ins_prev >= 0.1 {
            micro_ins += 1;
            micro_reads_bp += ins_median;
        }
        if clip_prev >= 0.1 {
            micro_clip += 1;
            micro_reads_bp += clip_median;
        }
    }

    // === Expansion_bp paired clips: Σ(distance) for pairs where both have prevalence ≥50% ===
    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum ClipType { Right, Left }

    let mut all_clips: Vec<(i32, ClipType, f64)> = Vec::new();

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

    all_clips.sort_by_key(|(pos, _, _)| *pos);

    // Find consecutive pairs: right-clipping followed by left-clipping, both ≥50%
    for i in 0..all_clips.len().saturating_sub(1) {
        let (pos_right, type_right, prev_right) = all_clips[i];
        let (pos_left, type_left, prev_left) = all_clips[i + 1];

        if type_right == ClipType::Right && type_left == ClipType::Left
            && prev_right >= 0.5 && prev_left >= 0.5
        {
            let distance = pos_left - pos_right;
            if distance > 0 {
                expansion_bp += distance as i64;
            }
        }
    }

    // === Microdiverse_bp_on_reference paired clips: Σ(distance) for pairs where both have prevalence ≥10% ===
    for i in 0..all_clips.len().saturating_sub(1) {
        let (pos_right, type_right, prev_right) = all_clips[i];
        let (pos_left, type_left, prev_left) = all_clips[i + 1];

        if type_right == ClipType::Right && type_left == ClipType::Left
            && prev_right >= 0.1 && prev_left >= 0.1
        {
            let distance = pos_left - pos_right;
            if distance > 0 {
                micro_ref_bp += distance as i64;
            }
        }
    }

    // === Side misassembly ===
    // Find left-side clipping events where median clipped length > distance to start AND prevalence ≥ 50%
    let left_result = left_clip_runs
        .iter()
        .filter_map(|run| {
            let pos = (run.start_pos - 1) as usize;
            let distance_to_start = run.start_pos - 1;
            let min_missing = if pos < left_clipping_lengths.len() {
                compute_median(&left_clipping_lengths[pos]) as i32
            } else {
                0
            };
            if min_missing > distance_to_start {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                if prevalence >= 0.5 {
                    Some((prevalence, distance_to_start, min_missing))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let right_result = right_clip_runs
        .iter()
        .filter_map(|run| {
            let pos = (run.start_pos - 1) as usize;
            let distance_to_end = contig_length as i32 - run.start_pos;
            let min_missing = if pos < right_clipping_lengths.len() {
                compute_median(&right_clipping_lengths[pos]) as i32
            } else {
                0
            };
            if min_missing > distance_to_end {
                let clipping_count = run.value as f64;
                let coverage = primary_reads.get(pos).copied().unwrap_or(0) as f64;
                let prevalence = if coverage > 0.0 { clipping_count / coverage } else { 0.0 };
                if prevalence >= 0.5 {
                    Some((prevalence, distance_to_end, min_missing))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let contig_end_misjoint = if contig_end_mates_mapped_on_another_contig > 0 {
        Some(contig_end_mates_mapped_on_another_contig)
    } else {
        None
    };

    let side = SideMisassemblyData {
        contig_name: contig_name.to_string(),
        coverage_first_position: primary_reads.first().copied().unwrap_or(0),
        contig_start_collapse_percentage: left_result.map(|(p, _, _)| (p * 1000.0).round() as i32),
        contig_start_collapse_bp: left_result.map(|(_, _, m)| m),
        contig_start_expansion_bp: left_result.map(|(_, d, _)| d),
        coverage_last_position: primary_reads.last().copied().unwrap_or(0),
        contig_end_collapse_percentage: right_result.map(|(p, _, _)| (p * 1000.0).round() as i32),
        contig_end_collapse_bp: right_result.map(|(_, _, m)| m),
        contig_end_expansion_bp: right_result.map(|(_, d, _)| d),
        contig_end_misjoint_mates: contig_end_misjoint,
    };

    // === Topology ===
    // Circularising reads percentage based on mean coverage at junction
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

    let circularising_inserts = if circularising_inserts_count > 0 {
        Some(circularising_inserts_count)
    } else {
        None
    };

    // Circularising insert size deviation: median(circ) - median(all_proper)
    let circ_insert_deviation = if !circularising_insert_sizes.is_empty() && !all_proper_insert_sizes.is_empty() {
        let circ_median = median_i32(circularising_insert_sizes);
        let all_median = median_i32(all_proper_insert_sizes);
        Some((circ_median - all_median).round() as i32)
    } else {
        None
    };

    let topology = TopologyData {
        contig_name: contig_name.to_string(),
        circularising_reads,
        circularising_reads_percentage,
        circularising_inserts,
        circularising_insert_size_deviation: circ_insert_deviation,
    };

    let misassembly = MisassemblyData {
        contig_name: contig_name.to_string(),
        mismatches_count: mis_mm,
        deletions_count: mis_del,
        insertions_count: mis_ins,
        clippings_count: mis_clip,
        collapse_bp,
        expansion_bp,
    };

    let microdiversity = MicrodiversityData {
        contig_name: contig_name.to_string(),
        mismatches_count: micro_mm,
        deletions_count: micro_del,
        insertions_count: micro_ins,
        clippings_count: micro_clip,
        microdiverse_bp_on_reference: micro_ref_bp,
        microdiverse_bp_on_reads: micro_reads_bp,
    };

    (misassembly, microdiversity, side, topology)
}
