//! Phage packaging mechanism classification and terminal repeat detection.

use crate::db::RepeatsData;
use crate::types::TerminusArea;
use statrs::distribution::{Poisson, DiscreteCDF, Normal, ContinuousCDF};

/// Peak area for phage termini detection.
/// Represents a merged region of nearby positions with significant signal.
#[derive(Clone, Debug)]
pub struct PeakArea {
    /// First position of the merged area
    pub start_pos: i32,
    /// Last position of the merged area
    pub end_pos: i32,
    /// Position with highest SPC (read starts + read ends)
    pub center_pos: i32,
    /// Total SPC (sum of read_starts + read_ends in area)
    pub total_spc: u32,
    /// Coverage at center position
    pub coverage: u32,
    /// Tau value at center position
    pub tau: f64,
    /// Number of positions merged into this area
    pub number_peaks: u32,
    /// Sum of pre-filtered clippings used for filtering
    pub sum_clippings: u64,
    /// Poisson test p-value (peak significance)
    pub pvalue: f64,
    /// Bonferroni-adjusted p-value
    pub adjusted_pvalue: f64,
    /// Whether this area passed the Poisson significance test
    pub passed_poisson_test: bool,
    /// Whether this area passed the clipping excess test
    pub passed_clipping_test: bool,
    /// Global clipped ratio for this contig/sample
    pub clipped_ratio: f64,
    /// Expected clippings (threshold from statistical test)
    pub expected_clippings: f64,
    /// Expected SPC from Poisson test (lambda_w)
    pub expected_spc: f64,
    /// Size of the area (distance between start_pos and end_pos + 1)
    pub size: i32,
    /// Median clipping length within this area (for diagnostics)
    pub median_clippings: f64,
}

/// Configuration for phage termini detection.
#[derive(Clone, Copy)]
pub struct PhageTerminiConfig {
    /// Minimum aligned fraction to evaluate the phage termini
    pub min_aligned_fraction: i32,
    /// Minimum identity (%) for considering duplications as DTR/ITR
    pub min_identity_dtr: i32,
    /// Maximum distance (bp) from reference ends for duplication regions to be considered valid.
    /// Duplications must have one region starting within this distance from position 1
    /// and the other region ending within this distance from the contig end.
    pub max_distance_duplication: i32,
    /// Minimum clipping length to count a clip as "real" for phagetermini.
    /// Reads with clips shorter than this are treated as clean (near-match).
    pub min_clipping_length: u32,
    /// Minimum event count to consider a signal significant (applies to both peaks and clippings)
    pub min_events: i32,
    /// Minimum frequence of events to consider a signal significant (applies to both peaks and clippings)
    pub min_frequency: f64,
    /// Maximum distance (bp) between peaks to merge them
    pub max_distance_peaks: i32,
    /// Significance threshold for Poisson and clipping test (p-value threshold)
    /// Lower values are more stringent (require higher confidence that peak is real terminus)
    pub threshold_significance: f64,
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_aligned_fraction: 90,
            min_identity_dtr: 90,
            max_distance_duplication: 100,
            min_clipping_length: 5,
            min_events: 10,
            min_frequency: 0.1, // 10% of coverage
            max_distance_peaks: 20,
            threshold_significance: 0.05, // 95% confidence interval
        }
    }
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
pub struct DtrRegion {
    pub first_start: i32,
    pub first_end: i32,
    pub second_start: i32,
    pub second_end: i32,
    /// true = DTR (direct), false = ITR (inverted)
    pub is_direct: bool,
}

/// Statistical test for excess clippings.
/// Returns (passes_test, expected_clippings).
fn is_area_clipping_acceptable(
    total_spc: u64,
    sum_clippings: u64,
    clipped_ratio: f64,
    threshold_significance: f64,
) -> (bool, f64) {
    let expected_clippings = clipped_ratio * total_spc as f64;

    if sum_clippings == 0 {
        return (true, expected_clippings);
    }

    if (sum_clippings as f64) <= expected_clippings {
        return (true, expected_clippings);
    }

    let std_err = expected_clippings.sqrt();
    if std_err == 0.0 {
        if sum_clippings > 0 {
            return (false, expected_clippings);
        }
        return (true, expected_clippings);
    }

    let z = (sum_clippings as f64 - expected_clippings) / std_err;

    let normal = Normal::new(0.0, 1.0).unwrap();
    let z_critical = normal.inverse_cdf(1.0 - threshold_significance);

    if z > z_critical {
        (false, expected_clippings)
    } else {
        (true, expected_clippings)
    }
}

/// Sum pre-filtered clippings within an area and max_distance_peaks of its edges.
/// Uses simple circular/linear distance - DTR deduplication happens at classification.
fn sum_prefiltered_clippings_for_area(
    area: &PeakArea,
    clippings: &[u64],           // already pre-filtered (zeros = insignificant)
    max_distance_peaks: i32,
    genome_length: usize,
    circular: bool,
) -> u64 {
    let mut sum = 0u64;

    for (idx, &clips) in clippings.iter().enumerate() {
        if clips == 0 {
            continue;
        }

        let pos = (idx + 1) as i32; // 1-indexed

        // Check if position is within area bounds
        let in_area = pos >= area.start_pos && pos <= area.end_pos;

        // Check if position is within max_distance_peaks of area edges
        let dist_to_start = if circular {
            circular_distance(pos, area.start_pos, genome_length)
        } else {
            (pos - area.start_pos).abs()
        };
        let dist_to_end = if circular {
            circular_distance(pos, area.end_pos, genome_length)
        } else {
            (pos - area.end_pos).abs()
        };

        let near_area = dist_to_start <= max_distance_peaks || dist_to_end <= max_distance_peaks;

        if in_area || near_area {
            sum += clips;
        }
    }

    sum
}

/// Poisson significance test for merged peak areas.
///
/// For each window W, computes:
/// - λ_W = sum(coverage_reduced[i] * p_rate) for positions i in W
/// - p-value = P(X >= observed | Poisson(λ_W))  (survival function)
/// - adjusted p-value = p-value * K (Bonferroni, K = number of areas tested)
///
/// Areas with adjusted_pvalue > 0.05 have `passed_poisson_test` set to false.
fn poisson_test_areas(
    areas: &mut [PeakArea],
    spc_array: &[u64],           // reads_starts or reads_ends
    coverage_reduced: &[u64],
    genome_length: usize,
    circular: bool,
    threshold_significance: f64
) {
    if areas.is_empty() {
        return;
    }

    // Compute background rate: p = sum(spc) / sum(coverage_reduced)
    let total_spc: f64 = spc_array.iter().sum::<u64>() as f64;
    let total_cov: f64 = coverage_reduced.iter().sum::<u64>() as f64;
    if total_cov == 0.0 || total_spc == 0.0 {
        return;
    }
    let p_rate = total_spc / total_cov;

    let k = areas.len() as f64; // number of windows for Bonferroni

    for area in areas.iter_mut() {
        // Compute size from area positions
        let dist = if circular {
            circular_distance(area.start_pos, area.end_pos, genome_length)
        } else {
            (area.start_pos - area.end_pos).abs()
        };
        area.size = dist + 1;

        // Compute expected λ_W = sum(coverage_reduced[i] * p_rate) for positions in window
        let lambda_w = if circular {
            // In circular mode, use the minimal-distance path between start_pos and end_pos
            let lo = area.start_pos.min(area.end_pos);
            let hi = area.start_pos.max(area.end_pos);
            let lo_idx = (lo - 1) as usize;
            let hi_idx = (hi - 1) as usize;

            // Path A (direct): lo_idx..=hi_idx — spans (hi - lo + 1) positions
            let direct_span = hi_idx - lo_idx + 1;
            // Path B (wrap): 0..=lo_idx + hi_idx..genome_length — spans (genome_length - direct_span + 2) positions
            let wrap_span = genome_length - direct_span + 2;

            if direct_span <= wrap_span {
                (lo_idx..=hi_idx.min(genome_length - 1))
                    .map(|i| coverage_reduced.get(i).copied().unwrap_or(0) as f64 * p_rate)
                    .sum::<f64>()
            } else {
                let part1: f64 = (0..=lo_idx)
                    .map(|i| coverage_reduced.get(i).copied().unwrap_or(0) as f64 * p_rate)
                    .sum();
                let part2: f64 = (hi_idx..genome_length)
                    .map(|i| coverage_reduced.get(i).copied().unwrap_or(0) as f64 * p_rate)
                    .sum();
                part1 + part2
            }
        } else {
            // Linear: use min/max so order doesn't matter
            let start_idx = (area.start_pos.min(area.end_pos) - 1) as usize;
            let end_idx = (area.start_pos.max(area.end_pos) - 1) as usize;
            (start_idx..=end_idx.min(genome_length - 1))
                .map(|i| coverage_reduced.get(i).copied().unwrap_or(0) as f64 * p_rate)
                .sum::<f64>()
        };

        // Store expected SPC (lambda_w) on the area
        area.expected_spc = lambda_w;

        if lambda_w <= 0.0 {
            area.pvalue = 1.0;
            area.adjusted_pvalue = 1.0;
            area.passed_poisson_test = false;
            continue;
        }

        let observed = area.total_spc as u64;

        // P(X >= observed) = 1 - P(X <= observed - 1) = 1 - CDF(observed - 1)
        // For observed == 0, p-value = 1.0
        let pvalue = if observed == 0 {
            1.0
        } else if let Ok(poisson) = Poisson::new(lambda_w) {
            1.0 - poisson.cdf(observed - 1)
        } else {
            1.0
        };

        let adjusted = (pvalue * k).min(1.0);
        area.pvalue = pvalue;
        area.adjusted_pvalue = adjusted;

        if adjusted > threshold_significance {
            area.passed_poisson_test = false;
        }
    }
}

/// Filter positions and merge nearby ones into peak areas, returning ALL areas with metadata.
///
/// This implements:
/// 1. Position Filtering: SPC >= 0.1*coverage_reduced AND SPC >= min_events
/// 2. Peak Merging: merge nearby positions into areas (DTR-aware)
/// 3. Clipping Pre-filter: identify significant clipping positions
/// 4. Area Clipping Aggregation: sum pre-filtered clippings for each area
/// 5. Statistical Test: filter areas with excess clippings
///
/// Returns ALL areas (both passed and discarded).
/// Each area has `passed_poisson_test` and `passed_clipping_test` set to indicate filtering results.
/// Sets filtering metadata on each area (sum_clippings, expected_clippings, clipped_ratio)
pub fn filter_and_merge_to_areas_with_diagnostics(
    reads_starts: &[u64],
    reads_ends: &[u64],
    left_clippings: &[u64],
    right_clippings: &[u64],
    coverage_reduced: &[u64],
    primary_reads: &[u64],
    clipped_ratio: f64,
    config: &PhageTerminiConfig,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (Vec<PeakArea>, Vec<PeakArea>) {
    let min_frequency = config.min_frequency;
    let min_events = config.min_events as u64;
    let max_distance_peaks = config.max_distance_peaks;
    let threshold_significance = config.threshold_significance;

    // Step 1: Position Filtering for starts (local criteria only)
    let mut filtered_start_positions: Vec<(i32, u64)> = Vec::new();
    for (idx, &spc) in reads_starts.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32;
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);

        let local_threshold = (cov as f64 * min_frequency) as u64;
        if spc < local_threshold || spc < min_events {
            continue;
        }

        filtered_start_positions.push((pos, spc));
    }

    // Step 1: Position Filtering for ends (local criteria only)
    let mut filtered_end_positions: Vec<(i32, u64)> = Vec::new();
    for (idx, &spc) in reads_ends.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32;
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);

        let local_threshold = (cov as f64 * min_frequency) as u64;
        if spc < local_threshold || spc < min_events {
            continue;
        }

        filtered_end_positions.push((pos, spc));
    }

    // Step 2: Peak Merging into areas
    let mut start_areas = merge_positions_to_areas(
        &filtered_start_positions,
        reads_starts,
        reads_ends,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
    );

    let mut end_areas = merge_positions_to_areas(
        &filtered_end_positions,
        reads_starts,
        reads_ends,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
    );

    // Step 2b: Poisson peak significance test
    poisson_test_areas(&mut start_areas, reads_starts, coverage_reduced, genome_length, circular, threshold_significance);
    poisson_test_areas(&mut end_areas, reads_ends, coverage_reduced, genome_length, circular, threshold_significance);

    // Step 3: Pre-filter clippings (zero out non-significant positions)
    let mut left_clippings_filtered: Vec<u64> = left_clippings.to_vec();
    for (idx, clips) in left_clippings_filtered.iter_mut().enumerate() {
        if *clips < min_events {
            *clips = 0;
            continue;
        }
        let prim = primary_reads.get(idx).copied().unwrap_or(0);
        let threshold = (prim as f64 * min_frequency) as u64;
        if *clips < threshold {
            *clips = 0;
        }
    }

    let mut right_clippings_filtered: Vec<u64> = right_clippings.to_vec();
    for (idx, clips) in right_clippings_filtered.iter_mut().enumerate() {
        if *clips < min_events {
            *clips = 0;
            continue;
        }
        let prim = primary_reads.get(idx).copied().unwrap_or(0);
        let threshold = (prim as f64 * min_frequency) as u64;
        if *clips < threshold {
            *clips = 0;
        }
    }

    // Step 3b: DTR both-copies confirmation
    // In a circular assembly with DTRs, boundary clippings are artifacts (reads
    // can't extend past the contig edge). Real biological clippings appear at
    // both DTR copies; artifacts appear at only one. Zero any clipping that
    // lacks a counterpart at the equivalent position in the other copy.
    for dtr in dtr_regions {
        if !dtr.is_direct {
            continue; // Only DTR, not ITR
        }

        let gl = genome_length as i32;

        // Snapshot to avoid cascade: zeroing copy A shouldn't cause copy B
        // to fail the check too.
        let left_snapshot = left_clippings_filtered.clone();
        let right_snapshot = right_clippings_filtered.clone();

        // Check first region against second
        for p in dtr.first_start..=dtr.first_end {
            let q = dtr.second_start + (p - dtr.first_start);
            if q < 1 || q > gl { continue; }
            let pi = (p - 1) as usize;
            let qi = (q - 1) as usize;
            if left_snapshot[pi] > 0 && left_snapshot[qi] == 0 {
                left_clippings_filtered[pi] = 0;
            }
            if right_snapshot[pi] > 0 && right_snapshot[qi] == 0 {
                right_clippings_filtered[pi] = 0;
            }
        }

        // Check second region against first
        for p in dtr.second_start..=dtr.second_end {
            let q = dtr.first_start + (p - dtr.second_start);
            if q < 1 || q > gl { continue; }
            let pi = (p - 1) as usize;
            let qi = (q - 1) as usize;
            if left_snapshot[pi] > 0 && left_snapshot[qi] == 0 {
                left_clippings_filtered[pi] = 0;
            }
            if right_snapshot[pi] > 0 && right_snapshot[qi] == 0 {
                right_clippings_filtered[pi] = 0;
            }
        }
    }

    // Step 4-5: Area Left Clipping Aggregation + Statistical Test for starts
    // Update each area with filtering metadata and kept status
    for area in &mut start_areas {
        let sum_clips = sum_prefiltered_clippings_for_area(
            area, &left_clippings_filtered,
            max_distance_peaks, genome_length, circular
        );
        let (passes, expected_clippings) = is_area_clipping_acceptable(
            area.total_spc as u64, sum_clips, clipped_ratio, threshold_significance
        );

        area.sum_clippings = sum_clips;
        area.clipped_ratio = clipped_ratio;
        area.expected_clippings = expected_clippings;
        area.passed_clipping_test = passes;
    }

    // Step 4-5: Area Right Clipping Aggregation + Statistical Test for ends
    for area in &mut end_areas {
        let sum_clips = sum_prefiltered_clippings_for_area(
            area, &right_clippings_filtered,
            max_distance_peaks, genome_length, circular
        );
        let (passes, expected_clippings) = is_area_clipping_acceptable(
            area.total_spc as u64, sum_clips, clipped_ratio, threshold_significance
        );

        area.sum_clippings = sum_clips;
        area.clipped_ratio = clipped_ratio;
        area.expected_clippings = expected_clippings;
        area.passed_clipping_test = passes;
    }

    (start_areas, end_areas)
}

/// Merge filtered positions into peak areas based on proximity.
fn merge_positions_to_areas(
    filtered_positions: &[(i32, u64)], // (pos, spc)
    reads_starts: &[u64],
    reads_ends: &[u64],
    coverage_reduced: &[u64],
    max_distance_peaks: i32,
    genome_length: usize,
    circular: bool,
) -> Vec<PeakArea> {
    if filtered_positions.is_empty() {
        return Vec::new();
    }

    // Sort by position
    let mut positions: Vec<(i32, u64)> = filtered_positions.to_vec();
    positions.sort_by_key(|(pos, _)| *pos);

    let mut areas: Vec<PeakArea> = Vec::new();

    // Start first area
    let (first_pos, first_spc) = positions[0];
    let mut current_start = first_pos;
    let mut current_end = first_pos;
    let mut current_center = first_pos;
    let mut current_max_spc = first_spc;
    let mut current_total_spc = first_spc;
    let mut current_number_peaks: u32 = 1; // Count of positions merged

    for &(pos, spc) in positions.iter().skip(1) {
        let dist = if circular {
            circular_distance(current_end, pos, genome_length)
        } else {
            (current_end - pos).abs()
        };

        if dist <= max_distance_peaks {
            // Extend current area
            current_end = pos;
            current_total_spc += spc;
            current_number_peaks += 1;
            if spc > current_max_spc {
                current_max_spc = spc;
                current_center = pos;
            }
        } else {
            // Finalize current area and start new one
            let idx = (current_center - 1) as usize;
            let cov = coverage_reduced.get(idx).copied().unwrap_or(0) as u32;
            let starts = reads_starts.get(idx).copied().unwrap_or(0);
            let ends = reads_ends.get(idx).copied().unwrap_or(0);
            let tau = if cov > 0 {
                (starts + ends) as f64 / cov as f64
            } else {
                0.0
            };

            // Store positions
            areas.push(PeakArea {
                start_pos: current_start,
                end_pos: current_end,
                center_pos: current_center,
                total_spc: current_total_spc as u32,
                coverage: cov,
                tau,
                number_peaks: current_number_peaks,
                sum_clippings: 0,      // Will be set later in filter step
                pvalue: 1.0,           // Will be set by poisson_test_areas
                adjusted_pvalue: 1.0,  // Will be set by poisson_test_areas
                passed_poisson_test: true,   // Will be set by poisson_test_areas
                passed_clipping_test: true,  // Will be set later in filter step
                clipped_ratio: 0.0,    // Will be set later in filter step
                expected_clippings: 0.0, // Will be set later in filter step
                expected_spc: 0.0,     // Will be set by poisson_test_areas
                size: 0,               // Will be set by poisson_test_areas
                median_clippings: 0.0, // Will be set later
            });

            // Start new area
            current_start = pos;
            current_end = pos;
            current_center = pos;
            current_max_spc = spc;
            current_total_spc = spc;
            current_number_peaks = 1;
        }
    }

    // Finalize last area
    let idx = (current_center - 1) as usize;
    let cov = coverage_reduced.get(idx).copied().unwrap_or(0) as u32;
    let starts = reads_starts.get(idx).copied().unwrap_or(0);
    let ends = reads_ends.get(idx).copied().unwrap_or(0);
    let tau = if cov > 0 {
        (starts + ends) as f64 / cov as f64
    } else {
        0.0
    };

    // Store positions (will be canonicalized before return)
    areas.push(PeakArea {
        start_pos: current_start,
        end_pos: current_end,
        center_pos: current_center,
        total_spc: current_total_spc as u32,
        coverage: cov,
        tau,
        number_peaks: current_number_peaks,
        sum_clippings: 0,      // Will be set later in filter step
        pvalue: 1.0,           // Will be set by poisson_test_areas
        adjusted_pvalue: 1.0,  // Will be set by poisson_test_areas
        passed_poisson_test: true,   // Will be set by poisson_test_areas
        passed_clipping_test: true,  // Will be set later in filter step
        clipped_ratio: 0.0,    // Will be set later in filter step
        expected_clippings: 0.0, // Will be set later in filter step
        expected_spc: 0.0,     // Will be set by poisson_test_areas
        size: 0,               // Will be set by poisson_test_areas
        median_clippings: 0.0, // Will be set later
    });

    // Handle circular wrap-around merge
    if circular && areas.len() > 1 {
        // Copy values we need before mutating
        let first_end_pos = areas[0].end_pos;
        let first_total_spc = areas[0].total_spc;
        let first_center_pos = areas[0].center_pos;
        let first_coverage = areas[0].coverage;
        let first_tau = areas[0].tau;
        let first_start_pos = areas[0].start_pos;
        let first_number_peaks = areas[0].number_peaks;

        let last_idx = areas.len() - 1;
        let last_end_pos = areas[last_idx].end_pos;
        let last_total_spc = areas[last_idx].total_spc;
        let last_center_pos = areas[last_idx].center_pos;
        let last_coverage = areas[last_idx].coverage;
        let last_tau = areas[last_idx].tau;
        let last_start_pos = areas[last_idx].start_pos;
        let last_number_peaks = areas[last_idx].number_peaks;

        let wrap_dist = circular_distance(last_end_pos, first_start_pos, genome_length);

        if wrap_dist <= max_distance_peaks {
            // Merge first and last areas
            let merged_total_spc = first_total_spc + last_total_spc;
            let merged_number_peaks = first_number_peaks + last_number_peaks;
            let (merged_center, merged_cov, merged_tau) = if first_total_spc >= last_total_spc {
                (first_center_pos, first_coverage, first_tau)
            } else {
                (last_center_pos, last_coverage, last_tau)
            };

            // Remove last area
            areas.pop();

            // Merge wrap-around: start from last area
            areas[0] = PeakArea {
                start_pos: last_start_pos,
                end_pos: first_end_pos,
                center_pos: merged_center,
                total_spc: merged_total_spc,
                coverage: merged_cov,
                tau: merged_tau,
                number_peaks: merged_number_peaks,
                sum_clippings: 0,
                pvalue: 1.0,
                adjusted_pvalue: 1.0,
                passed_poisson_test: true,
                passed_clipping_test: true,
                clipped_ratio: 0.0,
                expected_clippings: 0.0,
                expected_spc: 0.0,
                size: 0,
                median_clippings: 0.0,
            };
        }
    }

    areas
}

/// Check if a repeat is valid for phage termini analysis.
/// A valid repeat must have:
/// - One region starting within max_distance from position 1 (beginning)
/// - The other region ending within max_distance from contig_length (end)
pub fn is_valid_terminal_repeat(dup: &RepeatsData, contig_length: usize, max_distance: i32) -> bool {
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


/// Translate a position from second DTR region to first DTR region.
/// Only applies to DTR (direct), not ITR (inverted).
/// Returns the translated position (or original if not in a DTR second region).
fn translate_to_first_dtr_region(pos: i32, dtr_regions: &[DtrRegion]) -> i32 {
    for dtr in dtr_regions {
        // Only translate for DTR (direct), not ITR (inverted)
        if !dtr.is_direct {
            continue;
        }
        // Check if pos is in second region
        if pos >= dtr.second_start && pos <= dtr.second_end {
            return dtr.first_start + pos - dtr.second_start;
        }
    }
    // Not in any DTR second region, return as-is
    pos
}

/// Check if positions are in ITR regions: one in first, other in second.
/// Returns Some(repeat_size) if this is an ITR configuration, None otherwise.
fn check_itr_configuration(pos1: i32, pos2: i32, dtr_regions: &[DtrRegion]) -> Option<i32> {
    for dtr in dtr_regions {
        // Only check ITR regions (is_direct = false)
        if dtr.is_direct {
            continue;
        }

        let in_first_1 = pos1 >= dtr.first_start && pos1 <= dtr.first_end;
        let in_second_1 = pos1 >= dtr.second_start && pos1 <= dtr.second_end;
        let in_first_2 = pos2 >= dtr.first_start && pos2 <= dtr.first_end;
        let in_second_2 = pos2 >= dtr.second_start && pos2 <= dtr.second_end;

        // One position in first ITR region, other in second ITR region
        if (in_first_1 && in_second_2) || (in_second_1 && in_first_2) {
            let repeat_size = dtr.first_end - dtr.first_start + 1;
            return Some(repeat_size);
        }
    }
    None
}

/// Combine duplication status from start and end peaks.
/// Returns Some(true) if all in DTR, Some(false) if all in ITR, None otherwise.
fn combine_duplication_status(
    starts_status: Option<bool>,
    ends_status: Option<bool>,
    has_starts: bool,
    has_ends: bool,
) -> Option<bool> {
    match (starts_status, ends_status, has_starts, has_ends) {
        // Both have same status
        (Some(true), Some(true), _, _) => Some(true),   // All DTR
        (Some(false), Some(false), _, _) => Some(false), // All ITR
        // Only one side has peaks - use that status
        (s, _, true, false) => s,
        (_, s, false, true) => s,
        // Mixed or no peaks
        _ => None,
    }
}

/// Deduplicate peak areas that are DTR-equivalent and determine duplication status.
/// Groups areas by equivalence (using center_pos) and keeps the canonical (first region) position.
/// Only applies to DTR (direct repeats), not ITR (inverted repeats).
///
/// Returns (unique_areas, duplication_status):
/// - unique_areas: areas with center_pos translated to first region
/// - duplication_status: Some(true) if all in DTR, Some(false) if all in ITR, None otherwise
fn deduplicate_dtr_equivalent_areas(areas: &[PeakArea], dtr_regions: &[DtrRegion]) -> (Vec<PeakArea>, Option<bool>) {
    if areas.is_empty() {
        return (Vec::new(), None);
    }

    if dtr_regions.is_empty() {
        return (areas.to_vec(), None);
    }

    let mut unique: Vec<PeakArea> = Vec::new();
    let mut all_in_dtr = true;
    let mut all_in_itr = true;
    let mut any_in_repeat = false;

    for area in areas {
        let pos = area.center_pos;

        // Check which region type this area's center is in
        let mut in_dtr = false;
        let mut in_itr = false;

        for dtr in dtr_regions {
            let in_first = pos >= dtr.first_start && pos <= dtr.first_end;
            let in_second = pos >= dtr.second_start && pos <= dtr.second_end;

            if in_first || in_second {
                any_in_repeat = true;
                if dtr.is_direct {
                    in_dtr = true;
                } else {
                    in_itr = true;
                }
            }
        }

        if !in_dtr {
            all_in_dtr = false;
        }
        if !in_itr {
            all_in_itr = false;
        }

        // Translate center_pos to canonical (first region) position - only for DTR
        let canonical_pos = translate_to_first_dtr_region(pos, dtr_regions);

        // Check if this canonical position already exists (within tolerance)
        let is_duplicate = unique.iter().any(|existing| {
            (existing.center_pos - canonical_pos).abs() <= 20 // Use 20bp tolerance
        });

        if !is_duplicate {
            // Store area with canonical position for correct distance calculation
            let mut translated_area = area.clone();
            translated_area.center_pos = canonical_pos;
            unique.push(translated_area);
        }
    }

    // Determine duplication status
    let status = if !any_in_repeat {
        None
    } else if all_in_dtr && !all_in_itr {
        Some(true)  // All in DTR
    } else if all_in_itr && !all_in_dtr {
        Some(false) // All in ITR
    } else {
        None // Mixed or unclear
    };

    (unique, status)
}

/// Convert PeakArea to TerminusArea for database storage.
fn peak_area_to_terminus_area(area: &PeakArea) -> TerminusArea {
    TerminusArea {
        start_pos: area.start_pos,
        end_pos: area.end_pos,
        center_pos: area.center_pos,
        total_spc: area.total_spc,
        coverage: area.coverage,
        tau: area.tau,
        number_peaks: area.number_peaks,
        pvalue: area.pvalue,
        adjusted_pvalue: area.adjusted_pvalue,
        passed_poisson_test: area.passed_poisson_test,
        passed_clipping_test: area.passed_clipping_test,
        sum_clippings: area.sum_clippings,
        clipped_ratio: area.clipped_ratio,
        expected_clippings: area.expected_clippings,
        expected_spc: area.expected_spc,
        size: area.size,
        median_clippings: area.median_clippings,
    }
}

/// Classify phage packaging mechanism based on peak areas.
/// Returns (mechanism, all_left_termini, all_right_termini, duplication_status, repeat_length).
///
/// DTR equivalence is handled here: areas from both DTR regions are received,
/// and we deduplicate equivalent areas to determine the "real" area count
/// for classification.
///
/// Classification uses only KEPT areas (those that passed filtering),
/// but ALL areas (kept and discarded) are returned for database storage.
///
/// `repeat_length` is the distance between start and end peak centers when there is
/// exactly 1 unique start and 1 unique end area, None otherwise.
pub fn classify_packaging_areas(
    start_areas: &[PeakArea],
    end_areas: &[PeakArea],
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (String, Vec<TerminusArea>, Vec<TerminusArea>, Option<bool>, Option<i32>, Vec<f64>, Vec<f64>) {
    // Filter to only kept areas for classification
    let kept_start_areas: Vec<&PeakArea> = start_areas.iter().filter(|a| a.passed_poisson_test && a.passed_clipping_test).collect();
    let kept_end_areas: Vec<&PeakArea> = end_areas.iter().filter(|a| a.passed_poisson_test && a.passed_clipping_test).collect();

    // Deduplicate DTR-equivalent areas from KEPT areas only
    let kept_start_clones: Vec<PeakArea> = kept_start_areas.iter().map(|&a| a.clone()).collect();
    let kept_end_clones: Vec<PeakArea> = kept_end_areas.iter().map(|&a| a.clone()).collect();
    let (unique_starts, starts_status) = deduplicate_dtr_equivalent_areas(&kept_start_clones, dtr_regions);
    let (unique_ends, ends_status) = deduplicate_dtr_equivalent_areas(&kept_end_clones, dtr_regions);

    // Use unique counts from KEPT areas for classification logic
    let mut repeat_length: Option<i32> = None;
    let mechanism = match (unique_starts.len(), unique_ends.len()) {
        (0, 0) => "No_packaging".to_string(),

        (1, 0) | (0, 1) => "PAC".to_string(),

        (1, 1) => {
            // Two unique areas - classify based on distance and order
            // Use center positions for classification
            let start_pos = unique_starts[0].center_pos;
            let end_pos = unique_ends[0].center_pos;

            // Calculate actual genomic distance (shortest path in circular)
            let distance = if circular {
                circular_distance(start_pos, end_pos, genome_length)
            } else {
                (start_pos - end_pos).abs()
            };

            // Store repeat_length as the distance between start and end peak centers
            repeat_length = Some(distance);

            let genome_10pct = (genome_length as f64 * 0.1) as i32;

            // Determine order: which comes first going forward (clockwise in circular)
            let end_before_start = if circular {
                let gl = genome_length as i32;
                let forward_dist = if end_pos >= start_pos {
                    end_pos - start_pos
                } else {
                    (gl - start_pos) + end_pos
                };
                forward_dist != distance
            } else {
                end_pos < start_pos
            };

            // Check for ITR configuration
            let itr_config = if distance <= 20 {
                check_itr_configuration(start_pos, end_pos, dtr_regions)
            } else {
                None
            };

            let suffix = if end_before_start { "_3'" } else { "_5'" };

            if let Some(repeat_size) = itr_config {
                if repeat_size <= 1000 {
                    format!("ITR_short{}", suffix)
                } else if repeat_size <= genome_10pct {
                    format!("ITR_long{}", suffix)
                } else {
                    format!("ITR_outlier{}", suffix)
                }
            } else if distance < 2 {
                "COS".to_string()
            } else if distance <= 20 {
                format!("COS{}", suffix)
            } else if distance <= 1000 {
                format!("DTR_short{}", suffix)
            } else if distance <= genome_10pct {
                format!("DTR_long{}", suffix)
            } else {
                format!("DTR_outlier{}", suffix)
            }
        }

        _ => "Unknown_packaging".to_string(),
    };

    // Return ALL areas (both kept and discarded) for database storage
    let all_left: Vec<TerminusArea> = start_areas.iter().map(peak_area_to_terminus_area).collect();
    let all_right: Vec<TerminusArea> = end_areas.iter().map(peak_area_to_terminus_area).collect();
    let dup_status = combine_duplication_status(starts_status, ends_status, !unique_starts.is_empty(), !unique_ends.is_empty());

    // Collect median clippings from kept termini
    let median_left: Vec<f64> = kept_start_areas.iter().map(|a| a.median_clippings).collect();
    let median_right: Vec<f64> = kept_end_areas.iter().map(|a| a.median_clippings).collect();

    (mechanism, all_left, all_right, dup_status, repeat_length, median_left, median_right)
}

