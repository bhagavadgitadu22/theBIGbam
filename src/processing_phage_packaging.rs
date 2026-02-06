//! Phage packaging mechanism classification and terminal repeat detection.

use crate::db::RepeatsData;
use crate::types::TerminusArea;

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
    /// Total clippings in area (for informational purposes)
    pub total_clips: u32,
    /// Coverage at center position
    pub coverage: u32,
    /// Tau value at center position
    pub tau: f64,
    /// Number of positions merged into this area
    pub number_peaks: u32,
    /// Sum of pre-filtered clippings used for filtering
    pub sum_clippings: u64,
    /// Whether this area passed the clipping test
    pub kept: bool,
    /// Global clipped ratio for this contig/sample
    pub clipped_ratio: f64,
    /// Expected clippings (threshold from statistical test)
    pub expected_clippings: f64,
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
    /// Minimum event count to consider a signal significant (applies to both peaks and clippings)
    pub min_events: i32,
    /// Minimum frequence of events to consider a signal significant (applies to both peaks and clippings)
    pub min_frequency: i32,
    /// Maximum distance (bp) between peaks to merge them
    pub max_distance_peaks: i32,
    /// Significance threshold for clipping test (p-value threshold)
    /// Lower values are more stringent (require higher confidence that peak is real terminus)
    pub clipping_significance: f64,
    /// Decay factor for density-based peak merging (0.0-1.0).
    /// Controls how much density can drop when merging areas.
    pub decay_factor_for_merge: f64,
    /// Minimum SPC to consider a position for merging (loose pre-filter).
    pub min_spc_for_merge: i32,
    /// Minimum frequency (%) of coverage_reduced for merging (loose pre-filter).
    pub min_frequency_for_merge: i32,
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_aligned_fraction: 90,
            min_identity_dtr: 90,
            max_distance_duplication: 100,
            min_events: 10,
            min_frequency: 10,
            max_distance_peaks: 20,
            clipping_significance: 0.05, // 95% confidence interval
            decay_factor_for_merge: 0.8,
            min_spc_for_merge: 3,
            min_frequency_for_merge: 1,
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

/// Micro-area for density-based merging.
/// Tracks both original positions (for final area bounds) and canonical positions (for density calculation).
#[derive(Clone, Debug)]
struct MicroArea {
    /// Original start position (1-indexed)
    start_pos: i32,
    /// Original end position (1-indexed)
    end_pos: i32,
    /// Canonical start position (translated to first DTR region)
    canonical_start: i32,
    /// Canonical end position (translated to first DTR region)
    canonical_end: i32,
    /// Total SPC accumulated in this area
    total_spc: u64,
    /// Position with highest SPC (in original coordinates)
    max_spc_pos: i32,
    /// Highest SPC value
    max_spc: u64,
    /// Number of positions merged into this area
    num_positions: u32,
}

impl MicroArea {
    /// Calculate density using canonical span
    fn density(&self) -> f64 {
        let span = (self.canonical_end - self.canonical_start + 1) as f64;
        if span <= 0.0 {
            return self.total_spc as f64;
        }
        self.total_spc as f64 / span
    }
}

/// Check if two positions are close enough to merge.
/// Close enough means:
/// - Physically close: distance <= max_distance (including circularity), OR
/// - DTR-equivalent close: pos1 is in a DTR region and pos2 is within
///   max_distance of the equivalent position in the other DTR region (or vice versa)
fn dtr_aware_close_enough(
    pos1: i32,
    pos2: i32,
    max_distance: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> bool {
    // Check physical proximity
    let phys_dist = if circular {
        circular_distance(pos1, pos2, genome_length)
    } else {
        (pos1 - pos2).abs()
    };
    if phys_dist <= max_distance {
        return true;
    }

    // Check DTR-equivalent proximity (only for DTR, not ITR)
    for dtr in dtr_regions {
        if !dtr.is_direct {
            continue;
        }

        // If pos1 is in first region, check if pos2 is close to equivalent in second
        if pos1 >= dtr.first_start && pos1 <= dtr.first_end {
            let equiv = dtr.second_start + (pos1 - dtr.first_start);
            let dist = if circular {
                circular_distance(pos2, equiv, genome_length)
            } else {
                (pos2 - equiv).abs()
            };
            if dist <= max_distance {
                return true;
            }
        }

        // If pos1 is in second region, check if pos2 is close to equivalent in first
        if pos1 >= dtr.second_start && pos1 <= dtr.second_end {
            let equiv = dtr.first_start + (pos1 - dtr.second_start);
            let dist = if circular {
                circular_distance(pos2, equiv, genome_length)
            } else {
                (pos2 - equiv).abs()
            };
            if dist <= max_distance {
                return true;
            }
        }

        // If pos2 is in first region, check if pos1 is close to equivalent in second
        if pos2 >= dtr.first_start && pos2 <= dtr.first_end {
            let equiv = dtr.second_start + (pos2 - dtr.first_start);
            let dist = if circular {
                circular_distance(pos1, equiv, genome_length)
            } else {
                (pos1 - equiv).abs()
            };
            if dist <= max_distance {
                return true;
            }
        }

        // If pos2 is in second region, check if pos1 is close to equivalent in first
        if pos2 >= dtr.second_start && pos2 <= dtr.second_end {
            let equiv = dtr.first_start + (pos2 - dtr.second_start);
            let dist = if circular {
                circular_distance(pos1, equiv, genome_length)
            } else {
                (pos1 - equiv).abs()
            };
            if dist <= max_distance {
                return true;
            }
        }
    }

    false
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

/// Statistical test for excess clippings with detailed reason for rejection.
/// Returns (passes_test, expected_clippings, rejection_reason).
fn is_area_clipping_acceptable_with_reason(
    total_spc: u64,
    sum_clippings: u64,
    clipped_ratio: f64,
    clipping_significance: f64,
) -> (bool, f64, String) {
    let expected_clippings = clipped_ratio * total_spc as f64;
    
    if sum_clippings == 0 {
        return (true, expected_clippings, String::new());
    }

    if (sum_clippings as f64) <= expected_clippings {
        return (true, expected_clippings, String::new());
    }

    let std_err = expected_clippings.sqrt();
    if std_err == 0.0 {
        if sum_clippings > 0 {
            return (false, expected_clippings, "excess_clippings: expected=0".to_string());
        }
        return (true, expected_clippings, String::new());
    }

    let z = (sum_clippings as f64 - expected_clippings) / std_err;

    let z_critical = if clipping_significance >= 0.10 {
        1.282
    } else if clipping_significance >= 0.05 {
        1.645
    } else if clipping_significance >= 0.01 {
        2.326
    } else {
        2.576
    };

    if z > z_critical {
        (false, expected_clippings, format!("excess_clippings: z={:.2} > z_critical={:.2}", z, z_critical))
    } else {
        (true, expected_clippings, String::new())
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

/// Density-based merging of positions into areas.
/// 
/// Algorithm:
/// 1. Initialize micro-areas from positions passing loose thresholds
/// 2. Greedy merge: find closest pair (physically OR DTR-equivalent close),
///    merge if combined density >= min(d_A, d_B) * decay_factor
/// 3. Handle circular wrap-around
/// 4. Return merged micro-areas
fn density_based_merge(
    reads_array: &[u64],
    coverage_reduced: &[u64],
    config: &PhageTerminiConfig,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> Vec<MicroArea> {
    let min_spc = config.min_spc_for_merge as u64;
    let min_freq = config.min_frequency_for_merge as f64 / 100.0;
    let max_distance = config.max_distance_peaks;
    let decay_factor = config.decay_factor_for_merge;

    // Step 1: Initialize micro-areas from qualifying positions
    let mut areas: Vec<MicroArea> = Vec::new();
    for (idx, &spc) in reads_array.iter().enumerate() {
        if spc < min_spc {
            continue;
        }
        let pos = (idx + 1) as i32;
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);
        let threshold = (cov as f64 * min_freq) as u64;
        if spc < threshold {
            continue;
        }

        // Compute canonical position for span/density calculation
        let canonical = translate_to_first_dtr_region(pos, dtr_regions);

        areas.push(MicroArea {
            start_pos: pos,
            end_pos: pos,
            canonical_start: canonical,
            canonical_end: canonical,
            total_spc: spc,
            max_spc_pos: pos,
            max_spc: spc,
            num_positions: 1,
        });
    }

    if areas.is_empty() {
        return areas;
    }

    // Sort by canonical position for efficient scanning
    areas.sort_by_key(|a| a.canonical_start);

    // Step 2: Greedy density-based merge
    // Use a simple O(N^2) approach since N is typically small (< 100)
    let mut merged = true;
    while merged {
        merged = false;

        // Find best merge candidate
        let mut best_i: Option<usize> = None;
        let mut best_j: Option<usize> = None;
        let mut best_dist = i32::MAX;

        for i in 0..areas.len() {
            for j in (i + 1)..areas.len() {
                // Check if areas are close enough (physically OR via DTR equivalence)
                let close = dtr_aware_close_enough(
                    areas[i].end_pos,
                    areas[j].start_pos,
                    max_distance,
                    genome_length,
                    circular,
                    dtr_regions,
                );
                if !close {
                    continue;
                }

                // Use canonical distance for prioritization
                let dist = (areas[j].canonical_start - areas[i].canonical_end).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best_i = Some(i);
                    best_j = Some(j);
                }
            }
        }

        // Try to merge best candidate
        if let (Some(i), Some(j)) = (best_i, best_j) {
            let density_i = areas[i].density();
            let density_j = areas[j].density();
            let min_density = density_i.min(density_j);

            // Calculate combined density using canonical span
            let combined_spc = areas[i].total_spc + areas[j].total_spc;
            let combined_canonical_start = areas[i].canonical_start.min(areas[j].canonical_start);
            let combined_canonical_end = areas[i].canonical_end.max(areas[j].canonical_end);
            let combined_span = (combined_canonical_end - combined_canonical_start + 1) as f64;
            let combined_density = combined_spc as f64 / combined_span;

            if combined_density >= min_density * decay_factor {
                // Merge: keep area i, absorb area j
                let area_j = areas.remove(j);

                // Update bounds (use original positions for bounds)
                areas[i].start_pos = areas[i].start_pos.min(area_j.start_pos);
                areas[i].end_pos = areas[i].end_pos.max(area_j.end_pos);

                // Update canonical bounds
                areas[i].canonical_start = combined_canonical_start;
                areas[i].canonical_end = combined_canonical_end;

                // Update SPC
                areas[i].total_spc = combined_spc;
                areas[i].num_positions += area_j.num_positions;

                // Update max SPC position
                if area_j.max_spc > areas[i].max_spc {
                    areas[i].max_spc = area_j.max_spc;
                    areas[i].max_spc_pos = area_j.max_spc_pos;
                }

                merged = true;
            }
        }
    }

    // Step 3: Handle circular wrap-around merge
    if circular && areas.len() > 1 {
        // Check if first and last areas can merge (wrap around genome)
        let first_idx = 0;
        let last_idx = areas.len() - 1;

        let close = dtr_aware_close_enough(
            areas[last_idx].end_pos,
            areas[first_idx].start_pos,
            max_distance,
            genome_length,
            circular,
            dtr_regions,
        );

        if close {
            let density_first = areas[first_idx].density();
            let density_last = areas[last_idx].density();
            let min_density = density_first.min(density_last);

            let combined_spc = areas[first_idx].total_spc + areas[last_idx].total_spc;
            // For wrap-around, span is complex - use sum of individual spans
            let span_first = (areas[first_idx].canonical_end - areas[first_idx].canonical_start + 1) as f64;
            let span_last = (areas[last_idx].canonical_end - areas[last_idx].canonical_start + 1) as f64;
            let combined_span = span_first + span_last;
            let combined_density = combined_spc as f64 / combined_span;

            if combined_density >= min_density * decay_factor {
                // Merge: absorb first into last, then swap to position 0
                let first_area = areas.remove(first_idx);
                let last_idx = areas.len() - 1; // Recalculate after removal

                // Wrap-around merge: start from last area, end at first area
                areas[last_idx].start_pos = areas[last_idx].start_pos; // Keep last start
                areas[last_idx].end_pos = first_area.end_pos; // Use first end

                // Canonical: use minimum start, maximum end
                areas[last_idx].canonical_start = areas[last_idx].canonical_start.min(first_area.canonical_start);
                areas[last_idx].canonical_end = areas[last_idx].canonical_end.max(first_area.canonical_end);

                areas[last_idx].total_spc = combined_spc;
                areas[last_idx].num_positions += first_area.num_positions;

                if first_area.max_spc > areas[last_idx].max_spc {
                    areas[last_idx].max_spc = first_area.max_spc;
                    areas[last_idx].max_spc_pos = first_area.max_spc_pos;
                }
            }
        }
    }

    areas
}

/// Convert MicroAreas to PeakAreas, computing coverage and tau at the center.
fn convert_micro_to_peak_areas(
    micro_areas: Vec<MicroArea>,
    reads_starts: &[u64],
    reads_ends: &[u64],
    coverage_reduced: &[u64],
) -> Vec<PeakArea> {
    micro_areas
        .into_iter()
        .map(|ma| {
            let center_idx = (ma.max_spc_pos - 1) as usize;
            let cov = coverage_reduced.get(center_idx).copied().unwrap_or(0) as u32;
            let starts = reads_starts.get(center_idx).copied().unwrap_or(0);
            let ends = reads_ends.get(center_idx).copied().unwrap_or(0);
            let tau = if cov > 0 {
                (starts + ends) as f64 / cov as f64
            } else {
                0.0
            };

            PeakArea {
                start_pos: ma.start_pos,
                end_pos: ma.end_pos,
                center_pos: ma.max_spc_pos,
                total_spc: ma.total_spc as u32,
                total_clips: 0, // Will be computed later if needed
                coverage: cov,
                tau,
                number_peaks: ma.num_positions,
                sum_clippings: 0,
                kept: true, // Will be set by filtering
                clipped_ratio: 0.0,
                expected_clippings: 0.0,
            }
        })
        .collect()
}

/// Apply strict area-level filters after density-based merge.
/// Sets `kept = false` for areas that don't pass.
fn apply_area_filters(
    areas: &mut [PeakArea],
    coverage_reduced: &[u64],
    min_events: u64,
    min_frequency: f64,
) {
    for area in areas.iter_mut() {
        let center_idx = (area.center_pos - 1) as usize;
        let cov = coverage_reduced.get(center_idx).copied().unwrap_or(0);
        let threshold = (cov as f64 * min_frequency) as u64;

        let spc = area.total_spc as u64;
        if spc < min_events || spc < threshold {
            area.kept = false;
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
/// Returns ALL areas (both kept and discarded)
/// Each area has its `kept` field set to indicate whether it passed filtering
/// Sets filtering metadata on each area (kept, sum_clippings, expected_clippings, clipped_ratio)
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
    let min_frequency = (config.min_frequency as f64) / 100.0;
    let min_events = config.min_events as u64;
    let max_distance_peaks = config.max_distance_peaks;
    let clipping_significance = config.clipping_significance;

    // Step 1-2: Density-based merge (uses loose thresholds internally, then strict area-level filters)
    // This approach merges first, then filters - capturing scattered signals better than pre-filter→merge
    let start_micro = density_based_merge(
        reads_starts,
        coverage_reduced,
        config,
        genome_length,
        circular,
        dtr_regions,
    );
    let end_micro = density_based_merge(
        reads_ends,
        coverage_reduced,
        config,
        genome_length,
        circular,
        dtr_regions,
    );

    // Convert MicroAreas to PeakAreas
    let mut start_areas = convert_micro_to_peak_areas(&start_micro, reads_starts, reads_ends, coverage_reduced);
    let mut end_areas = convert_micro_to_peak_areas(&end_micro, reads_starts, reads_ends, coverage_reduced);

    // Apply strict area-level filters (min_events, min_frequency on total_spc)
    apply_area_filters(&mut start_areas, coverage_reduced, min_events, min_frequency);
    apply_area_filters(&mut end_areas, coverage_reduced, min_events, min_frequency);

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
    // In a doubled assembly with DTRs, boundary clippings are artifacts (reads
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
        let (passes, expected_clippings, _reason) = is_area_clipping_acceptable_with_reason(
            area.total_spc as u64, sum_clips, clipped_ratio, clipping_significance
        );

        area.sum_clippings = sum_clips;
        area.clipped_ratio = clipped_ratio;
        area.expected_clippings = expected_clippings;
        area.kept = passes;
    }

    // Step 4-5: Area Right Clipping Aggregation + Statistical Test for ends
    for area in &mut end_areas {
        let sum_clips = sum_prefiltered_clippings_for_area(
            area, &right_clippings_filtered,
            max_distance_peaks, genome_length, circular
        );
        let (passes, expected_clippings, _reason) = is_area_clipping_acceptable_with_reason(
            area.total_spc as u64, sum_clips, clipped_ratio, clipping_significance
        );

        area.sum_clippings = sum_clips;
        area.clipped_ratio = clipped_ratio;
        area.expected_clippings = expected_clippings;
        area.kept = passes;
    }

    (start_areas, end_areas)
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
        kept: area.kept,
        sum_clippings: area.sum_clippings,
        clipped_ratio: area.clipped_ratio,
        expected_clippings: area.expected_clippings,
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
) -> (String, Vec<TerminusArea>, Vec<TerminusArea>, Option<bool>, Option<i32>) {
    // Filter to only kept areas for classification
    let kept_start_areas: Vec<&PeakArea> = start_areas.iter().filter(|a| a.kept).collect();
    let kept_end_areas: Vec<&PeakArea> = end_areas.iter().filter(|a| a.kept).collect();

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

    (mechanism, all_left, all_right, dup_status, repeat_length)
}

