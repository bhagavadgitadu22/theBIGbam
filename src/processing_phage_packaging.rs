//! Phage packaging mechanism classification and terminal repeat detection.

use crate::compress::Run;
use crate::db::RepeatsData;
use crate::features::FeatureArrays;

/// Peak for phage termini detection.
#[derive(Clone, Debug)]
pub struct Peak {
    pub position: i32,
    pub value: f64,
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
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_aligned_fraction: 90,
            min_identity_dtr: 90,
            max_distance_duplication: 100,
            min_events: 10,
            min_frequency: 10,
            max_distance_peaks: 20
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
/// Returns Some(clip_pos) if a nearby clip is found, None otherwise.
fn has_nearby_clip(
    peak_pos: i32,
    clip_positions: &[u32],
    min_distance: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> Option<i32> {
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
            return Some(clip_pos);
        }
    }
    None
}

pub fn filter_peaks_by_clippings(
    reads_starts: &[Run],
    reads_ends: &[Run],
    left_clip_runs: &[Run],
    right_clip_runs: &[Run],
    min_events: f64,
    min_distance: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (Vec<Run>, Vec<Run>) {
    // Filter clipping runs by min_events threshold, then extract positions
    let left_clips: Vec<u32> = left_clip_runs
        .iter()
        .filter(|r| r.value as f64 >= min_events)
        .map(|r| r.start_pos as u32)
        .collect();
    let right_clips: Vec<u32> = right_clip_runs
        .iter()
        .filter(|r| r.value as f64 >= min_events)
        .map(|r| r.start_pos as u32)
        .collect();

    // Filter reads_starts: keep peaks above threshold without nearby left clipping
    let filtered_starts: Vec<Run> = reads_starts
        .iter()
        .filter(|run| {
            (run.value as f64) >= min_events
                && has_nearby_clip(run.start_pos, &left_clips, min_distance, genome_length, circular, dtr_regions).is_none()
        })
        .cloned()
        .collect();

    // Filter reads_ends: keep peaks above threshold without nearby right clipping
    let filtered_ends: Vec<Run> = reads_ends
        .iter()
        .filter(|run| {
            (run.value as f64) >= min_events
                && has_nearby_clip(run.start_pos, &right_clips, min_distance, genome_length, circular, dtr_regions).is_none()
        })
        .cloned()
        .collect();

    (filtered_starts, filtered_ends)
}

/// Merge nearby peaks within max_distance, keeping the highest value peak.
/// Considers DTR-equivalent positions when calculating distances.
pub fn merge_peaks(
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

/// Classify phage packaging mechanism based on peak configuration.
/// Returns (mechanism, left_termini, right_termini, duplication_status).
pub fn classify_packaging(
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

            // Check for ITR configuration: COS-like peaks in ITR regions
            // For ITR, we use the repeat size from autoblast, not the peak distance
            let itr_config = if distance <= 20 {
                check_itr_configuration(start_pos, end_pos, dtr_regions)
            } else {
                None
            };

            let (mechanism, left_pos, right_pos) = if let Some(repeat_size) = itr_config {
                // ITR detected - classify based on repeat size from autoblast
                let suffix = if end_before_start { "_3'" } else { "_5'" };
                let (left, right) = if end_before_start {
                    (end_pos, start_pos)
                } else {
                    (start_pos, end_pos)
                };

                let itr_type = if repeat_size <= 1000 {
                    format!("ITR_short{}", suffix)
                } else if repeat_size <= genome_10pct {
                    format!("ITR_long{}", suffix)
                } else {
                    format!("ITR_outlier{}", suffix)
                };
                (itr_type, left, right)
            } else if distance < 2 {
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

/// Apply terminal repeat merging for phagetermini analysis.
///
/// For duplications with ≥90% identity (both DTR and ITR):
/// - Merge signals from second region into first region (position i → position i)
/// - Zero out second region
///
/// This ensures termini detection isn't confused by duplicated regions.
pub fn apply_dtr_merging(
    arrays: &mut FeatureArrays,
    dtr_regions: &[DtrRegion],
) {
    if dtr_regions.is_empty() {
        return;
    }

    // Process each duplication pair using pre-built DtrRegion
    // Only merge DTR (direct repeats), not ITR (inverted repeats)
    for dtr in dtr_regions {
        // Skip ITR - merging doesn't make sense for inverted repeats
        if !dtr.is_direct {
            continue;
        }

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

        // Merge signals from second region to first region (DTR only)
        for i in 0..region_len {
            let second_idx = second_start + i;
            let first_idx = first_start + i;

            // Merge coverage_primary_reads
            if first_idx < arrays.primary_reads.len() && second_idx < arrays.primary_reads.len() {
                arrays.primary_reads[first_idx] += arrays.primary_reads[second_idx];
                arrays.primary_reads[second_idx] = 0;
            }

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
