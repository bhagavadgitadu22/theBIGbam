//! Adaptive Run-Length Encoding (RLE) compression for genomic feature data.
//!
//! Compresses signals by grouping consecutive positions with similar values into runs.
//! - **Curves**: Uses adaptive RLE with self-referential threshold (compress_ratio * run_value).
//! - **Bars**: Filters positions where `value > coverage * compress_ratio`. Consecutive
//!   significant positions are grouped into runs with average values.

use crate::processing::ProcessConfig;
use crate::types::{get_plot_type, PlotType};

// ============================================================================
// RUN STRUCTURE
// ============================================================================

/// A run of consecutive positions with similar values.
#[derive(Debug, Clone)]
pub struct Run {
    /// First position in the run (1-indexed for database compatibility)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// The value for this run
    pub value: f32,
    /// Optional relative value (value / coverage * 1000) for coverage-dependent features
    pub value_relative: Option<f32>,
}

// ============================================================================
// MAIN COMPRESSION FUNCTION
// ============================================================================

/// Filter positions by coverage threshold and compute relative values.
///
/// Returns Vec of (position, value, value_relative) tuples for positions where
/// value > coverage * threshold.
///
/// # Arguments
/// * `values` - The signal values (e.g., mismatches, non_inward_pairs)
/// * `coverage` - The coverage reference signal
/// * `threshold` - Fraction threshold (e.g., 0.1 for 10%)
///
/// # Returns
/// Vector of (0-indexed position, absolute value, relative value * 1000)
fn filter_by_coverage(values: &[f64], coverage: &[f64], threshold: f64) -> Vec<(usize, f32, f32)> {
    let n = values.len().min(coverage.len());
    let mut filtered = Vec::new();
    
    for i in 0..n {
        let val = values[i];
        let cov = coverage[i];
        let thresh = cov * threshold;
        
        if val > thresh {
            let value_relative = if cov > 0.0 {
                ((val / cov) * 1000.0) as f32
            } else {
                0.0
            };
            filtered.push((i, val as f32, value_relative));
        }
    }
    
    filtered
}

/// Compress signal using adaptive RLE with optional coverage reference.
///
/// # Parameters
/// - `values`: The signal to compress (e.g., coverage at each position)
/// - `reference`: Optional reference signal (e.g., coverage for filtering). Use None for standard curves.
/// - `plot_type`: Curve (with/without filtering) or Bars (threshold filtering)
/// - `curve_ratio`: Relative tolerance for adaptive RLE curves (e.g., 10 for 10% range threshold)
/// - `bar_ratio`: Threshold for coverage-dependent filtering (e.g., 10 for 10%)
///
/// # Algorithm
/// - **Bars with reference**: Coverage-filtered single-position spikes. Each position where
///   `value > coverage * bar_ratio` is saved as a separate single-position run.
/// - **Curves with reference**: Coverage-filtered then exact-value RLE. Filter positions where
///   `value > coverage * bar_ratio`, then merge consecutive positions with identical values.
/// - **Curves without reference**: Standard adaptive RLE. Range-based grouping with
///   `(max - min) <= ratio * min(|min|, |max|)` threshold.
///
/// # Example (Bars)
/// ```
/// coverage =     [1000, 1000, 1000, 50,  50,  1000, 1000]
/// mismatches =   [5,    5,    5,    5,   5,   5,    5]
/// bar_ratio = 0.1
///
/// Position 1-3: 5 <= 0.1*1000 → filtered (insignificant)
/// Position 4-5: 5 > 0.1*50 → saved as runs (significant at low coverage)
/// Position 6-7: 5 <= 0.1*1000 → filtered (insignificant)
/// ```
pub fn compress_signal_with_reference(
    values: &[f64],
    reference: Option<&[f64]>,
    plot_type: PlotType,
    mut curve_ratio: f64,
    mut bar_ratio: f64,
) -> Vec<Run> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }

    let mut runs = Vec::new();
    curve_ratio *= 0.01; // Convert percentage to fraction
    bar_ratio *= 0.01;   // Convert percentage to fraction

    // Case 1: Bars with reference - coverage-filtered single-position spikes
    if matches!(plot_type, PlotType::Bars) && reference.is_some() {
        let coverage = reference.unwrap();
        let filtered = filter_by_coverage(values, coverage, bar_ratio);
        
        for (pos, val, val_rel) in filtered {
            runs.push(Run {
                start_pos: (pos + 1) as i32,
                end_pos: (pos + 1) as i32,
                value: val,
                value_relative: Some(val_rel),
            });
        }
        return runs;
    }

    // Case 2: Curves with reference - coverage-filtered then exact-value RLE
    if matches!(plot_type, PlotType::Curve) && reference.is_some() {
        let coverage = reference.unwrap();
        let filtered = filter_by_coverage(values, coverage, bar_ratio);
        
        if filtered.is_empty() {
            return runs;
        }

        // Start first run
        let mut run_start_idx = 0;
        let mut current_value = filtered[0].1;
        let mut rel_sum = filtered[0].2 as f64;
        let mut rel_count = 1;

        for i in 1..filtered.len() {
            let (pos, val, val_rel) = filtered[i];
            let prev_pos = filtered[i - 1].0;

            // Check if consecutive positions with same value
            if pos == prev_pos + 1 && (val - current_value).abs() < 1e-6 {
                // Extend run
                rel_sum += val_rel as f64;
                rel_count += 1;
            } else {
                // Close current run
                let run_start_pos = filtered[run_start_idx].0;
                let run_end_pos = filtered[i - 1].0;
                runs.push(Run {
                    start_pos: (run_start_pos + 1) as i32,
                    end_pos: (run_end_pos + 1) as i32,
                    value: current_value,
                    value_relative: Some((rel_sum / rel_count as f64) as f32),
                });

                // Start new run
                run_start_idx = i;
                current_value = val;
                rel_sum = val_rel as f64;
                rel_count = 1;
            }
        }

        // Save last run
        let run_start_pos = filtered[run_start_idx].0;
        let run_end_pos = filtered[filtered.len() - 1].0;
        runs.push(Run {
            start_pos: (run_start_pos + 1) as i32,
            end_pos: (run_end_pos + 1) as i32,
            value: current_value,
            value_relative: Some((rel_sum / rel_count as f64) as f32),
        });

        return runs;
    }

    // Case 3: Curves without reference - standard adaptive RLE
    // Track min/max within each run to prevent drift on gradual changes
    let mut run_start = 0;
    let mut run_sum = values[0];
    let mut run_count = 1;
    let mut run_min = values[0];
    let mut run_max = values[0];

    for i in 1..n {
        let val = values[i];
        let new_min = run_min.min(val);
        let new_max = run_max.max(val);
        let range = new_max - new_min;

        // Threshold relative to the smaller absolute value in the run
        let min_abs = new_min.abs().min(new_max.abs());
        let threshold = if min_abs < 1e-9 {
            curve_ratio // absolute threshold for near-zero values
        } else {
            curve_ratio * min_abs
        };

        if range <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_min = new_min;
            run_max = new_max;
        } else {
            // Close current run (store average as value)
            runs.push(Run {
                start_pos: (run_start + 1) as i32,
                end_pos: (run_start + run_count) as i32,
                value: (run_sum / run_count as f64) as f32,
                value_relative: None,
            });

            // Start new run
            run_start = i;
            run_sum = val;
            run_count = 1;
            run_min = val;
            run_max = val;
        }
    }

    // Save the last run
    runs.push(Run {
        start_pos: (run_start + 1) as i32,
        end_pos: (run_start + run_count) as i32,
        value: (run_sum / run_count as f64) as f32,
        value_relative: None,
    });

    runs
}

/// Merge consecutive runs with identical values (0% tolerance RLE).
///
/// Applied to features that should be constant along entire reads (deletions, mate flags).
/// Only merges runs that are BOTH adjacent (end+1 == start) AND have the same value.
#[inline]
pub fn merge_identical_runs(runs: Vec<Run>) -> Vec<Run> {
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


/// Compress a feature signal using adaptive RLE.
///
/// Convenience wrapper around `add_compressed_feature_with_reference` with no reference signal.
#[inline]
pub fn add_compressed_feature(
    values: &[f64],
    feature: &str,
    _contig_name: &str,
    config: &ProcessConfig,
) {
    add_compressed_feature_with_reference(values, None, feature, config);
}

/// Compress a feature signal with optional coverage reference.
///
/// Returns the compressed runs (used downstream for metrics like clipping runs).
#[inline]
pub fn add_compressed_feature_with_reference(
    values: &[f64],
    reference: Option<&[f64]>,
    feature: &str,
    config: &ProcessConfig,
) -> Vec<Run> {
    let plot_type = get_plot_type(feature);
    compress_signal_with_reference(values, reference, plot_type, config.curve_ratio, config.bar_ratio)
}

/// Compress a feature signal with median statistics (for reads_starts/reads_ends).
///
/// Returns the compressed runs.
#[inline]
pub fn add_compressed_feature_with_median(
    counts: &[f64],
    _medians: &[f64],
    reference: Option<&[f64]>,
    feature: &str,
    config: &ProcessConfig,
) -> Vec<Run> {
    let plot_type = get_plot_type(feature);
    compress_signal_with_reference(counts, reference, plot_type, config.curve_ratio, config.bar_ratio)
}

/// Compress a feature signal with statistics (for clippings/insertions).
///
/// Returns the compressed runs (used for clipping run metrics).
#[inline]
pub fn add_compressed_feature_with_stats(
    counts: &[f64],
    _means: &[f64],
    _medians: &[f64],
    _stds: &[f64],
    reference: Option<&[f64]>,
    feature: &str,
    config: &ProcessConfig,
) -> Vec<Run> {
    let plot_type = get_plot_type(feature);
    compress_signal_with_reference(counts, reference, plot_type, config.curve_ratio, config.bar_ratio)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adaptive_rle_simple() {
        // Test case from the documentation example
        let values = vec![100.0, 102.0, 98.0, 200.0, 205.0, 95.0, 100.0];
        let curve_ratio = 10.0; // 10% tolerance
        let bar_ratio = 10.0; // 10% tolerance
        
        let runs = compress_signal_with_reference(&values, None, PlotType::Curve, curve_ratio, bar_ratio);
        
        // Expected runs (approximately):
        // Run 1: positions 1-3, value ~100
        // Run 2: positions 4-5, value ~202.5
        // Run 3: positions 6-7, value ~97.5
        assert_eq!(runs.len(), 3);
        
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 3);
        assert!((runs[0].value - 100.0).abs() < 5.0); // ~100
        
        assert_eq!(runs[1].start_pos, 4);
        assert_eq!(runs[1].end_pos, 5);
        assert!((runs[1].value - 202.5).abs() < 5.0); // ~202.5
        
        assert_eq!(runs[2].start_pos, 6);
        assert_eq!(runs[2].end_pos, 7);
        assert!((runs[2].value - 97.5).abs() < 5.0); // ~97.5
    }

    #[test]
    fn test_adaptive_rle_constant() {
        // All values the same - should produce one run
        let values = vec![50.0; 1000];
        let curve_ratio = 10.0;
        let bar_ratio = 10.0;
        
        let runs = compress_signal_with_reference(&values, None, PlotType::Curve, curve_ratio, bar_ratio);
        
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 1000);
        assert_eq!(runs[0].value, 50.0);
    }

    #[test]
    fn test_adaptive_rle_spike() {
        // Constant with a spike - spike should be its own run
        let values = vec![100.0, 100.0, 500.0, 100.0, 100.0];
        let curve_ratio = 10.0;
        let bar_ratio = 10.0;
        
        let runs = compress_signal_with_reference(&values, None, PlotType::Curve, curve_ratio, bar_ratio);
        
        // Should have 3 runs: [1-2], [3], [4-5]
        assert_eq!(runs.len(), 3);
        
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 2);
        
        assert_eq!(runs[1].start_pos, 3);
        assert_eq!(runs[1].end_pos, 3); // Singleton run
        assert_eq!(runs[1].value, 500.0);
        
        assert_eq!(runs[2].start_pos, 4);
        assert_eq!(runs[2].end_pos, 5);
    }

    #[test]
    fn test_adaptive_rle_empty() {
        let values: Vec<f64> = vec![];
        let curve_ratio = 10.0;
        let bar_ratio = 10.0;
        
        let runs = compress_signal_with_reference(&values, None, PlotType::Curve, curve_ratio, bar_ratio);
        assert_eq!(runs.len(), 0);
    }

    #[test]
    fn test_context_aware_bars() {
        // Test bar threshold filtering with coverage reference
        // Scenario: constant mismatches (5) but varying coverage
        let coverage = vec![1000.0, 1000.0, 1000.0, 50.0, 50.0, 1000.0, 1000.0];
        let mismatches = vec![5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0];
        let curve_ratio = 10.0; // 10% threshold
        let bar_ratio = 10.0;
        
        // Positions 0-2: 5 <= 0.1*1000 (100) → filtered out (insignificant at high coverage)
        // Positions 3-4: 5 > 0.1*50 (5) → saved as run (significant at low coverage!)
        // Positions 5-6: 5 <= 0.1*1000 (100) → filtered out (back to insignificant)
        
        let runs = compress_signal_with_reference(&mismatches, Some(&coverage), PlotType::Bars, curve_ratio, bar_ratio);
        
        // Should have only 1 run: the low-coverage spike (others filtered out)
        assert_eq!(runs.len(), 1);
        
        assert_eq!(runs[0].start_pos, 4);
        assert_eq!(runs[0].end_pos, 5);
        assert_eq!(runs[0].value, 5.0);
    }

    #[test]
    fn test_adaptive_rle_gradual_decline() {
        // Gradual decline from 1000 to 500 in steps of 10
        // Old algorithm: single run (running average drifts along)
        // New algorithm: must break into multiple runs (range exceeds 10%)
        let values: Vec<f64> = (0..=50).map(|i| 1000.0 - i as f64 * 10.0).collect();
        let runs = compress_signal_with_reference(&values, None, PlotType::Curve, 10.0, 10.0);
        // Must NOT be a single run (that would mean infinite drift)
        assert!(runs.len() > 1, "gradual decline must break into multiple runs, got {}", runs.len());
    }

    #[test]
    fn test_bars_without_reference_fallback() {
        // Test that bars without reference fall back to curve-like behavior
        let values = vec![5.0, 5.0, 50.0, 50.0, 5.0];
        let curve_ratio = 10.0;
        let bar_ratio = 10.0;
        
        let runs = compress_signal_with_reference(&values, None, PlotType::Bars, curve_ratio, bar_ratio);
        
        // Should create runs based on value changes (10x jump)
        assert!(runs.len() >= 2); // At least spike separated
    }
}
