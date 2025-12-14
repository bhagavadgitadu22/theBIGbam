//! Adaptive Run-Length Encoding (RLE) compression for genomic feature data.
//!
//! Compresses signals by grouping consecutive positions with similar values into runs.
//! - **Curves**: Uses adaptive RLE with self-referential threshold (compress_ratio * run_value).
//! - **Bars**: Filters positions where `value > coverage * compress_ratio`. Consecutive 
//!   significant positions are grouped into runs with average values.

use crate::types::PlotType;

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
}

// ============================================================================
// MAIN COMPRESSION FUNCTION
// ============================================================================

/// Compress signal using adaptive run-length encoding.
///
/// # Algorithm
/// 1. Start first run with position 0 and its value
/// 2. For each subsequent position:
///    - If |value - current_run_value| <= compress_ratio * |current_run_value|, extend run
///    - Otherwise, close current run and start new one
/// 3. Return runs as (start_pos, end_pos, value) tuples
///
/// Compress signal using adaptive RLE with optional coverage reference.
///
/// # Parameters
/// - `values`: The signal to compress (e.g., coverage at each position)
/// - `reference`: Optional reference signal (e.g., coverage for bar plots). Use None for curves.
/// - `plot_type`: Curve (self-referential RLE) or Bars (threshold filtering)
/// - `compress_ratio`: Relative tolerance (e.g., 10 for 10% change threshold)
///
/// # Algorithm
/// - **Curves** (reference=None): Adaptive RLE. Compresses based on value changes.
/// - **Bars** (reference=Some): Threshold filtering. Saves positions where `value > coverage * compress_ratio`.
///   Consecutive significant positions are grouped into runs with average values.
///
/// # Example
/// ```
/// coverage =     [1000, 1000, 1000, 50,  50,  1000, 1000]
/// mismatches =   [5,    5,    5,    5,   5,   5,    5]
/// compress_ratio = 0.1
///
/// Position 1-3: 5 <= 0.1*1000 → filtered (insignificant)
/// Position 4-5: 5 > 0.1*50 → saved as run (significant at low coverage)
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

    // For bars with reference: save positions where value > coverage * bar_ratio
    curve_ratio *= 0.01; // Convert percentage to fraction
    bar_ratio *= 0.01;   // Convert percentage to fraction
    if matches!(plot_type, PlotType::Bars) && reference.is_some() {
        let coverage = reference.unwrap();
        let mut i = 0;
        while i < n {
            let val = values[i];
            let threshold = coverage[i] * bar_ratio;
            
            if val > threshold {
                // Found a significant position - check if consecutive positions also qualify
                let run_start = i;
                let mut run_sum = val;
                let mut run_count = 1;
                
                // Extend run while consecutive positions are also significant
                while i + 1 < n {
                    let next_val = values[i + 1];
                    let next_threshold = coverage[i + 1] * bar_ratio;
                    if next_val > next_threshold {
                        i += 1;
                        run_sum += next_val;
                        run_count += 1;
                    } else {
                        break;
                    }
                }
                
                // Save the run with average value
                runs.push(Run {
                    start_pos: (run_start + 1) as i32,
                    end_pos: (i + 1) as i32,
                    value: (run_sum / run_count as f64) as f32,
                });
            }
            i += 1;
        }
        return runs;
    }

    // For curves: use adaptive RLE compression
    let mut run_start = 0;
    let mut run_value = values[0];
    let mut run_sum = values[0];
    let mut run_count = 1;

    for i in 1..n {
        let val = values[i];
        
        // Self-referential threshold for curves
        let threshold = if run_value.abs() < 1e-9 {
            curve_ratio
        } else {
            curve_ratio * run_value.abs()
        };

        if (val - run_value).abs() <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_value = run_sum / run_count as f64;
        } else {
            // Close current run and save it
            runs.push(Run {
                start_pos: (run_start + 1) as i32,
                end_pos: (run_start + run_count) as i32,
                value: run_value as f32,
            });

            // Start new run
            run_start = i;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Save the last run
    runs.push(Run {
        start_pos: (run_start + 1) as i32,
        end_pos: (run_start + run_count) as i32,
        value: run_value as f32,
    });

    runs
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
