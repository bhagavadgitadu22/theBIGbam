//! GC content computation for genomic sequences.
//!
//! Computes GC content (percentage of G and C bases) using a sliding window approach
//! and applies RLE compression for efficient storage.

/// Configuration parameters for GC content and GC skew computation.
#[derive(Debug, Clone, Copy)]
pub struct GCParams {
    /// Window size for GC content computation (default 500bp)
    pub gc_content_window_size: usize,
    /// Window size for GC skew computation (default 1000bp)
    pub gc_skew_window_size: usize,
}

impl Default for GCParams {
    fn default() -> Self {
        Self {
            gc_content_window_size: 500,
            gc_skew_window_size: 1000,
        }
    }
}

/// A run of consecutive positions with similar GC percentage.
#[derive(Debug, Clone)]
pub struct GCContentRun {
    /// First position in the run (1-indexed for database compatibility)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// GC percentage for this run (0-100)
    pub gc_percentage: u8,
}

/// Statistics for GC content across a contig.
#[derive(Debug, Clone)]
pub struct GCStats {
    /// Mean GC percentage (0-100)
    pub average: f32,
    /// Standard deviation of windowed GC percentages
    pub sd: f32,
}

/// Statistics for GC skew across a contig.
#[derive(Debug, Clone)]
pub struct GCSkewStats {
    /// Amplitude of GC skew: max(skew) - min(skew)
    /// Measures overall strand asymmetry (range 0 to 2)
    pub amplitude: f32,
    /// Percentage of positions with positive GC skew (more G than C)
    pub percent_positive: f32,
}

/// Run-length encoded GC skew data for storage.
/// GC skew ranges from -1.0 to +1.0, stored as i16 (-100 to +100).
#[derive(Debug, Clone)]
pub struct GCSkewRun {
    /// First position in the run (1-indexed for database compatibility)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// GC skew × 100 (range: -100 to +100)
    pub gc_skew: i16,
}

/// Compute GC content using non-overlapping windows with RLE compression.
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of each non-overlapping window in base pairs (typically 500)
/// - `_contig_variation_percentage`: RLE compression tolerance for contig-level features (default 0.1%)
///
/// # Algorithm
/// 1. Divide the sequence into non-overlapping windows of window_size
/// 2. For each window, compute GC% (excluding N bases)
/// 3. Apply RLE compression: merge consecutive windows with similar GC%
///
/// # Returns
/// Tuple of (runs, stats) where:
/// - runs: Vector of GC content runs with (start_pos, end_pos, gc_percentage)
/// - stats: GCStats with average and sd GC percentages
pub fn compute_gc_content(sequence: &[u8], window_size: usize, contig_variation_percentage: f64) -> (Vec<GCContentRun>, GCStats) {
    let n = sequence.len();
    if n == 0 {
        return (Vec::new(), GCStats { average: 0.0, sd: 0.0 });
    }

    let mut gc_values: Vec<u8> = Vec::new();

    // Compute GC content for each non-overlapping window
    let mut pos = 0;
    while pos < n {
        let window_end = (pos + window_size).min(n);
        let window = &sequence[pos..window_end];

        let (gc_count, valid_count) = count_gc_in_window(window);

        let gc_pct = if valid_count > 0 {
            ((gc_count as f64 / valid_count as f64) * 100.0).round() as u8
        } else {
            50 // Default to 50% if window contains only Ns
        };

        gc_values.push(gc_pct);
        pos = window_end;
    }

    // Compute statistics from raw values before compression
    let stats = compute_gc_stats(&gc_values);

    // Build runs with RLE compression applied to consecutive windows with similar values
    let runs = build_gc_runs(&gc_values, window_size, n, contig_variation_percentage);

    (runs, stats)
}

/// Compute GC statistics (average, sd) from raw GC values.
fn compute_gc_stats(gc_values: &[u8]) -> GCStats {
    if gc_values.is_empty() {
        return GCStats { average: 0.0, sd: 0.0 };
    }

    let n = gc_values.len() as f64;

    // Compute average
    let sum: f64 = gc_values.iter().map(|&v| v as f64).sum();
    let average = sum / n;

    // Compute standard deviation
    let variance: f64 = gc_values.iter()
        .map(|&v| {
            let diff = v as f64 - average;
            diff * diff
        })
        .sum::<f64>() / n;
    let sd = variance.sqrt();

    GCStats {
        average: average as f32,
        sd: sd as f32,
    }
}

/// Compute GC skew using non-overlapping windows with RLE compression.
///
/// GC skew = (G - C) / (G + C) at each position
/// - Range: -1 to +1
/// - Positive = more G than C (leading strand bias)
/// - Negative = more C than G (lagging strand bias)
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of each non-overlapping window in base pairs (typically 1000)
/// - `_contig_variation_percentage`: RLE compression tolerance (default 0.1%)
///
/// # Returns
/// Tuple of (runs, stats) where:
/// - runs: Vector of GC skew runs with (start_pos, end_pos, gc_skew × 100)
/// - stats: GCSkewStats with amplitude and percent_positive
pub fn compute_gc_skew(sequence: &[u8], window_size: usize, contig_variation_percentage: f64) -> (Vec<GCSkewRun>, GCSkewStats) {
    let n = sequence.len();
    if n == 0 {
        return (Vec::new(), GCSkewStats { amplitude: 0.0, percent_positive: 0.0 });
    }

    let mut skew_values: Vec<i16> = Vec::new();

    // Compute GC skew for each non-overlapping window
    let mut pos = 0;
    while pos < n {
        let window_end = (pos + window_size).min(n);
        let window = &sequence[pos..window_end];

        let (g_count, c_count) = count_g_c_in_window(window);

        let skew = if g_count + c_count > 0 {
            (g_count as f64 - c_count as f64) / (g_count as f64 + c_count as f64)
        } else {
            0.0 // Default to 0 if window contains no G or C
        };

        // Store as i16 × 100 (range: -100 to +100)
        skew_values.push((skew * 100.0).round() as i16);
        pos = window_end;
    }

    // Compute statistics from raw values before compression
    if skew_values.is_empty() {
        return (Vec::new(), GCSkewStats { amplitude: 0.0, percent_positive: 0.0 });
    }

    let min_skew = skew_values.iter().cloned().min().unwrap_or(0) as f64 / 100.0;
    let max_skew = skew_values.iter().cloned().max().unwrap_or(0) as f64 / 100.0;
    let amplitude = (max_skew - min_skew) as f32;

    let positive_count = skew_values.iter().filter(|&&s| s > 0).count();
    let percent_positive = (positive_count as f64 / skew_values.len() as f64 * 100.0) as f32;

    let stats = GCSkewStats {
        amplitude,
        percent_positive,
    };

    // Build runs with RLE compression applied to consecutive windows with similar values
    let runs = build_gc_skew_runs(&skew_values, window_size, n, contig_variation_percentage);

    (runs, stats)
}

/// Count G and C bases in a window, excluding N bases.
#[inline]
fn count_gc_in_window(window: &[u8]) -> (usize, usize) {
    let mut gc_count = 0;
    let mut valid_count = 0;

    for &base in window {
        match base {
            b'G' | b'g' | b'C' | b'c' => {
                gc_count += 1;
                valid_count += 1;
            }
            b'A' | b'a' | b'T' | b't' => {
                valid_count += 1;
            }
            // N and other characters are excluded
            _ => {}
        }
    }

    (gc_count, valid_count)
}

/// Count G and C bases separately in a window.
#[inline]
fn count_g_c_in_window(window: &[u8]) -> (usize, usize) {
    let mut g_count = 0;
    let mut c_count = 0;

    for &base in window {
        match base {
            b'G' | b'g' => g_count += 1,
            b'C' | b'c' => c_count += 1,
            _ => {}
        }
    }

    (g_count, c_count)
}

/// Build GC content runs from non-overlapping window values with RLE compression.
/// Applies compression to consecutive windows with similar GC percentages.
fn build_gc_runs(gc_values: &[u8], window_size: usize, total_length: usize, variation_percentage: f64) -> Vec<GCContentRun> {
    if gc_values.is_empty() {
        return Vec::new();
    }

    // Apply RLE compression to consecutive windows with similar values
    let ratio = variation_percentage * 0.01; // Convert percentage to fraction
    let mut runs = Vec::new();

    let mut run_start_idx = 0;
    let mut run_value = gc_values[0] as f64;
    let mut run_sum = gc_values[0] as f64;
    let mut run_count = 1;

    for (idx, &gc_pct) in gc_values.iter().enumerate().skip(1) {
        let val = gc_pct as f64;

        // RLE formula: |x[i] - x[i-1]| <= ratio × min(x[i], x[i-1])
        let min_val = val.min(run_value);
        let threshold = ratio * min_val.max(1.0); // Use at least 1.0 to handle zero values

        if (val - run_value).abs() <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_value = run_sum / run_count as f64;
        } else {
            // Close current run and save it
            let window_start = run_start_idx * window_size;
            let window_end = ((run_start_idx + run_count) * window_size).min(total_length);

            runs.push(GCContentRun {
                start_pos: (window_start + 1) as i32,           // Convert to 1-indexed
                end_pos: window_end as i32,
                gc_percentage: run_value.round() as u8,
            });

            // Start new run
            run_start_idx = idx;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Save the last run
    let window_start = run_start_idx * window_size;
    let window_end = ((run_start_idx + run_count) * window_size).min(total_length);
    runs.push(GCContentRun {
        start_pos: (window_start + 1) as i32,           // Convert to 1-indexed
        end_pos: window_end as i32,
        gc_percentage: run_value.round() as u8,
    });

    runs
}

/// Build GC skew runs from non-overlapping window values with RLE compression.
/// Applies compression to consecutive windows with similar GC skew values.
fn build_gc_skew_runs(skew_values: &[i16], window_size: usize, total_length: usize, variation_percentage: f64) -> Vec<GCSkewRun> {
    if skew_values.is_empty() {
        return Vec::new();
    }

    // Apply RLE compression to consecutive windows with similar values
    let ratio = variation_percentage * 0.01; // Convert percentage to fraction
    let mut runs = Vec::new();

    let mut run_start_idx = 0;
    let mut run_value = skew_values[0] as f64;
    let mut run_sum = skew_values[0] as f64;
    let mut run_count = 1;

    for (idx, &skew) in skew_values.iter().enumerate().skip(1) {
        let val = skew as f64;

        // RLE formula: |x[i] - x[i-1]| <= ratio × |min(|x[i]|, |x[i-1]|)|
        // Use absolute values since skew can be negative
        let min_abs = val.abs().min(run_value.abs());
        let threshold = ratio * min_abs.max(1.0); // Use at least 1.0 to handle zero values

        if (val - run_value).abs() <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_value = run_sum / run_count as f64;
        } else {
            // Close current run and save it
            let window_start = run_start_idx * window_size;
            let window_end = ((run_start_idx + run_count) * window_size).min(total_length);

            runs.push(GCSkewRun {
                start_pos: (window_start + 1) as i32,           // Convert to 1-indexed
                end_pos: window_end as i32,
                gc_skew: run_value.round() as i16,
            });

            // Start new run
            run_start_idx = idx;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Save the last run
    let window_start = run_start_idx * window_size;
    let window_end = ((run_start_idx + run_count) * window_size).min(total_length);
    runs.push(GCSkewRun {
        start_pos: (window_start + 1) as i32,           // Convert to 1-indexed
        end_pos: window_end as i32,
        gc_skew: run_value.round() as i16,
    });

    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_simple() {
        // ATGC sequence - 50% GC
        let sequence = b"ATGCATGCATGC";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with ~50% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert!(run.gc_percentage >= 40 && run.gc_percentage <= 60);
        }
        // Stats should also be around 50%
        assert!(stats.average >= 40.0 && stats.average <= 60.0);
    }

    #[test]
    fn test_gc_content_high_gc() {
        // All GC
        let sequence = b"GGGGCCCCGGGGCCCC";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 100% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 100);
        }
        // Stats should be 100% with 0 sd
        assert_eq!(stats.average, 100.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_low_gc() {
        // All AT
        let sequence = b"AAAATTTTAAAATTTT";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 0% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 0);
        }
        // Stats should be 0% with 0 sd
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_with_n() {
        // Sequence with N bases - N should be excluded
        let sequence = b"ATGCNNNNATGC";
        let (runs, _stats) = compute_gc_content(sequence, 4, 10.0);

        // Should still compute ~50% GC from valid bases
        assert!(!runs.is_empty());
    }

    #[test]
    fn test_gc_content_empty() {
        let sequence: &[u8] = b"";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);
        assert!(runs.is_empty());
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_full_coverage() {
        // Verify that the runs cover the entire sequence
        let sequence = b"ATGCATGCATGCATGCATGC"; // 20 bases
        let (runs, _stats) = compute_gc_content(sequence, 4, 10.0);

        // First run should start at position 1
        assert_eq!(runs.first().unwrap().start_pos, 1);
        // Last run should end at position 20
        assert_eq!(runs.last().unwrap().end_pos, 20);

        // Verify continuous coverage (no gaps)
        for i in 1..runs.len() {
            assert_eq!(runs[i].start_pos, runs[i - 1].end_pos + 1);
        }
    }

    #[test]
    fn test_compression() {
        // Long sequence with uniform GC should compress well
        let sequence: Vec<u8> = b"ATGC".repeat(1000);
        let (runs, stats) = compute_gc_content(&sequence, 100, 10.0);

        // Should compress but still cover full length
        assert!(!runs.is_empty());
        assert_eq!(runs.first().unwrap().start_pos, 1);
        assert_eq!(runs.last().unwrap().end_pos, 4000);
        // Stats should be around 50%
        assert!(stats.average >= 45.0 && stats.average <= 55.0);
    }

    #[test]
    fn test_gc_stats() {
        // Test stats computation with known values
        let gc_values: Vec<u8> = vec![40, 50, 60];
        let stats = compute_gc_stats(&gc_values);

        // Average should be 50
        assert_eq!(stats.average, 50.0);
        // SD should be sqrt((100+0+100)/3) = sqrt(200/3) ≈ 8.16
        assert!(stats.sd > 8.0 && stats.sd < 9.0);
    }

    #[test]
    fn test_gc_skew_balanced() {
        // Equal G and C - skew should be 0 everywhere
        let sequence = b"GCGCGCGCGCGCGCGC";
        let (_runs, skew) = compute_gc_skew(sequence, 4, 10.0);

        // Amplitude should be 0 (or very close)
        assert!(skew.amplitude < 0.1);
        // Percent positive should be around 50% (or 0% if exactly 0)
        // With equal G and C, skew is 0, so no positions are strictly positive
    }

    #[test]
    fn test_gc_skew_g_rich() {
        // All G - skew should be +1 everywhere
        let sequence = b"GGGGGGGGGGGGGGGG";
        let (_runs, skew) = compute_gc_skew(sequence, 4, 10.0);

        // Amplitude should be 0 (constant +1)
        assert_eq!(skew.amplitude, 0.0);
        // Percent positive should be 100%
        assert_eq!(skew.percent_positive, 100.0);
    }

    #[test]
    fn test_gc_skew_c_rich() {
        // All C - skew should be -1 everywhere
        let sequence = b"CCCCCCCCCCCCCCCC";
        let (_runs, skew) = compute_gc_skew(sequence, 4, 10.0);

        // Amplitude should be 0 (constant -1)
        assert_eq!(skew.amplitude, 0.0);
        // Percent positive should be 0%
        assert_eq!(skew.percent_positive, 0.0);
    }

    #[test]
    fn test_gc_skew_empty() {
        let sequence: &[u8] = b"";
        let (runs, skew) = compute_gc_skew(sequence, 4, 10.0);
        assert!(runs.is_empty());
        assert_eq!(skew.amplitude, 0.0);
        assert_eq!(skew.percent_positive, 0.0);
    }
}
