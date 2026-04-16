//! GC content computation for genomic sequences.
//!
//! Computes GC content (percentage of G and C bases) using a sliding window approach
//! and applies RLE compression for efficient storage.

/// Default window size for GC content computation (500bp).
pub const DEFAULT_GC_CONTENT_WINDOW_SIZE: usize = 500;
/// Default window size for GC skew computation (1000bp).
pub const DEFAULT_GC_SKEW_WINDOW_SIZE: usize = 1000;

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
            gc_content_window_size: DEFAULT_GC_CONTENT_WINDOW_SIZE,
            gc_skew_window_size: DEFAULT_GC_SKEW_WINDOW_SIZE,
        }
    }
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

/// Compute GC content using non-overlapping windows.
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of each non-overlapping window in base pairs (typically 500)
///
/// # Algorithm
/// 1. Divide the sequence into non-overlapping windows of window_size
/// 2. For each window, compute GC% (excluding N bases)
///
/// # Returns
/// Tuple of (values, stats) where:
/// - values: Vector of GC percentages (0-100), one per window
/// - stats: GCStats with average and sd GC percentages
pub fn compute_gc_content(sequence: &[u8], window_size: usize) -> (Vec<u8>, GCStats) {
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

    (gc_values, stats)
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

/// Compute GC skew using non-overlapping windows.
///
/// GC skew = (G - C) / (G + C) at each position
/// - Range: -1 to +1
/// - Positive = more G than C (leading strand bias)
/// - Negative = more C than G (lagging strand bias)
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of each non-overlapping window in base pairs (typically 1000)
///
/// # Returns
/// Tuple of (values, stats) where:
/// - values: Vector of GC skew × 100 (range: -100 to +100), one per window
/// - stats: GCSkewStats with amplitude and percent_positive
pub fn compute_gc_skew(sequence: &[u8], window_size: usize) -> (Vec<i16>, GCSkewStats) {
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

    (skew_values, stats)
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_simple() {
        // ATGC sequence - 50% GC per window
        let sequence = b"ATGCATGCATGC";
        let (values, stats) = compute_gc_content(sequence, 4);

        // Should produce 3 values (window size 4, 12bp sequence)
        assert_eq!(values.len(), 3);
        // Each window (ATGC, ATGC, ATGC) should be 50% GC
        for value in &values {
            assert_eq!(*value, 50);
        }
        // Stats should also be around 50%
        assert_eq!(stats.average, 50.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_high_gc() {
        // All GC
        let sequence = b"GGGGCCCCGGGGCCCC";
        let (values, stats) = compute_gc_content(sequence, 4);

        // Should produce 4 values (window size 4, 16bp sequence)
        assert_eq!(values.len(), 4);
        // Each window should be 100% GC
        for value in &values {
            assert_eq!(*value, 100);
        }
        // Stats should be 100% with 0 sd
        assert_eq!(stats.average, 100.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_low_gc() {
        // All AT
        let sequence = b"AAAATTTTAAAATTTT";
        let (values, stats) = compute_gc_content(sequence, 4);

        // Should produce 4 values (window size 4, 16bp sequence)
        assert_eq!(values.len(), 4);
        // Each window should be 0% GC
        for value in &values {
            assert_eq!(*value, 0);
        }
        // Stats should be 0% with 0 sd
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_with_n() {
        // Sequence with N bases - N should be excluded
        let sequence = b"ATGCNNNNATGC";
        let (values, _stats) = compute_gc_content(sequence, 4);

        // Should produce 3 values (window size 4, 12bp sequence)
        assert_eq!(values.len(), 3);
    }

    #[test]
    fn test_gc_content_empty() {
        let sequence: &[u8] = b"";
        let (values, stats) = compute_gc_content(sequence, 4);
        assert!(values.is_empty());
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_individual_windows() {
        // Verify that each window value is computed correctly
        let sequence = b"ATGCATGCATGCATGCATGC"; // 20 bases, 5 windows of 4bp
        let (values, _stats) = compute_gc_content(sequence, 4);

        // Should have 5 windows
        assert_eq!(values.len(), 5);
        // Each window (ATGC) should be 50% GC
        for value in &values {
            assert_eq!(*value, 50);
        }
    }

    #[test]
    fn test_varying_gc() {
        // Long sequence with uniform GC should give consistent values
        let sequence: Vec<u8> = b"ATGC".repeat(100);
        let (values, stats) = compute_gc_content(&sequence, 100);

        // Should produce 4 windows (400bp sequence, 100bp windows)
        assert_eq!(values.len(), 4);
        // All values should be 50% GC
        for value in &values {
            assert_eq!(*value, 50);
        }
        // Stats should be around 50%
        assert_eq!(stats.average, 50.0);
        assert_eq!(stats.sd, 0.0);
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
        let (_runs, skew) = compute_gc_skew(sequence, 4);

        // Amplitude should be 0 (or very close)
        assert!(skew.amplitude < 0.1);
        // Percent positive should be around 50% (or 0% if exactly 0)
        // With equal G and C, skew is 0, so no positions are strictly positive
    }

    #[test]
    fn test_gc_skew_g_rich() {
        // All G - skew should be +1 everywhere
        let sequence = b"GGGGGGGGGGGGGGGG";
        let (_runs, skew) = compute_gc_skew(sequence, 4);

        // Amplitude should be 0 (constant +1)
        assert_eq!(skew.amplitude, 0.0);
        // Percent positive should be 100%
        assert_eq!(skew.percent_positive, 100.0);
    }

    #[test]
    fn test_gc_skew_c_rich() {
        // All C - skew should be -1 everywhere
        let sequence = b"CCCCCCCCCCCCCCCC";
        let (_runs, skew) = compute_gc_skew(sequence, 4);

        // Amplitude should be 0 (constant -1)
        assert_eq!(skew.amplitude, 0.0);
        // Percent positive should be 0%
        assert_eq!(skew.percent_positive, 0.0);
    }

    #[test]
    fn test_gc_skew_empty() {
        let sequence: &[u8] = b"";
        let (runs, skew) = compute_gc_skew(sequence, 4);
        assert!(runs.is_empty());
        assert_eq!(skew.amplitude, 0.0);
        assert_eq!(skew.percent_positive, 0.0);
    }
}
