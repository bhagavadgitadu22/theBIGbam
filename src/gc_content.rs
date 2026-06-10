//! GC content computation for genomic sequences.
//!
//! Computes GC content (percentage of G and C bases) and GC skew using
//! non-overlapping windows in a single pass over the sequence.

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

/// Compute GC content and GC skew in a single pass over the sequence.
///
/// Uses two independent sets of non-overlapping windows (typically 500bp for
/// GC content and 1000bp for GC skew). Each byte is classified once and both
/// window accumulators are updated simultaneously.
pub fn compute_gc_content_and_skew(
    sequence: &[u8],
    gc_window_size: usize,
    skew_window_size: usize,
) -> (Vec<u8>, GCStats, Vec<i16>, GCSkewStats) {
    let n = sequence.len();
    let empty_gc = GCStats { average: 0.0, sd: 0.0 };
    let empty_skew = GCSkewStats { amplitude: 0.0, percent_positive: 0.0 };
    if n == 0 {
        return (Vec::new(), empty_gc, Vec::new(), empty_skew);
    }

    let mut gc_values: Vec<u8> = Vec::with_capacity(n / gc_window_size + 1);
    let mut skew_values: Vec<i16> = Vec::with_capacity(n / skew_window_size + 1);
    let skew_scale = crate::types::get_value_scale("gc_skew").multiplier();

    // GC content window accumulators
    let mut gc_g: usize = 0;
    let mut gc_c: usize = 0;
    let mut gc_at: usize = 0;
    let mut gc_pos: usize = 0;

    // GC skew window accumulators
    let mut skew_g: usize = 0;
    let mut skew_c: usize = 0;
    let mut skew_pos: usize = 0;

    for &base in sequence {
        match base {
            b'G' | b'g' => { gc_g += 1; skew_g += 1; }
            b'C' | b'c' => { gc_c += 1; skew_c += 1; }
            b'A' | b'a' | b'T' | b't' => { gc_at += 1; }
            _ => {}
        }
        gc_pos += 1;
        skew_pos += 1;

        if gc_pos == gc_window_size {
            let gc_count = gc_g + gc_c;
            let valid = gc_count + gc_at;
            let pct = if valid > 0 {
                ((gc_count as f64 / valid as f64) * 100.0).round() as u8
            } else {
                50
            };
            gc_values.push(pct);
            gc_g = 0; gc_c = 0; gc_at = 0; gc_pos = 0;
        }

        if skew_pos == skew_window_size {
            let gc_sum = skew_g + skew_c;
            let skew = if gc_sum > 0 {
                (skew_g as f64 - skew_c as f64) / gc_sum as f64
            } else {
                0.0
            };
            skew_values.push((skew * skew_scale).round() as i16);
            skew_g = 0; skew_c = 0; skew_pos = 0;
        }
    }

    // Trailing partial GC content window
    if gc_pos > 0 {
        let gc_count = gc_g + gc_c;
        let valid = gc_count + gc_at;
        let pct = if valid > 0 {
            ((gc_count as f64 / valid as f64) * 100.0).round() as u8
        } else {
            50
        };
        gc_values.push(pct);
    }

    // Trailing partial GC skew window
    if skew_pos > 0 {
        let gc_sum = skew_g + skew_c;
        let skew = if gc_sum > 0 {
            (skew_g as f64 - skew_c as f64) / gc_sum as f64
        } else {
            0.0
        };
        skew_values.push((skew * skew_scale).round() as i16);
    }

    // GC content stats
    let gc_stats = compute_gc_stats(&gc_values);

    // GC skew stats
    let skew_stats = if skew_values.is_empty() {
        empty_skew
    } else {
        let min_skew = skew_values.iter().cloned().min().unwrap_or(0) as f64 / skew_scale;
        let max_skew = skew_values.iter().cloned().max().unwrap_or(0) as f64 / skew_scale;
        let amplitude = (max_skew - min_skew) as f32;
        let positive_count = skew_values.iter().filter(|&&s| s > 0).count();
        let percent_positive = (positive_count as f64 / skew_values.len() as f64 * 100.0) as f32;
        GCSkewStats { amplitude, percent_positive }
    };

    (gc_values, gc_stats, skew_values, skew_stats)
}

fn compute_gc_stats(gc_values: &[u8]) -> GCStats {
    if gc_values.is_empty() {
        return GCStats { average: 0.0, sd: 0.0 };
    }

    let n = gc_values.len() as f64;
    let sum: f64 = gc_values.iter().map(|&v| v as f64).sum();
    let average = sum / n;

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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_simple() {
        let sequence = b"ATGCATGCATGC";
        let (values, stats, _, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(values.len(), 3);
        for value in &values {
            assert_eq!(*value, 50);
        }
        assert_eq!(stats.average, 50.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_high_gc() {
        let sequence = b"GGGGCCCCGGGGCCCC";
        let (values, stats, _, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(values.len(), 4);
        for value in &values {
            assert_eq!(*value, 100);
        }
        assert_eq!(stats.average, 100.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_low_gc() {
        let sequence = b"AAAATTTTAAAATTTT";
        let (values, stats, _, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(values.len(), 4);
        for value in &values {
            assert_eq!(*value, 0);
        }
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_content_with_n() {
        let sequence = b"ATGCNNNNATGC";
        let (values, _, _, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(values.len(), 3);
    }

    #[test]
    fn test_gc_content_empty() {
        let sequence: &[u8] = b"";
        let (values, stats, skew_values, skew_stats) = compute_gc_content_and_skew(sequence, 4, 4);
        assert!(values.is_empty());
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
        assert!(skew_values.is_empty());
        assert_eq!(skew_stats.amplitude, 0.0);
        assert_eq!(skew_stats.percent_positive, 0.0);
    }

    #[test]
    fn test_individual_windows() {
        let sequence = b"ATGCATGCATGCATGCATGC"; // 20 bases, 5 windows of 4bp
        let (values, _, _, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(values.len(), 5);
        for value in &values {
            assert_eq!(*value, 50);
        }
    }

    #[test]
    fn test_varying_gc() {
        let sequence: Vec<u8> = b"ATGC".repeat(100);
        let (values, stats, _, _) = compute_gc_content_and_skew(&sequence, 100, 100);

        assert_eq!(values.len(), 4);
        for value in &values {
            assert_eq!(*value, 50);
        }
        assert_eq!(stats.average, 50.0);
        assert_eq!(stats.sd, 0.0);
    }

    #[test]
    fn test_gc_stats() {
        let gc_values: Vec<u8> = vec![40, 50, 60];
        let stats = compute_gc_stats(&gc_values);

        assert_eq!(stats.average, 50.0);
        assert!(stats.sd > 8.0 && stats.sd < 9.0);
    }

    #[test]
    fn test_gc_skew_balanced() {
        let sequence = b"GCGCGCGCGCGCGCGC";
        let (_, _, _, skew) = compute_gc_content_and_skew(sequence, 4, 4);

        assert!(skew.amplitude < 0.1);
    }

    #[test]
    fn test_gc_skew_g_rich() {
        let sequence = b"GGGGGGGGGGGGGGGG";
        let (_, _, _, skew) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(skew.amplitude, 0.0);
        assert_eq!(skew.percent_positive, 100.0);
    }

    #[test]
    fn test_gc_skew_c_rich() {
        let sequence = b"CCCCCCCCCCCCCCCC";
        let (_, _, _, skew) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(skew.amplitude, 0.0);
        assert_eq!(skew.percent_positive, 0.0);
    }

    #[test]
    fn test_gc_skew_empty() {
        let sequence: &[u8] = b"";
        let (_, _, skew_values, skew) = compute_gc_content_and_skew(sequence, 4, 4);
        assert!(skew_values.is_empty());
        assert_eq!(skew.amplitude, 0.0);
        assert_eq!(skew.percent_positive, 0.0);
    }

    #[test]
    fn test_different_window_sizes() {
        // 16 bases: GC window=4 → 4 values, skew window=8 → 2 values
        let sequence = b"GGGGCCCCGGGGCCCC";
        let (gc_values, _, skew_values, _) = compute_gc_content_and_skew(sequence, 4, 8);

        assert_eq!(gc_values.len(), 4);
        assert_eq!(skew_values.len(), 2);
    }

    #[test]
    fn test_trailing_partial_windows() {
        // 10 bases with window=4 → 2 full windows + 1 partial (2 bases)
        let sequence = b"ATGCATGCAT";
        let (gc_values, _, skew_values, _) = compute_gc_content_and_skew(sequence, 4, 4);

        assert_eq!(gc_values.len(), 3); // 4+4+2
        assert_eq!(skew_values.len(), 3);
    }
}
