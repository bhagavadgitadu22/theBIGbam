//! GC content computation for genomic sequences.
//!
//! Computes GC content (percentage of G and C bases) using a sliding window approach
//! and applies RLE compression for efficient storage.

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

/// Compute GC content using a sliding window approach with RLE compression.
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of the sliding window in base pairs (typically 100)
/// - `contig_variation_percentage`: RLE compression tolerance for contig-level features (default 0.1%)
///
/// # Algorithm
/// 1. For each position, compute GC% in a centered window
/// 2. N bases are excluded from both numerator and denominator
/// 3. Apply RLE compression: merge consecutive positions with similar GC%
///
/// # Returns
/// Vector of GC content runs with (start_pos, end_pos, gc_percentage)
pub fn compute_gc_content(sequence: &[u8], window_size: usize, contig_variation_percentage: f64) -> Vec<GCContentRun> {
    let n = sequence.len();
    if n == 0 {
        return Vec::new();
    }

    let half_window = window_size / 2;
    let mut gc_values: Vec<u8> = Vec::with_capacity(n);

    // Compute GC content for each position using a centered sliding window
    for i in 0..n {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window + 1).min(n);

        let (gc_count, valid_count) = count_gc_in_window(&sequence[start..end]);

        let gc_pct = if valid_count > 0 {
            ((gc_count as f64 / valid_count as f64) * 100.0).round() as u8
        } else {
            50 // Default to 50% if window contains only Ns
        };

        gc_values.push(gc_pct);
    }

    // Apply RLE compression with user-configurable tolerance
    compress_gc_values(&gc_values, contig_variation_percentage)
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

/// Compress GC values using adaptive RLE.
/// Uses the formula: |x[i] - x[i-1]| <= ratio × min(x[i], x[i-1])
/// Input is percentage (e.g., 10 for 10%), converted to fraction internally.
fn compress_gc_values(values: &[u8], percentage: f64) -> Vec<GCContentRun> {
    if values.is_empty() {
        return Vec::new();
    }

    let ratio = percentage * 0.01; // Convert percentage to fraction
    let mut runs = Vec::new();

    let mut run_start = 0;
    let mut run_value = values[0] as f64;
    let mut run_sum = values[0] as f64;
    let mut run_count = 1;

    for i in 1..values.len() {
        let val = values[i] as f64;

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
            runs.push(GCContentRun {
                start_pos: (run_start + 1) as i32,
                end_pos: (run_start + run_count) as i32,
                gc_percentage: run_value.round() as u8,
            });

            // Start new run
            run_start = i;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Save the last run
    runs.push(GCContentRun {
        start_pos: (run_start + 1) as i32,
        end_pos: (run_start + run_count) as i32,
        gc_percentage: run_value.round() as u8,
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
        let runs = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with ~50% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert!(run.gc_percentage >= 40 && run.gc_percentage <= 60);
        }
    }

    #[test]
    fn test_gc_content_high_gc() {
        // All GC
        let sequence = b"GGGGCCCCGGGGCCCC";
        let runs = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 100% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 100);
        }
    }

    #[test]
    fn test_gc_content_low_gc() {
        // All AT
        let sequence = b"AAAATTTTAAAATTTT";
        let runs = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 0% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 0);
        }
    }

    #[test]
    fn test_gc_content_with_n() {
        // Sequence with N bases - N should be excluded
        let sequence = b"ATGCNNNNATGC";
        let runs = compute_gc_content(sequence, 4, 10.0);

        // Should still compute ~50% GC from valid bases
        assert!(!runs.is_empty());
    }

    #[test]
    fn test_gc_content_empty() {
        let sequence: &[u8] = b"";
        let runs = compute_gc_content(sequence, 4, 10.0);
        assert!(runs.is_empty());
    }

    #[test]
    fn test_full_coverage() {
        // Verify that the runs cover the entire sequence
        let sequence = b"ATGCATGCATGCATGCATGC"; // 20 bases
        let runs = compute_gc_content(sequence, 4, 10.0);

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
        let runs = compute_gc_content(&sequence, 100, 10.0);

        // Should compress but still cover full length
        assert!(!runs.is_empty());
        assert_eq!(runs.first().unwrap().start_pos, 1);
        assert_eq!(runs.last().unwrap().end_pos, 4000);
    }
}
