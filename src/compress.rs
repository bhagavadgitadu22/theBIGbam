//! Signal compression functions.
//!
//! # Bioinformatics Background
//!
//! Raw genomic signals (coverage, tau, etc.) can have millions of data points.
//! Storing and visualizing all points would be inefficient. This module compresses
//! signals by keeping:
//! 1. Regular samples every N positions (subsampling)
//! 2. Outlier points with extreme values (z-score threshold)
//! 3. Points with large derivatives (rapid changes)
//!
//! This preserves the signal's shape while reducing storage dramatically.
//!
//! # Python Equivalent
//!
//! This module corresponds to `compress_signal()` in `calculating_data.py:91-145`:
//! ```python
//! def compress_signal(type_picked, feature_values, ref_length, step, z_thresh, deriv_thresh, max_points):
//!     # Value outliers using z-score
//!     y_mean = np.mean(feature_values)
//!     y_std = np.std(feature_values) or 1e-9
//!     val_outliers = np.abs(feature_values - y_mean) > z_thresh * y_std
//!
//!     # Regular subsampling + derivative outliers
//!     if type_picked == "curve":
//!         regular_idx = np.arange(0, n, step)
//!         dy = np.diff(feature_values, prepend=feature_values[0])
//!         der_outliers = np.abs(dy) > deriv_thresh * dy_std
//!         # ... merge and keep points
//! ```
//!
//! # Rust Concepts
//!
//! ## Slices vs Vecs
//! The function takes `&[f64]` (a slice reference) as input - this borrows data
//! without copying. It returns `Vec<i32>` and `Vec<f32>` (owned vectors) because
//! the caller needs to own the new compressed data.
//!
//! ## Iterator Chains
//! Rust iterators chain operations like `.iter().enumerate().filter().map().collect()`.
//! This is similar to Python's list comprehensions but evaluated lazily. The compiler
//! optimizes these chains into efficient loops.

use crate::types::{mean_std, PlotType};

// ============================================================================
// MAIN COMPRESSION FUNCTION
// ============================================================================

/// Compress signal for storage (subsampling + outlier detection).
///
/// # Algorithm
/// 1. Find value outliers: points where |value - mean| > z_thresh * std
/// 2. For curve type: also keep regularly spaced samples and derivative outliers
/// 3. For bars type: only keep value outliers
/// 4. Enforce max_points limit by further subsampling if needed
///
/// # Parameters
/// - `values`: The signal to compress (e.g., coverage at each position)
/// - `plot_type`: Curve (line plot) or Bars (bar chart)
/// - `step`: Keep every Nth point for regular subsampling
/// - `z_thresh`: Z-score threshold for value outliers (e.g., 3.0)
/// - `deriv_thresh`: Z-score threshold for derivative outliers
/// - `max_points`: Hard limit on output points
/// - `python_compat`: If true, replicate Python's buggy behavior (for testing)
///
/// # Returns
/// Tuple of (x_positions, y_values) where x is 1-indexed (for database storage)
///
/// # Python Equivalent
/// This is a direct port of `compress_signal()` from `calculating_data.py:91-145`.
pub fn compress_signal(
    values: &[f64],
    plot_type: PlotType,
    step: usize,
    z_thresh: f64,
    deriv_thresh: f64,
    max_points: usize,
    python_compat: bool,
) -> (Vec<i32>, Vec<f32>) {
    let n = values.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // -----------------------------------------------------------------------
    // Step 1: Calculate value statistics for outlier detection
    // Python: y_mean = np.mean(feature_values); y_std = np.std(feature_values)
    // -----------------------------------------------------------------------
    let (mean, std) = mean_std(values);

    // -----------------------------------------------------------------------
    // Step 2: Find value outliers: |value - mean| > z_thresh * std
    // Python: val_outliers = np.abs(feature_values - y_mean) > z_thresh * y_std
    //
    // Rust Concept: Iterator chain
    // .iter() - create iterator over references
    // .enumerate() - add index to each item: (index, &value)
    // .filter() - keep only items matching predicate
    // .map() - transform items (here: extract just the index)
    // .collect() - gather into a Vec
    // -----------------------------------------------------------------------
    let val_outliers: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, &v)| (v - mean).abs() > z_thresh * std)
        .map(|(i, _)| i)
        .collect();

    // -----------------------------------------------------------------------
    // Step 3: Compute indices to keep based on plot type
    // Python: if type_picked == "curve": ... elif type_picked == "bars": ...
    // -----------------------------------------------------------------------
    let mut keep_idx: Vec<usize> = match plot_type {
        PlotType::Curve => compute_curve_indices(values, n, step, deriv_thresh, &val_outliers, python_compat),
        PlotType::Bars => val_outliers,  // Bars only keep outliers
    };

    // -----------------------------------------------------------------------
    // Step 4: Apply max_points limit if needed
    // Python: if len(keep_idx) > max_points: ...
    // -----------------------------------------------------------------------
    apply_max_points_limit(&mut keep_idx, max_points, n);

    // -----------------------------------------------------------------------
    // Step 5: Build output with 1-indexed positions (for database compatibility)
    // Python: x = np.arange(1, ref_length + 1); return {"x": x[keep_idx], ...}
    // -----------------------------------------------------------------------
    let xs: Vec<i32> = keep_idx.iter().map(|&i| (i + 1) as i32).collect();
    let ys: Vec<f32> = keep_idx.iter().map(|&i| values[i] as f32).collect();

    (xs, ys)
}

// ============================================================================
// CURVE-TYPE COMPRESSION
// ============================================================================

/// Compute indices to keep for curve-type features.
///
/// # Python Equivalent
/// ```python
/// # From compress_signal() in calculating_data.py:111-130
/// regular_idx = np.arange(0, n, step, dtype=int)
/// dy = np.diff(feature_values, prepend=feature_values[0])
/// dy_std = np.std(dy) or 1e-9
/// der_outliers = np.abs(dy) > deriv_thresh * dy_std
/// der_outliers = np.unique(np.clip(np.concatenate([der_outliers - 1, der_outliers]), 0, n - 1))
/// keep_idx = merge_sorted_unique(regular_idx, outlier_idx, der_outliers, last_idx)
/// ```
fn compute_curve_indices(
    values: &[f64],
    n: usize,
    step: usize,
    deriv_thresh: f64,
    val_outliers: &[usize],
    python_compat: bool,
) -> Vec<usize> {
    // Regular subsampling: keep every `step`th point
    // Python: regular_idx = np.arange(0, n, step, dtype=int)
    let mut indices: Vec<usize> = (0..n).step_by(step).collect();

    // Compute derivatives (first differences) for derivative outlier detection
    // Python: dy = np.diff(feature_values, prepend=feature_values[0])
    let derivatives: Vec<f64> = std::iter::once(0.0)
        .chain(values.windows(2).map(|w| w[1] - w[0]))
        .collect();

    let (_, deriv_std) = mean_std(&derivatives);

    // Find derivative outliers using the appropriate method
    let deriv_outliers = if python_compat {
        compute_derivative_outliers_python_compat(&derivatives, deriv_thresh, deriv_std, n)
    } else {
        compute_derivative_outliers(&derivatives, deriv_thresh, deriv_std, n)
    };

    // Merge all indices: regular samples + value outliers + derivative outliers
    // Python: keep_idx = merge_sorted_unique(regular_idx, outlier_idx, ...)
    indices.extend(val_outliers);
    indices.extend(deriv_outliers);
    if n > 0 {
        indices.push(n - 1); // Always include last point
    }

    // Sort and deduplicate
    indices.sort_unstable();
    indices.dedup();
    indices
}

// ============================================================================
// DERIVATIVE OUTLIER DETECTION
// ============================================================================

/// Compute derivative outlier indices (correct implementation).
///
/// For each position where |derivative| > threshold, keeps both that position
/// and the position before it (to capture the change point).
fn compute_derivative_outliers(derivatives: &[f64], deriv_thresh: f64, deriv_std: f64, n: usize) -> Vec<usize> {
    derivatives
        .iter()
        .enumerate()
        .filter(|(_, &d)| d.abs() > deriv_thresh * deriv_std)
        .flat_map(|(i, _)| {
            // Include both the point before and at the derivative spike
            if i > 0 {
                vec![i - 1, i]
            } else {
                vec![i]
            }
        })
        .filter(|&i| i < n)
        .collect()
}

/// Compute derivative outlier indices (Python-compatible version).
///
/// Replicates a bug in the original Python implementation where
/// `np.concatenate([der_outliers - 1, der_outliers])` operates on a boolean
/// array instead of indices. This results in: True-1=0, False-1=-1 (clipped to 0).
///
/// # Python Equivalent
/// ```python
/// der_outliers = np.abs(dy) > deriv_thresh * dy_std  # Boolean array!
/// der_outliers = np.unique(np.clip(np.concatenate([der_outliers - 1, der_outliers]), 0, n - 1))
/// ```
fn compute_derivative_outliers_python_compat(
    derivatives: &[f64],
    deriv_thresh: f64,
    deriv_std: f64,
    n: usize,
) -> Vec<usize> {
    // Create boolean array matching Python behavior
    let bools: Vec<bool> = derivatives
        .iter()
        .map(|&d| d.abs() > deriv_thresh * deriv_std)
        .collect();

    // Simulate: np.concatenate([bools - 1, bools]) where bools is boolean
    // True-1=0 (True treated as 1), False-1=-1
    let mut indices: Vec<i64> = Vec::with_capacity(bools.len() * 2);

    for &b in &bools {
        indices.push(if b { 0 } else { -1 }); // True-1=0, False-1=-1
    }
    for &b in &bools {
        indices.push(if b { 1 } else { 0 }); // True=1, False=0
    }

    // Clip to [0, n-1] and deduplicate (np.unique)
    let mut result: Vec<usize> = indices
        .into_iter()
        .map(|i| i.max(0).min((n - 1) as i64) as usize)
        .collect();

    result.sort_unstable();
    result.dedup();
    result
}

// ============================================================================
// POINT LIMIT ENFORCEMENT
// ============================================================================

/// Apply max_points limit if needed.
///
/// If the number of kept indices exceeds max_points, subsample further
/// while ensuring the last point is always included.
///
/// # Python Equivalent
/// ```python
/// if len(keep_idx) > max_points:
///     step_lim = len(keep_idx) // max_points
///     keep_idx = keep_idx[::step_lim]
///     if keep_idx[-1] != n - 1:
///         keep_idx = np.append(keep_idx, n - 1)
/// ```
fn apply_max_points_limit(keep_idx: &mut Vec<usize>, max_points: usize, n: usize) {
    if keep_idx.len() > max_points {
        let step_lim = keep_idx.len() / max_points;
        *keep_idx = keep_idx.iter().copied().step_by(step_lim).collect();

        // Ensure last point is included
        if let Some(&last) = keep_idx.last() {
            if last != n - 1 && n > 0 {
                keep_idx.push(n - 1);
            }
        }
    }
}
