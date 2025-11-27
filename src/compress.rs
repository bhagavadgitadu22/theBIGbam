//! Signal compression functions.
//!
//! This module handles compressing feature signals for efficient storage.
//! Python equivalent: `compress_signal()` in calculating_data.py:91-145

use crate::types::{mean_std, PlotType};

/// Compress signal for storage (subsampling + outlier detection).
/// Python equivalent: `compress_signal()` in calculating_data.py:91-145
///
/// If `python_compat` is true, replicates the Python bug where derivative outlier
/// detection uses boolean arithmetic instead of indices, resulting in fewer kept points.
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

    // Calculate value statistics
    let (mean, std) = mean_std(values);

    // Find value outliers: |value - mean| > z_thresh * std
    let val_outliers: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, &v)| (v - mean).abs() > z_thresh * std)
        .map(|(i, _)| i)
        .collect();

    let mut keep_idx: Vec<usize> = match plot_type {
        PlotType::Curve => compute_curve_indices(values, n, step, deriv_thresh, &val_outliers, python_compat),
        PlotType::Bars => val_outliers,
    };

    // Apply max_points limit if needed
    apply_max_points_limit(&mut keep_idx, max_points, n);

    // Build output with 1-indexed positions
    let xs: Vec<i32> = keep_idx.iter().map(|&i| (i + 1) as i32).collect();
    let ys: Vec<f32> = keep_idx.iter().map(|&i| values[i] as f32).collect();

    (xs, ys)
}

/// Compute indices to keep for curve-type features.
fn compute_curve_indices(
    values: &[f64],
    n: usize,
    step: usize,
    deriv_thresh: f64,
    val_outliers: &[usize],
    python_compat: bool,
) -> Vec<usize> {
    // Regular subsampling
    let mut indices: Vec<usize> = (0..n).step_by(step).collect();

    // Compute derivatives and find derivative outliers
    let derivatives: Vec<f64> = std::iter::once(0.0)
        .chain(values.windows(2).map(|w| w[1] - w[0]))
        .collect();

    let (_, deriv_std) = mean_std(&derivatives);

    let deriv_outliers = if python_compat {
        compute_derivative_outliers_python_compat(&derivatives, deriv_thresh, deriv_std, n)
    } else {
        compute_derivative_outliers(&derivatives, deriv_thresh, deriv_std, n)
    };

    // Merge all indices
    indices.extend(val_outliers);
    indices.extend(deriv_outliers);
    if n > 0 {
        indices.push(n - 1); // Always include last point
    }

    indices.sort_unstable();
    indices.dedup();
    indices
}

/// Compute derivative outlier indices (correct implementation).
fn compute_derivative_outliers(derivatives: &[f64], deriv_thresh: f64, deriv_std: f64, n: usize) -> Vec<usize> {
    derivatives
        .iter()
        .enumerate()
        .filter(|(_, &d)| d.abs() > deriv_thresh * deriv_std)
        .flat_map(|(i, _)| {
            if i > 0 {
                vec![i - 1, i]
            } else {
                vec![i]
            }
        })
        .filter(|&i| i < n)
        .collect()
}

/// Compute derivative outlier indices (Python-compatible buggy version).
///
/// Replicates Python bug: np.concatenate([der_outliers - 1, der_outliers])
/// where der_outliers is boolean, so True-1=0, False-1=-1 (clipped to 0).
fn compute_derivative_outliers_python_compat(
    derivatives: &[f64],
    deriv_thresh: f64,
    deriv_std: f64,
    n: usize,
) -> Vec<usize> {
    let bools: Vec<bool> = derivatives
        .iter()
        .map(|&d| d.abs() > deriv_thresh * deriv_std)
        .collect();

    // Simulate: np.concatenate([bools - 1, bools]) then clip to [0, n-1]
    let mut indices: Vec<i64> = Vec::with_capacity(bools.len() * 2);

    for &b in &bools {
        indices.push(if b { 0 } else { -1 }); // True-1=0, False-1=-1
    }
    for &b in &bools {
        indices.push(if b { 1 } else { 0 }); // True=1, False=0
    }

    // Clip to [0, n-1] and deduplicate
    let mut result: Vec<usize> = indices
        .into_iter()
        .map(|i| i.max(0).min((n - 1) as i64) as usize)
        .collect();

    result.sort_unstable();
    result.dedup();
    result
}

/// Apply max_points limit if needed.
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
