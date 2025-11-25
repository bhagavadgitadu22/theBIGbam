//! Signal compression functions.
//!
//! This module handles compressing feature signals for efficient storage.
//! Python equivalent: `compress_signal()` in calculating_data.py:91-145

use crate::types::PlotType;

/// Compress signal for storage (subsampling + outlier detection).
/// Python equivalent: `compress_signal()` in calculating_data.py:91-145
/// NOTE: Rust correctly identifies derivative outliers; Python has a bug (see PYTHON_COMPARISON.md)
pub fn compress_signal(
    values: &[f64],
    plot_type: PlotType,
    step: usize,
    z_thresh: f64,
    deriv_thresh: f64,
    max_points: usize,
) -> (Vec<i32>, Vec<f32>) {
    // calculating_data.py:91-140 - compress_signal()
    let n = values.len();  // py:99
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // py:102-103 - Calculate mean and std
    let mean: f64 = (values.iter().sum::<f64>() / n as f64 * 1e9).round() / 1e9;  // py:102 np.mean()
    let variance: f64 = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    let std = (variance.sqrt().max(1e-9) * 1e9).round() / 1e9;  // py:103 np.std()

    // py:104 - Find value outliers: np.abs(feature_values - y_mean) > z_thresh * y_std
    let val_outliers: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, &v)| {
            let z = ((v - mean).abs() * 1e6).round() / 1e6;
            let thresh = (z_thresh * std * 1e6).round() / 1e6;
            z > thresh
        })
        .map(|(i, _)| i)
        .collect();

    let mut keep_idx: Vec<usize> = match plot_type {
        PlotType::Curve => {
            // py:107-108 - Regular subsampling: np.arange(0, n, step)
            let mut regular: Vec<usize> = (0..n).step_by(step).collect();

            // py:110-113 - Derivative outliers
            let mut derivatives = vec![0.0f64; n];
            for i in 1..n {
                derivatives[i] = values[i] - values[i - 1];  // py:111 np.diff()
            }
            let deriv_std = {
                let deriv_mean: f64 = (derivatives.iter().sum::<f64>() / n as f64 * 1e9).round() / 1e9;
                let deriv_var: f64 = derivatives.iter().map(|x| (x - deriv_mean).powi(2)).sum::<f64>() / n as f64;
                (deriv_var.sqrt().max(1e-9) * 1e9).round() / 1e9  // py:112 np.std(dy)
            };

            // py:113 - der_outliers = np.abs(dy) > deriv_thresh * dy_std
            // NOTE: Rust correctly keeps positions i-1 and i for each outlier
            // Python py:116 has a bug with boolean array arithmetic (see PYTHON_COMPARISON.md)
            let deriv_outliers: Vec<usize> = derivatives
                .iter()
                .enumerate()
                .filter(|(_, &d)| {
                    let z = (d.abs() * 1e6).round() / 1e6;
                    let thresh = (deriv_thresh * deriv_std * 1e6).round() / 1e6;
                    z > thresh
                })
                .flat_map(|(i, _)| if i > 0 { vec![i - 1, i] } else { vec![i] })  // correct behavior
                .filter(|&i| i < n)
                .collect();

            // py:123-128 - Combine indices: merge_sorted_unique()
            regular.extend(val_outliers);   // py:124-125
            regular.extend(deriv_outliers);
            if n > 0 {
                regular.push(n - 1);  // py:122-123 - always include last point
            }
            regular.sort_unstable();
            regular.dedup();
            regular
        }
        PlotType::Bars => {
            // py:126-127 - For bars: only keep value outliers
            val_outliers
        }
    };

    // py:137-141 - Apply max_points limit if needed
    if keep_idx.len() > max_points {
        let step_lim = keep_idx.len() / max_points;  // py:138
        keep_idx = keep_idx.into_iter().step_by(step_lim).collect();  // py:139
        if let Some(&last) = keep_idx.last() {
            if last != n - 1 && n > 0 {
                keep_idx.push(n - 1);  // py:140-141 - ensure last point included
            }
        }
    }

    // py:144-145 - Build output: x = np.arange(1, ref_length+1)[keep_idx], y = feature_values[keep_idx]
    // Positions are 1-indexed to match Python's np.arange(1, ref_length+1)
    let xs: Vec<i32> = keep_idx.iter().map(|&i| (i + 1) as i32).collect();  // py:144
    let ys: Vec<f32> = keep_idx.iter().map(|&i| values[i] as f32).collect();  // py:145

    (xs, ys)
}
