//! Shared per-contig result construction for contig and MAG processing.
use crate::features::FeatureArrays;
use crate::types::PresenceData;

/// The caller supplies the mean captured at its existing processing stage.
/// MAG summaries precede feature encoding; contig summaries follow it.
pub(super) fn build_presence(
    arrays: &FeatureArrays,
    contig_id: i64,
    contig_name: &str,
    coverage_pct: f64,
    primary_count: u64,
    coverage_mean: f32,
) -> PresenceData {
    let (median, trimmed_mean) = arrays.coverage_median_and_trimmed_mean(0.05);
    let coverage_median = median as f32;
    let coverage_trimmed_mean = trimmed_mean as f32;
    let expected_af = (1.0 - (-0.883 * coverage_mean as f64).exp()) * 100.0;
    let coverage_relative_coverage_roughness = if arrays.primary_reads.len() > 1
        && coverage_mean > 0.0
    {
        let nn = arrays.primary_reads.len();
        let mean_cov = coverage_mean as f64;
        let sq_diff: f64 = arrays
            .primary_reads
            .windows(2)
            .map(|w| {
                let d = w[1] as f64 - w[0] as f64;
                d * d
            })
            .sum::<f64>()
            / (nn - 1) as f64;
        (sq_diff.sqrt() / mean_cov
            * crate::types::get_column_scale("Coverage", "Coverage_relative_coverage_roughness"))
            as f32
    } else {
        0.0
    };
    let coverage_coefficient_of_variation =
        if !arrays.primary_reads.is_empty() && coverage_mean > 0.0 {
            let mean_cov = coverage_mean as f64;
            let nn = arrays.primary_reads.len() as f64;
            let var: f64 = arrays
                .primary_reads
                .iter()
                .map(|&x| {
                    let d = x as f64 - mean_cov;
                    d * d
                })
                .sum::<f64>()
                / nn;
            (var.sqrt() / mean_cov
                * crate::types::get_column_scale("Coverage", "Coverage_coefficient_of_variation"))
                as f32
        } else {
            0.0
        };

    PresenceData {
        contig_id: contig_id,
        contig_name: contig_name.to_owned(),
        coverage_pct: coverage_pct as f32,
        expected_aligned_fraction: expected_af as f32,
        read_count: primary_count,
        coverage_relative_coverage_roughness,
        coverage_coefficient_of_variation,
        coverage_mean,
        coverage_median,
        coverage_trimmed_mean,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_presence_preserves_summary_and_empty_coverage() {
        let mut arrays = FeatureArrays::new(4);
        arrays.primary_reads.copy_from_slice(&[0, 2, 2, 4]);
        let presence = build_presence(&arrays, 7, "a", 75.0, 3, 2.0);
        assert_eq!((presence.contig_id, presence.read_count), (7, 3));
        assert_eq!(
            (
                presence.coverage_pct,
                presence.coverage_mean,
                presence.coverage_median,
                presence.coverage_trimmed_mean
            ),
            (75.0, 2.0, 2.0, 2.0)
        );
        let expected_roughness = (8.0_f64 / 3.0).sqrt() / 2.0
            * crate::types::get_column_scale("Coverage", "Coverage_relative_coverage_roughness");
        let expected_cv = 2.0_f64.sqrt() / 2.0
            * crate::types::get_column_scale("Coverage", "Coverage_coefficient_of_variation");
        assert_eq!(
            presence.coverage_relative_coverage_roughness,
            expected_roughness as f32
        );
        assert_eq!(
            presence.coverage_coefficient_of_variation,
            expected_cv as f32
        );
        arrays.primary_reads.fill(0);
        let empty = build_presence(&arrays, 7, "a", 0.0, 0, 0.0);
        assert_eq!(
            (
                empty.coverage_mean,
                empty.coverage_relative_coverage_roughness,
                empty.coverage_coefficient_of_variation,
                empty.expected_aligned_fraction
            ),
            (0.0, 0.0, 0.0, 0.0)
        );
    }
}
