//! MAG-level helpers.
//!
//! - `build_mag_member_lookup` — builds contig→MAG mapping for the processing
//!   side to group jobs by MAG.
//! - `compute_mag_coverage_stats_from_raw` — MAG coverage stats from raw
//!   accumulated per-position values (no decode needed).

use std::collections::HashMap;

use anyhow::Result;

use crate::db::DbWriter;

/// Build contig_name → (mag_name, mag_id, offset_in_mag, mag_length) lookup.
/// Used by the processing side to group jobs by MAG before par_iter.
pub fn build_mag_member_lookup(db_writer: &DbWriter) -> Result<crate::processing::MagMemberLookup> {
    if !db_writer.is_mag_mode() {
        return Ok(HashMap::new());
    }
    let conn = db_writer.lock_conn()?;
    let mut stmt = conn.prepare(
        "SELECT c.Contig_name, m.MAG_name, mca.MAG_id, mca.Offset_in_MAG,
                mag_len.MAG_length
         FROM MAG_contigs_association mca
         JOIN Contig c ON c.Contig_id = mca.Contig_id
         JOIN MAG m ON m.MAG_id = mca.MAG_id
         JOIN (
             SELECT mca2.MAG_id, SUM(c2.Contig_length) AS MAG_length
             FROM MAG_contigs_association mca2
             JOIN Contig c2 ON c2.Contig_id = mca2.Contig_id
             GROUP BY mca2.MAG_id
         ) mag_len ON mag_len.MAG_id = mca.MAG_id"
    )?;
    let rows = stmt.query_map([], |r| Ok((
        r.get::<_, String>(0)?,
        r.get::<_, String>(1)?,
        r.get::<_, i64>(2)?,
        r.get::<_, i64>(3)?,
        r.get::<_, i64>(4)?,
    )))?;
    let mut out = HashMap::new();
    for row in rows {
        let (contig_name, mag_name, mag_id, offset, mag_length) = row?;
        out.insert(contig_name, (mag_name, mag_id, offset as u32, mag_length as u32));
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// MAG-level coverage metrics (recomputed from raw per-position values)
// ---------------------------------------------------------------------------

/// Pre-scaled integer fields ready for MAG_coverage INSERT.
/// Scaling matches the per-contig Coverage table so the view can reuse
/// the same /10 and /1_000_000 divisors.
pub struct MagCoverageStats {
    pub mag_length: u32,
    pub aligned_fraction_pct_x10: i32,
    pub expected_aligned_fraction_pct_x10: i32,
    pub read_count: u64,
    pub coverage_mean_x10: i32,
    pub coverage_median_x10: i32,
    pub coverage_trimmed_mean_x10: i32,
    pub coverage_cv_x1e6: i32,
    pub coverage_roughness_x1e6: i32,
}

/// Compute MAG-level coverage stats from pre-accumulated raw per-position
/// i32 values, split by contig boundaries. `contig_lengths` gives the length
/// of each member contig in MAG order so cross-contig jumps are excluded from
/// the roughness metric.
pub fn compute_mag_coverage_stats_from_raw(
    mag_values: &[i32],
    mag_length: u32,
    contig_lengths: &[u32],
    total_read_count: u64,
) -> MagCoverageStats {
    if mag_values.is_empty() {
        return MagCoverageStats {
            mag_length, aligned_fraction_pct_x10: 0, expected_aligned_fraction_pct_x10: 0,
            read_count: total_read_count, coverage_mean_x10: 0, coverage_median_x10: 0,
            coverage_trimmed_mean_x10: 0, coverage_cv_x1e6: 0, coverage_roughness_x1e6: 0,
        };
    }
    let mut sum_x: u128 = 0;
    let mut sum_x2: u128 = 0;
    let mut sum_sq_diff: u128 = 0;
    let mut n_covered: u64 = 0;
    let n_positions = mag_values.len() as u64;
    let mut n_contigs_with_data: u32 = 0;

    let mut offset = 0usize;
    for &clen in contig_lengths {
        let end = (offset + clen as usize).min(mag_values.len());
        let slice = &mag_values[offset..end];
        if !slice.is_empty() {
            let has_data = slice.iter().any(|&v| v > 0);
            if has_data { n_contigs_with_data += 1; }
            let mut prev: Option<i32> = None;
            for &v in slice {
                let vv = v.max(0);
                sum_x += vv as u128;
                sum_x2 += (vv as u128) * (vv as u128);
                if vv >= 1 { n_covered += 1; }
                if let Some(p) = prev {
                    let d = (vv as i64 - p as i64) as i128;
                    sum_sq_diff += (d * d) as u128;
                }
                prev = Some(vv);
            }
        }
        offset = end;
    }

    let n_f = n_positions as f64;
    let mean_raw = sum_x as f64 / n_f;
    let af_raw = (n_covered as f64 / n_f) * 100.0;
    let expected_af_raw = if mean_raw > 0.0 {
        (1.0 - (-0.883 * mean_raw).exp()) * 100.0
    } else { 0.0 };

    let mut sorted = mag_values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    let median_raw = sorted[n / 2] as f64;
    let lo = (n as f64 * 0.05).floor() as usize;
    let hi = ((n as f64 * 0.95).ceil() as usize).min(n);
    let trimmed_raw = if hi > lo {
        sorted[lo..hi].iter().map(|&v| v as f64).sum::<f64>() / (hi - lo) as f64
    } else { 0.0 };

    let variance = ((sum_x2 as f64 / n_f) - (mean_raw * mean_raw)).max(0.0);
    let cv_raw = if mean_raw > 0.0 { variance.sqrt() / mean_raw } else { 0.0 };
    let rough_denom = (n_positions as i64) - (n_contigs_with_data as i64);
    let roughness_raw = if mean_raw > 0.0 && rough_denom > 0 {
        (sum_sq_diff as f64 / rough_denom as f64).sqrt() / mean_raw
    } else { 0.0 };

    let s_m = crate::types::get_column_scale("Coverage", "Coverage_mean");
    let s_af = crate::types::get_column_scale("Coverage", "Aligned_fraction_percentage");
    let s_cv = crate::types::get_column_scale("Coverage", "Coverage_coefficient_of_variation");
    let s_rcr = crate::types::get_column_scale("Coverage", "Coverage_relative_coverage_roughness");
    MagCoverageStats {
        mag_length,
        aligned_fraction_pct_x10: (af_raw * s_af).round() as i32,
        expected_aligned_fraction_pct_x10: (expected_af_raw * s_af).round() as i32,
        read_count: total_read_count,
        coverage_mean_x10: (mean_raw * s_m).round() as i32,
        coverage_median_x10: (median_raw * s_m).round() as i32,
        coverage_trimmed_mean_x10: (trimmed_raw * s_m).round() as i32,
        coverage_cv_x1e6: (cv_raw * s_cv).round() as i32,
        coverage_roughness_x1e6: (roughness_raw * s_rcr).round() as i32,
    }
}
