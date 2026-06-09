//! MAG-scale blob aggregation.
//!
//! Builds MAG-wide zoom-only BLOBs and persists them to the `MAG_blob` /
//! `MAG_contig_blob` tables. Downstream the Python plotting layer fetches
//! one row per subplot instead of stitching per-contig blobs together —
//! preserving the DNAFeaturesViewer y_range in the gene map and guaranteeing
//! every member contig contributes at least one zoom bin regardless of length.
//!
//! Entry points:
//! - `build_mag_contig_blobs_from_raw` — sample-independent. Builds
//!   MAG_contig_blob rows directly from raw GC / repeat / hit data.
//! - `build_mag_member_lookup` — builds contig→MAG mapping for the
//!   processing side to group jobs by MAG.
//! - `compute_mag_coverage_stats_from_raw` — MAG coverage stats from raw
//!   accumulated per-position values (no decode needed).

use std::collections::HashMap;

use anyhow::Result;
use duckdb::Connection;

use crate::blob::{
    encode_dense_blob, encode_sparse_blob,
    EncodedBlob, EventMeta, MetadataFlags,
};
use crate::db::{DbWriter, GCContentData};
use crate::gc_content::{DEFAULT_GC_CONTENT_WINDOW_SIZE, DEFAULT_GC_SKEW_WINDOW_SIZE};
use crate::types::get_value_scale;

/// How `EventMeta.partner` is interpreted for a sparse Contig_blob feature.
/// Drives offset adjustment during MAG aggregation.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum PartnerKind {
    /// No partner metadata carried.
    None,
    /// partner is a global `Contig_id` — leave as-is during MAG aggregation.
    ContigId,
    /// partner is a contig-local position — shift by `MagMember.offset` so
    /// the partner position lands in MAG coordinates.
    Position,
}

fn contig_blob_partner_kind(name: &str) -> PartnerKind {
    match name {
        "direct_repeat_identity"   => PartnerKind::Position,
        "inverted_repeat_identity" => PartnerKind::Position,
        "hit_identity_within_mag"  => PartnerKind::ContigId,
        _ => PartnerKind::None,
    }
}

/// Raw per-contig sparse feature data returned by DB writers so MAG blobs
/// can be built directly without decode→re-encode.
pub struct ContigSparseRaw {
    pub feature_name: &'static str,
    pub contig_id: i64,
    pub contig_length: u32,
    pub positions: Vec<u32>,
    pub values: Vec<i32>,
    pub metadata: Option<Vec<EventMeta>>,
    pub flags: MetadataFlags,
}

pub struct MagMember {
    pub contig_id: i64,
    pub offset: u32,
    pub length: u32,
}

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

/// Build MAG_contig_blob rows directly from raw GC, repeat, and hit data.
/// Runs once in pre-sample setup, after per-contig Contig_blob rows have been
/// written. No decode→re-encode — uses raw arrays/sparse data directly.
pub fn build_mag_contig_blobs_from_raw(
    db_writer: &DbWriter,
    gc_data: &[GCContentData],
    sparse_raw: Vec<ContigSparseRaw>,
) -> Result<()> {
    if !db_writer.is_mag_mode() {
        return Ok(());
    }
    let conn = db_writer.lock_conn()?;
    let mags = load_mag_members(&conn)?;
    if mags.is_empty() {
        return Ok(());
    }

    // Build contig_name → contig_id map from DbWriter for GC data lookup.
    let contig_name_to_id = db_writer.contig_name_to_id_map();

    // Index sparse raw data by contig_id for fast lookup.
    let mut sparse_by_cid: HashMap<i64, Vec<&ContigSparseRaw>> = HashMap::new();
    for entry in &sparse_raw {
        sparse_by_cid.entry(entry.contig_id).or_default().push(entry);
    }

    // Index GC data by contig_id.
    let mut gc_by_cid: HashMap<i64, &GCContentData> = HashMap::new();
    for gd in gc_data {
        if let Some(&cid) = contig_name_to_id.get(&gd.contig_name) {
            gc_by_cid.insert(cid, gd);
        }
    }

    for (mag_id, members) in &mags {
        let mag_length: u32 = members.iter().map(|m| m.length).sum();
        let mut blobs: Vec<(String, EncodedBlob)> = Vec::new();

        // --- GC content (dense, windowed) ---
        build_gc_mag_blob(
            &gc_by_cid, members, mag_length,
            "gc_content", DEFAULT_GC_CONTENT_WINDOW_SIZE as u32,
            |gd| gd.gc_values.iter().map(|&v| v as i32).collect(),
            &mut blobs,
        );

        // --- GC skew (dense, windowed) ---
        build_gc_mag_blob(
            &gc_by_cid, members, mag_length,
            "gc_skew", DEFAULT_GC_SKEW_WINDOW_SIZE as u32,
            |gd| gd.skew_values.iter().map(|&v| v as i32).collect(),
            &mut blobs,
        );

        // --- Sparse contig features (repeats, hits) ---
        // Group per-MAG sparse data by feature_name.
        let mut sparse_by_feature: HashMap<&str, Vec<(&ContigSparseRaw, &MagMember)>> = HashMap::new();
        for m in members {
            if let Some(entries) = sparse_by_cid.get(&m.contig_id) {
                for entry in entries {
                    sparse_by_feature.entry(entry.feature_name).or_default().push((entry, m));
                }
            }
        }

        for (feature_name, entries) in sparse_by_feature {
            let partner_kind = contig_blob_partner_kind(feature_name);

            let mut all_pos: Vec<u32> = Vec::new();
            let mut all_val: Vec<i32> = Vec::new();
            let mut all_meta: Vec<EventMeta> = Vec::new();
            let mut union_flags = MetadataFlags { sparse: true, ..Default::default() };

            for (entry, m) in &entries {
                union_flags.has_stats    |= entry.flags.has_stats;
                union_flags.has_sequence |= entry.flags.has_sequence;
                union_flags.has_codons   |= entry.flags.has_codons;
                union_flags.has_partner  |= entry.flags.has_partner;
                let meta_slice = entry.metadata.as_deref();
                for (i, &p) in entry.positions.iter().enumerate() {
                    all_pos.push(p + m.offset);
                    all_val.push(entry.values[i]);
                    let mut em = meta_slice.and_then(|ms| ms.get(i)).cloned().unwrap_or_default();
                    if partner_kind == PartnerKind::Position {
                        if let Some(ref mut pv) = em.partner {
                            if *pv >= 0 { *pv += m.offset as i32; }
                        }
                    }
                    all_meta.push(em);
                }
            }

            if all_pos.is_empty() { continue; }
            let mut order: Vec<usize> = (0..all_pos.len()).collect();
            order.sort_unstable_by_key(|&i| all_pos[i]);
            let sorted_pos: Vec<u32> = order.iter().map(|&i| all_pos[i]).collect();
            let sorted_val: Vec<i32> = order.iter().map(|&i| all_val[i]).collect();
            let sorted_meta: Vec<EventMeta> = order.iter().map(|&i| all_meta[i].clone()).collect();

            let scale = get_value_scale(feature_name);
            let has_any_meta = union_flags.has_stats || union_flags.has_sequence
                || union_flags.has_codons || union_flags.has_partner;
            let meta_ref = if has_any_meta { Some(sorted_meta.as_slice()) } else { None };
            let mut encoded = encode_sparse_blob(
                &sorted_pos, &sorted_val, meta_ref, union_flags, scale, mag_length,
            );
            encoded.chunks = Vec::new();
            blobs.push((feature_name.to_string(), encoded));
        }

        if !blobs.is_empty() {
            db_writer.write_mag_contig_blobs(&conn, *mag_id, &blobs)?;
        }
    }
    Ok(())
}

fn build_gc_mag_blob(
    gc_by_cid: &HashMap<i64, &GCContentData>,
    members: &[MagMember],
    mag_length: u32,
    feature_name: &str,
    window_size: u32,
    extract_values: impl Fn(&GCContentData) -> Vec<i32>,
    blobs: &mut Vec<(String, EncodedBlob)>,
) {
    let mut mag_values = vec![0i32; mag_length as usize];
    let mut any = false;
    let w = window_size as usize;

    for m in members {
        let gd = match gc_by_cid.get(&m.contig_id) {
            Some(gd) => gd,
            None => continue,
        };
        let values = extract_values(gd);
        if values.is_empty() { continue; }

        let off = m.offset as usize;
        let contig_end = (off + m.length as usize).min(mag_values.len());
        for (i, &v) in values.iter().enumerate() {
            let s = off + i * w;
            if s >= contig_end { break; }
            let e = (s + w).min(contig_end);
            for slot in &mut mag_values[s..e] { *slot = v; }
        }
        any = true;
    }

    if !any { return; }
    let scale = get_value_scale(feature_name);
    let mut encoded = encode_dense_blob(&mag_values, scale, mag_length);
    encoded.chunks = Vec::new();
    blobs.push((feature_name.to_string(), encoded));
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

/// Returns `[(MAG_id, [MagMember ordered by Offset_in_MAG])]`.
fn load_mag_members(conn: &Connection) -> Result<Vec<(i64, Vec<MagMember>)>> {
    let has_mag_tables: bool = conn.query_row(
        "SELECT COUNT(*) > 0 FROM information_schema.tables
         WHERE table_name = 'MAG_contigs_association'",
        [],
        |r| r.get(0),
    ).unwrap_or(false);
    if !has_mag_tables {
        return Ok(Vec::new());
    }

    let mut stmt = conn.prepare(
        "SELECT mca.MAG_id, mca.Contig_id, mca.Offset_in_MAG, c.Contig_length
         FROM MAG_contigs_association mca
         JOIN Contig c ON c.Contig_id = mca.Contig_id
         ORDER BY mca.MAG_id, mca.Offset_in_MAG"
    )?;
    let rows = stmt.query_map([], |r| Ok((
        r.get::<_, i64>(0)?,
        r.get::<_, i64>(1)?,
        r.get::<_, i64>(2)?,
        r.get::<_, i64>(3)?,
    )))?;

    let mut out: Vec<(i64, Vec<MagMember>)> = Vec::new();
    for row in rows {
        let (mag_id, cid, off, len) = row?;
        if out.last().map(|(m, _)| *m != mag_id).unwrap_or(true) {
            out.push((mag_id, Vec::new()));
        }
        out.last_mut().unwrap().1.push(MagMember {
            contig_id: cid,
            offset: off as u32,
            length: len as u32,
        });
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

// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------
