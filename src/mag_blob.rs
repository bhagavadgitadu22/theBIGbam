//! MAG-scale blob aggregation.
//!
//! Aggregates per-contig `Feature_blob` / `Contig_blob` data into a single
//! MAG-wide blob pair (zoom + chunks) and persists it to the `MAG_blob` /
//! `MAG_contig_blob` tables. Downstream the Python plotting layer fetches
//! one row per subplot instead of stitching per-contig blobs together —
//! preserving the DNAFeaturesViewer y_range in the gene map and guaranteeing
//! every member contig contributes at least one zoom bin regardless of length.
//!
//! Entry points:
//! - `build_mag_contig_blobs_all` — sample-independent. Runs once in
//!   pre-sample setup, after every Contig_blob source has been written.
//! - `build_mag_member_lookup` — builds contig→MAG mapping for the
//!   processing side to group jobs by MAG.
//! - `aggregate_feature_dense_inmem` / `aggregate_feature_sparse_inmem` —
//!   used by `process_mag_group` (in `processing.rs`) for remaining features.
//! - `compute_mag_coverage_stats_from_raw` — MAG coverage stats from raw
//!   accumulated per-position values (no decode needed).

use std::collections::HashMap;

use anyhow::Result;
use duckdb::{params, Connection};

use crate::blob::{
    decode_dense_from_chunks, decode_sparse_from_chunks, encode_dense_blob, encode_sparse_blob,
    EncodedBlob, EventMeta, MetadataFlags,
};
use crate::db::DbWriter;
use crate::gc_content::{DEFAULT_GC_CONTENT_WINDOW_SIZE, DEFAULT_GC_SKEW_WINDOW_SIZE};
use crate::types::{feature_name_to_id, get_encoding, Encoding, ValueScale, VARIABLES};

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

/// Per-feature window size + partner semantics for features stored in
/// `Contig_blob` (i.e. sample-independent contig-scoped data).
///
/// The Dense/Sparse routing itself comes from `types::get_encoding`, so
/// `VARIABLES` stays authoritative. `window_size == 1` means one value per
/// base-pair (the common case for sparse events and for the repeat RLE
/// re-expansion).
struct ContigBlobConfig {
    window_size: u32,
    partner_kind: PartnerKind,
}

fn contig_blob_config(name: &str) -> Option<ContigBlobConfig> {
    match name {
        "gc_content" => Some(ContigBlobConfig { window_size: DEFAULT_GC_CONTENT_WINDOW_SIZE as u32, partner_kind: PartnerKind::None }),
        "gc_skew"    => Some(ContigBlobConfig { window_size: DEFAULT_GC_SKEW_WINDOW_SIZE as u32, partner_kind: PartnerKind::None }),
        "direct_repeat_count"      => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "inverted_repeat_count"    => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "direct_repeat_identity"   => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::Position }),
        "inverted_repeat_identity" => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::Position }),
        "hit_count_within_mag"     => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "hit_identity_within_mag"  => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::ContigId }),
        _ => None,
    }
}

pub fn contig_blob_config_pub(name: &str) -> Option<()> {
    contig_blob_config(name).map(|_| ())
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

/// Sample-independent aggregation: Contig_blob → MAG_contig_blob.
/// Runs once in pre-sample setup, after GC content / GC skew / repeat blob
/// conversion / MAG hit features have all written their Contig_blob rows.
pub fn build_mag_contig_blobs_all(db_writer: &DbWriter) -> Result<()> {
    if !db_writer.is_mag_mode() {
        return Ok(());
    }
    let conn = db_writer.lock_conn()?;
    let mags = load_mag_members(&conn)?;
    if mags.is_empty() {
        return Ok(());
    }
    for (mag_id, members) in &mags {
        let mag_length: u32 = members.iter().map(|m| m.length).sum();
        build_mag_contig_blobs(db_writer, &conn, *mag_id, members, mag_length)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

/// Returns `[(MAG_id, [MagMember ordered by Offset_in_MAG])]`.
fn load_mag_members(conn: &Connection) -> Result<Vec<(i64, Vec<MagMember>)>> {
    // Ensure MAG table exists — non-MAG DBs won't have it.
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
// Sample-feature blobs (Feature_blob → MAG_blob) — in-memory aggregators
// ---------------------------------------------------------------------------

/// Convert `EncodedBlob.chunks` (`Vec<Vec<u8>>`, indexed by chunk_idx) into the
/// `(chunk_idx, data)` shape the existing chunk decoders expect.
fn chunks_with_idx(blob: &EncodedBlob) -> Vec<(i32, Vec<u8>)> {
    blob.chunks
        .iter()
        .enumerate()
        .map(|(i, data)| (i as i32, data.clone()))
        .collect()
}

pub fn aggregate_feature_dense_inmem(
    members_and_blobs: &[(MagMember, &EncodedBlob)],
    mag_length: u32,
    scale: ValueScale,
) -> EncodedBlob {
    let mut mag_values = vec![0i32; mag_length as usize];
    let mut any = false;

    for (m, blob) in members_and_blobs {
        if blob.chunks.is_empty() { continue; }
        let chunks = chunks_with_idx(blob);
        let values = decode_dense_from_chunks(&chunks, m.length);
        let off = m.offset as usize;
        let end = (off + m.length as usize).min(mag_values.len());
        let copy_n = end - off;
        if copy_n > 0 && copy_n <= values.len() {
            mag_values[off..end].copy_from_slice(&values[..copy_n]);
        }
        any = true;
    }

    if !any {
        return empty_blob();
    }
    let mut encoded = encode_dense_blob(&mag_values, scale, mag_length);
    encoded.chunks = Vec::new(); // MAG chunks not stored — assembled on the fly from per-contig chunks
    encoded
}

pub fn aggregate_feature_sparse_inmem(
    members_and_blobs: &[(MagMember, &EncodedBlob)],
    mag_length: u32,
    scale: ValueScale,
) -> EncodedBlob {
    let mut all_pos: Vec<u32> = Vec::new();
    let mut all_val: Vec<i32> = Vec::new();
    let mut all_meta: Vec<EventMeta> = Vec::new();
    let mut union_flags = MetadataFlags { sparse: true, ..Default::default() };

    for (m, blob) in members_and_blobs {
        if blob.chunks.is_empty() { continue; }
        let chunks = chunks_with_idx(blob);
        let (positions, values, meta, flags) = decode_sparse_from_chunks(&chunks);
        union_flags.has_stats    |= flags.has_stats;
        union_flags.has_sequence |= flags.has_sequence;
        union_flags.has_codons   |= flags.has_codons;
        union_flags.has_partner  |= flags.has_partner;
        for (i, p) in positions.iter().enumerate() {
            all_pos.push(p + m.offset);
            all_val.push(values[i]);
            all_meta.push(meta.get(i).cloned().unwrap_or_default());
        }
    }

    if all_pos.is_empty() {
        return empty_blob();
    }
    // Sort by position — zoom level computation is order-independent, but
    // encode_sparse_blob builds the zoom from sorted positions internally.
    let mut order: Vec<usize> = (0..all_pos.len()).collect();
    order.sort_unstable_by_key(|&i| all_pos[i]);
    let all_pos: Vec<u32> = order.iter().map(|&i| all_pos[i]).collect();
    let all_val: Vec<i32> = order.iter().map(|&i| all_val[i]).collect();
    let all_meta: Vec<EventMeta> = order.iter().map(|&i| all_meta[i].clone()).collect();
    let mut encoded = encode_sparse_blob(
        &all_pos, &all_val,
        Some(&all_meta), union_flags, scale, mag_length,
    );
    encoded.chunks = Vec::new(); // MAG chunks not stored — assembled on the fly from per-contig chunks
    encoded
}

// ---------------------------------------------------------------------------
// Contig-level blobs (Contig_blob → MAG_contig_blob)
// ---------------------------------------------------------------------------

fn build_mag_contig_blobs(
    db_writer: &DbWriter,
    conn: &Connection,
    mag_id: i64,
    members: &[MagMember],
    mag_length: u32,
) -> Result<()> {
    let mut blobs: Vec<(String, EncodedBlob)> = Vec::new();

    for var in VARIABLES.iter() {
        let cfg = match contig_blob_config(var.name) {
            Some(c) => c,
            None => continue,
        };
        let fid = match feature_name_to_id(var.name) { Some(v) => v, None => continue };

        let encoded = match get_encoding(var.name) {
            Encoding::Dense => aggregate_contig_dense(
                conn, fid, members, mag_length, var.value_scale, cfg.window_size,
            )?,
            Encoding::Sparse => aggregate_contig_sparse(
                conn, fid, members, mag_length, var.value_scale, cfg.partner_kind,
            )?,
        };

        if !encoded.zoom.is_empty() {
            blobs.push((var.name.to_string(), encoded));
        }
    }

    if !blobs.is_empty() {
        db_writer.write_mag_contig_blobs(conn, mag_id, &blobs)?;
    }
    Ok(())
}

fn aggregate_contig_dense(
    conn: &Connection,
    fid: i16,
    members: &[MagMember],
    mag_length: u32,
    scale: ValueScale,
    window_size: u32,
) -> Result<EncodedBlob> {
    let mut mag_values = vec![0i32; mag_length as usize];
    let mut any = false;

    for m in members {
        let chunks = load_contig_chunks(conn, m.contig_id, fid)?;
        if chunks.is_empty() { continue; }

        // For windowed features the source array has ceil(contig_length/window_size)
        // values; for per-bp features it's contig_length.
        let num_values = if window_size <= 1 {
            m.length
        } else {
            m.length.div_ceil(window_size)
        };
        let values = decode_dense_from_chunks(&chunks, num_values);

        // Expand to per-bp at MAG offset (no-op when window_size==1).
        let off = m.offset as usize;
        let contig_end = (off + m.length as usize).min(mag_values.len());
        if window_size <= 1 {
            let copy_n = contig_end - off;
            mag_values[off..contig_end].copy_from_slice(&values[..copy_n.min(values.len())]);
        } else {
            let w = window_size as usize;
            for (i, v) in values.iter().enumerate() {
                let s = off + i * w;
                if s >= contig_end { break; }
                let e = (s + w).min(contig_end);
                for slot in &mut mag_values[s..e] { *slot = *v; }
            }
        }
        any = true;
    }

    if !any {
        return Ok(empty_blob());
    }
    let mut encoded = encode_dense_blob(&mag_values, scale, mag_length);
    encoded.chunks = Vec::new(); // MAG chunks not stored — assembled on the fly from per-contig chunks
    Ok(encoded)
}

fn aggregate_contig_sparse(
    conn: &Connection,
    fid: i16,
    members: &[MagMember],
    mag_length: u32,
    scale: ValueScale,
    partner_kind: PartnerKind,
) -> Result<EncodedBlob> {
    let mut all_pos: Vec<u32> = Vec::new();
    let mut all_val: Vec<i32> = Vec::new();
    let mut all_meta: Vec<EventMeta> = Vec::new();
    let mut union_flags = MetadataFlags { sparse: true, ..Default::default() };

    for m in members {
        let chunks = load_contig_chunks(conn, m.contig_id, fid)?;
        if chunks.is_empty() { continue; }
        let (positions, values, meta, flags) = decode_sparse_from_chunks(&chunks);
        union_flags.has_stats    |= flags.has_stats;
        union_flags.has_sequence |= flags.has_sequence;
        union_flags.has_codons   |= flags.has_codons;
        union_flags.has_partner  |= flags.has_partner;
        for (i, p) in positions.iter().enumerate() {
            all_pos.push(p + m.offset);
            all_val.push(values[i]);
            let mut em = meta.get(i).cloned().unwrap_or_default();
            // Partner positions are contig-local and must be shifted into
            // MAG coordinates. ContigId partners are global ids — leave as-is.
            if partner_kind == PartnerKind::Position {
                if let Some(p) = em.partner {
                    if p >= 0 {
                        em.partner = Some(p + m.offset as i32);
                    }
                }
            }
            all_meta.push(em);
        }
    }

    if all_pos.is_empty() {
        return Ok(empty_blob());
    }
    // Sort by position — encode_sparse_blob uses binary search (partition_point)
    // to assign events to chunks, which requires sorted input.
    let mut order: Vec<usize> = (0..all_pos.len()).collect();
    order.sort_unstable_by_key(|&i| all_pos[i]);
    let all_pos: Vec<u32> = order.iter().map(|&i| all_pos[i]).collect();
    let all_val: Vec<i32> = order.iter().map(|&i| all_val[i]).collect();
    let all_meta: Vec<EventMeta> = order.iter().map(|&i| all_meta[i].clone()).collect();
    let mut encoded = encode_sparse_blob(
        &all_pos, &all_val,
        Some(&all_meta), union_flags, scale, mag_length,
    );
    encoded.chunks = Vec::new();
    Ok(encoded)
}

fn load_contig_chunks(
    conn: &Connection,
    contig_id: i64,
    fid: i16,
) -> Result<Vec<(i32, Vec<u8>)>> {
    let mut stmt = conn.prepare(
        "SELECT Chunk_idx, Data FROM Contig_blob_chunk
         WHERE Contig_id = ? AND Feature_id = ?
         ORDER BY Chunk_idx"
    )?;
    let rows = stmt.query_map(params![contig_id, fid], |r| {
        Ok((r.get::<_, i32>(0)?, r.get::<_, Vec<u8>>(1)?))
    })?;
    Ok(rows.collect::<std::result::Result<Vec<_>, _>>()?)
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
    pub above_expected: bool,
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
            mag_length, aligned_fraction_pct_x10: 0, above_expected: false,
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
    let above_expected = mean_raw > 0.0
        && af_raw >= (1.0 - (-0.883 * mean_raw).exp()) * 100.0;

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

    MagCoverageStats {
        mag_length,
        aligned_fraction_pct_x10: (af_raw * 10.0).round() as i32,
        above_expected,
        read_count: total_read_count,
        coverage_mean_x10: (mean_raw * 10.0).round() as i32,
        coverage_median_x10: (median_raw * 10.0).round() as i32,
        coverage_trimmed_mean_x10: (trimmed_raw * 10.0).round() as i32,
        coverage_cv_x1e6: (cv_raw * 1_000_000.0).round() as i32,
        coverage_roughness_x1e6: (roughness_raw * 1_000_000.0).round() as i32,
    }
}

// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------

fn empty_blob() -> EncodedBlob {
    EncodedBlob { zoom: Vec::new(), chunks: Vec::new() }
}
