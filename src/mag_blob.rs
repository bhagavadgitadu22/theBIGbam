//! MAG-scale blob aggregation.
//!
//! Aggregates per-contig `Feature_blob` / `Contig_blob` data into a single
//! MAG-wide blob pair (zoom + chunks) and persists it to the `MAG_blob` /
//! `MAG_contig_blob` tables. Downstream the Python plotting layer fetches
//! one row per subplot instead of stitching per-contig blobs together —
//! preserving the DNAFeaturesViewer y_range in the gene map and guaranteeing
//! every member contig contributes at least one zoom bin regardless of length.
//!
//! Two entry points, both running inline with normal sample processing:
//! - `build_mag_contig_blobs_all` — sample-independent. Runs once in
//!   pre-sample setup, after every Contig_blob source has been written.
//! - `build_mag_blobs_for_sample` — sample-specific. Runs in the writer
//!   thread immediately after each sample's Feature_blob rows are written.
//!   Decodes from the in-memory `EncodedBlob.chunks` rather than re-querying.

use std::collections::HashMap;

use anyhow::Result;
use duckdb::{params, Connection};

use crate::blob::{
    decode_dense_from_chunks, decode_sparse_from_chunks, encode_dense_blob, encode_sparse_blob,
    EncodedBlob, EventMeta, MetadataFlags,
};
use crate::db::DbWriter;
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
        "gc_content" => Some(ContigBlobConfig { window_size:  500, partner_kind: PartnerKind::None }),
        "gc_skew"    => Some(ContigBlobConfig { window_size: 1000, partner_kind: PartnerKind::None }),
        "direct_repeat_count"      => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "inverted_repeat_count"    => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "direct_repeat_identity"   => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::Position }),
        "inverted_repeat_identity" => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::Position }),
        "hit_count_within_mag"     => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::None }),
        "hit_identity_within_mag"  => Some(ContigBlobConfig { window_size: 1, partner_kind: PartnerKind::ContigId }),
        _ => None,
    }
}

pub struct MagMember {
    pub contig_id: i64,
    pub offset: u32,
    pub length: u32,
}

/// Per-contig MAG attribution, built once in the writer thread and shared
/// across every per-sample MAG-aggregation call. `mag_length` is duplicated
/// across each contig entry of a MAG for O(1) lookup without a second map.
#[derive(Clone, Default)]
pub struct MagLookup {
    /// contig_id → (mag_id, offset_in_mag, mag_length, contig_length)
    by_contig: HashMap<i64, (i64, u32, u32, u32)>,
    /// contig_name → contig_id
    name_to_id: HashMap<String, i64>,
}

impl MagLookup {
    pub fn is_empty(&self) -> bool {
        self.by_contig.is_empty()
    }

    fn contig_length(&self, cid: i64) -> u32 {
        self.by_contig.get(&cid).map(|t| t.3).unwrap_or(0)
    }
}

/// Build the contig→MAG lookup once per writer-thread session. Safe to call
/// outside MAG mode (returns empty).
pub fn build_mag_lookup(db_writer: &DbWriter) -> Result<MagLookup> {
    if !db_writer.is_mag_mode() {
        return Ok(MagLookup::default());
    }
    let conn = db_writer.lock_conn()?;
    build_mag_lookup_from_conn(&conn)
}

fn build_mag_lookup_from_conn(conn: &Connection) -> Result<MagLookup> {
    let has_assoc: bool = conn.query_row(
        "SELECT COUNT(*) > 0 FROM information_schema.tables
         WHERE table_name = 'MAG_contigs_association'",
        [],
        |r| r.get(0),
    ).unwrap_or(false);
    if !has_assoc {
        return Ok(MagLookup::default());
    }

    let mut stmt = conn.prepare(
        "SELECT mca.MAG_id, mca.Contig_id, mca.Offset_in_MAG,
                c.Contig_length, c.Contig_name,
                mag_len.MAG_length
         FROM MAG_contigs_association mca
         JOIN Contig c ON c.Contig_id = mca.Contig_id
         JOIN (
             SELECT MAG_id, SUM(Contig_length) AS MAG_length
             FROM MAG_contigs_association mca2
             JOIN Contig c2 ON c2.Contig_id = mca2.Contig_id
             GROUP BY MAG_id
         ) mag_len ON mag_len.MAG_id = mca.MAG_id"
    )?;
    let rows = stmt.query_map([], |r| Ok((
        r.get::<_, i64>(0)?,
        r.get::<_, i64>(1)?,
        r.get::<_, i64>(2)?,
        r.get::<_, i64>(3)?,
        r.get::<_, String>(4)?,
        r.get::<_, i64>(5)?,
    )))?;

    let mut out = MagLookup::default();
    for row in rows {
        let (mag_id, cid, off, clen, cname, mlen) = row?;
        out.by_contig.insert(cid, (mag_id, off as u32, mlen as u32, clen as u32));
        out.name_to_id.insert(cname, cid);
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

/// Sample-specific aggregation: Feature_blob → MAG_blob, scoped to one
/// sample's just-computed in-memory blobs. Decodes chunks straight from
/// `feature_blobs[i].2.chunks` (no DuckDB round-trip).
pub fn build_mag_blobs_for_sample(
    db_writer: &DbWriter,
    sample_id: i64,
    feature_blobs: &[(String, String, EncodedBlob)],
    lookup: &MagLookup,
) -> Result<()> {
    if !db_writer.is_mag_mode() || lookup.is_empty() || feature_blobs.is_empty() {
        return Ok(());
    }

    // Group (feature_name → mag_id → Vec<(member, &EncodedBlob)>).
    // We materialize MagMember per group so the existing aggregator helpers
    // work unchanged.
    type FeatureGroups<'a> = HashMap<String, HashMap<i64, Vec<(MagMember, &'a EncodedBlob)>>>;
    let mut groups: FeatureGroups = HashMap::new();
    let mut mag_lengths: HashMap<i64, u32> = HashMap::new();

    for (feature_name, contig_name, blob) in feature_blobs {
        // Skip Contig_blob features — handled by build_mag_contig_blobs_all.
        if contig_blob_config(feature_name).is_some() {
            continue;
        }
        let cid = match lookup.name_to_id.get(contig_name) {
            Some(v) => *v,
            None => continue,
        };
        let (mag_id, offset, mag_length, _clen) = match lookup.by_contig.get(&cid) {
            Some(v) => *v,
            None => continue,
        };
        mag_lengths.insert(mag_id, mag_length);
        groups
            .entry(feature_name.clone())
            .or_insert_with(HashMap::new)
            .entry(mag_id)
            .or_insert_with(Vec::new)
            .push((MagMember { contig_id: cid, offset, length: lookup.contig_length(cid) }, blob));
    }

    // Encode + write per (mag_id, feature). Group writes per (mag_id) so we
    // call write_mag_blobs once per MAG with the full feature batch.
    let mut by_mag: HashMap<i64, Vec<(String, EncodedBlob)>> = HashMap::new();
    // Side-channel: MAG-level coverage stats, recomputed from raw per-position
    // values while we're already decoding the `coverage` feature.
    let mut mag_cov_stats: HashMap<i64, MagCoverageStats> = HashMap::new();

    for (feature_name, per_mag) in groups {
        // Variable lookup for encoding metadata.
        let var = match VARIABLES.iter().find(|v| v.name == feature_name) {
            Some(v) => v,
            None => continue,
        };
        for (mag_id, members_and_blobs) in per_mag {
            let mag_length = *mag_lengths.get(&mag_id).unwrap_or(&0);
            if mag_length == 0 {
                continue;
            }
            if feature_name == "coverage" {
                mag_cov_stats.insert(
                    mag_id,
                    compute_mag_coverage_stats(&members_and_blobs, mag_length),
                );
            }
            let encoded = match var.encoding {
                Encoding::Dense => aggregate_feature_dense_inmem(
                    &members_and_blobs, mag_length, var.value_scale,
                ),
                Encoding::Sparse => aggregate_feature_sparse_inmem(
                    &members_and_blobs, mag_length, var.value_scale,
                ),
            };
            if !encoded.zoom.is_empty() {
                by_mag.entry(mag_id).or_insert_with(Vec::new)
                    .push((feature_name.clone(), encoded));
            }
        }
    }

    if by_mag.is_empty() && mag_cov_stats.is_empty() {
        return Ok(());
    }

    let conn = db_writer.lock_conn()?;
    for (mag_id, blobs) in &by_mag {
        if !blobs.is_empty() {
            db_writer.write_mag_blobs(&conn, *mag_id, sample_id, blobs)?;
        }
    }

    // Pull Read_count for every contig in this sample once, then sum per MAG.
    if !mag_cov_stats.is_empty() {
        let read_counts = load_sample_read_counts(&conn, sample_id)?;
        for (mag_id, stats) in mag_cov_stats.iter_mut() {
            let mut total: u64 = 0;
            for (cid, (mid, _, _, _)) in lookup.by_contig.iter() {
                if *mid == *mag_id {
                    if let Some(rc) = read_counts.get(cid) {
                        total += *rc;
                    }
                }
            }
            stats.read_count = total;
            db_writer.write_mag_coverage(&conn, *mag_id, sample_id, stats)?;
        }
    }
    Ok(())
}

fn load_sample_read_counts(
    conn: &Connection,
    sample_id: i64,
) -> Result<HashMap<i64, u64>> {
    let mut stmt = conn.prepare(
        "SELECT Contig_id, Read_count FROM Coverage WHERE Sample_id = ?"
    )?;
    let rows = stmt.query_map(params![sample_id], |r| Ok((
        r.get::<_, i64>(0)?,
        r.get::<_, i64>(1)?,
    )))?;
    let mut out: HashMap<i64, u64> = HashMap::new();
    for row in rows {
        let (cid, rc) = row?;
        out.insert(cid, rc.max(0) as u64);
    }
    Ok(out)
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

fn aggregate_feature_dense_inmem(
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
    encode_dense_blob(&mag_values, scale, mag_length)
}

fn aggregate_feature_sparse_inmem(
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
    encode_sparse_blob(
        &all_pos, &all_val,
        Some(&all_meta), union_flags, scale, mag_length,
    )
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
            (m.length + window_size - 1) / window_size
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
    Ok(encode_dense_blob(&mag_values, scale, mag_length))
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
    Ok(encode_sparse_blob(
        &all_pos, &all_val,
        Some(&all_meta), union_flags, scale, mag_length,
    ))
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

/// Compute MAG-level coverage metrics from per-member in-memory blobs.
///
/// Cross-contig boundary jumps are excluded from Σ(Δ²) by accumulating per
/// contig. Median and trimmed mean are true order statistics computed over
/// the concatenated per-position array (contig ordering does not affect
/// order statistics).
fn compute_mag_coverage_stats(
    members_and_blobs: &[(MagMember, &EncodedBlob)],
    mag_length: u32,
) -> MagCoverageStats {
    let mut sum_x: u128 = 0;
    let mut sum_x2: u128 = 0;
    let mut sum_sq_diff: u128 = 0;
    let mut n_positions: u64 = 0;
    let mut n_covered: u64 = 0;
    let mut n_contigs_with_data: u32 = 0;
    let mut all_values: Vec<i32> = Vec::new();

    for (m, blob) in members_and_blobs {
        if blob.chunks.is_empty() || m.length == 0 { continue; }
        let chunks = chunks_with_idx(blob);
        let values = decode_dense_from_chunks(&chunks, m.length);
        if values.is_empty() { continue; }
        n_contigs_with_data += 1;
        n_positions += values.len() as u64;

        // Single pass: Σx, Σx², n_covered; Σ(Δx)² over windows(2) stays
        // within this contig so junctions don't contaminate.
        let mut prev: Option<i32> = None;
        for &v in &values {
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

        all_values.extend_from_slice(&values);
    }

    if n_positions == 0 || all_values.is_empty() {
        return MagCoverageStats {
            mag_length,
            aligned_fraction_pct_x10: 0,
            above_expected: false,
            read_count: 0,
            coverage_mean_x10: 0,
            coverage_median_x10: 0,
            coverage_trimmed_mean_x10: 0,
            coverage_cv_x1e6: 0,
            coverage_roughness_x1e6: 0,
        };
    }

    let n_f = n_positions as f64;
    let mean_raw = sum_x as f64 / n_f;
    let af_raw = (n_covered as f64 / n_f) * 100.0;
    let above_expected = mean_raw > 0.0
        && af_raw >= (1.0 - (-0.883 * mean_raw).exp()) * 100.0;

    // Sort once, then read median + trimmed-mean slice as order statistics.
    // n is bounded by the MAG length so this is cheap per (MAG, sample).
    all_values.sort_unstable();
    let n = all_values.len();
    let median_raw = all_values[n / 2] as f64;

    let lo = (n as f64 * 0.05).floor() as usize;
    let hi = ((n as f64 * 0.95).ceil() as usize).min(n);
    let trimmed_raw = if hi > lo {
        let slice = &all_values[lo..hi];
        slice.iter().map(|&v| v as f64).sum::<f64>() / slice.len() as f64
    } else {
        0.0
    };

    let variance = (sum_x2 as f64 / n_f) - (mean_raw * mean_raw);
    let variance = variance.max(0.0);
    let cv_raw = if mean_raw > 0.0 { variance.sqrt() / mean_raw } else { 0.0 };

    let rough_denom = (n_positions as i64) - (n_contigs_with_data as i64);
    let roughness_raw = if mean_raw > 0.0 && rough_denom > 0 {
        (sum_sq_diff as f64 / rough_denom as f64).sqrt() / mean_raw
    } else {
        0.0
    };

    MagCoverageStats {
        mag_length,
        aligned_fraction_pct_x10: (af_raw * 10.0).round() as i32,
        above_expected,
        read_count: 0, // filled in by caller after bulk Coverage lookup
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
    EncodedBlob { data: Vec::new(), zoom: Vec::new(), chunks: Vec::new() }
}
