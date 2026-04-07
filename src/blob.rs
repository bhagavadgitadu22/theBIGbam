//! Compressed BLOB encoding for genomic feature data.
//!
//! Replaces per-row RLE storage with BigWig-inspired compressed BLOBs.
//! Each BLOB stores one feature for one contig/sample combination.
//!
//! Binary format:
//! ```text
//! BLOB = Header (32 bytes) + BaseResolution + ZoomLevels[0..N] + SparseMetadata?
//! ```
//!
//! Key advantages:
//! - Dense BLOBs use array index as position (no position columns needed)
//! - Delta + zigzag + varint + zstd compression outperforms DuckDB column compression on RLE data
//! - Sparse metadata packed columnar within BLOB compresses better than per-row TEXT columns
//! - 1 table with ~3K rows vs 20 tables with ~1.7M rows

// ============================================================================
// Constants
// ============================================================================

/// Magic bytes identifying a theBIGbam feature BLOB
const MAGIC: &[u8; 4] = b"TBB\x01";

/// Current format version (simplified zoom levels: mean for dense, max for sparse + 10kbp zoom)
const FORMAT_VERSION: u8 = 1;

/// Positions per chunk for large genomes
const CHUNK_SIZE: u32 = 65536;

/// Zoom level bin sizes (100bp, 1000bp, 10000bp)
const ZOOM_BIN_SIZES: &[u32] = &[100, 1000, 10000];

/// Number of zoom levels
const NUM_ZOOM_LEVELS: u8 = 3;

/// Zoom level bin sizes for contig features (GC content, GC skew, repeats)
/// Only 10kbp zoom - used for contigs > 1 Mbp; base resolution used otherwise
const CONTIG_ZOOM_BIN_SIZES: &[u32] = &[10000];

/// Number of zoom levels for contig features
const NUM_CONTIG_ZOOM_LEVELS: u8 = 1;

/// Zstd compression level (3 = good balance of speed and compression)
const ZSTD_LEVEL: i32 = 3;

/// Magic bytes for standalone zoom BLOB
const ZOOM_MAGIC: &[u8; 4] = b"TBZ\x01";

// ============================================================================
// Types
// ============================================================================

/// Result of encoding a feature BLOB: the full data blob and a lightweight zoom-only blob.
/// The zoom blob can be fetched independently for large-window views, avoiding the cost
/// of transferring the full base-resolution data from the database.
pub struct EncodedBlob {
    /// Full BLOB (header + base resolution + zoom levels + metadata)
    pub data: Vec<u8>,
    /// Standalone zoom-only BLOB (small header + zoom levels only)
    pub zoom: Vec<u8>,
    /// Individual base-resolution chunks (one per CHUNK_SIZE bp region).
    /// Each entry is the raw zstd-compressed chunk data, stored as separate DB rows
    /// so Python can fetch only the 1-2 chunks overlapping a zoomed-in view.
    /// Empty for sparse BLOBs (sparse data is small enough to fetch whole).
    pub chunks: Vec<Vec<u8>>,
}

// ValueScale is defined in types.rs (single source of truth) and re-exported here for convenience.
pub use crate::types::ValueScale;

/// Flags indicating what metadata is present in a sparse BLOB.
#[derive(Clone, Copy, Debug, Default)]
pub struct MetadataFlags {
    /// BLOB contains sparse (event) data rather than dense per-position data
    pub sparse: bool,
    /// BLOB contains mean/median/std statistics per event
    pub has_stats: bool,
    /// BLOB contains sequence data per event
    pub has_sequence: bool,
    /// BLOB contains codon change data per event
    pub has_codons: bool,
}

impl MetadataFlags {
    fn to_byte(&self) -> u8 {
        let mut flags = 0u8;
        if self.sparse { flags |= 0x01; }
        if self.has_stats { flags |= 0x02; }
        if self.has_sequence { flags |= 0x04; }
        if self.has_codons { flags |= 0x08; }
        flags
    }

    #[allow(dead_code)]
    fn from_byte(b: u8) -> Self {
        Self {
            sparse: b & 0x01 != 0,
            has_stats: b & 0x02 != 0,
            has_sequence: b & 0x04 != 0,
            has_codons: b & 0x08 != 0,
        }
    }
}

/// Per-event metadata for sparse features.
#[derive(Clone, Debug, Default)]
pub struct EventMeta {
    /// Mean of event lengths (e.g., clip lengths) — scaled ×100
    pub mean: Option<i32>,
    /// Median of event lengths — scaled ×100
    pub median: Option<i32>,
    /// Standard deviation of event lengths — scaled ×100
    pub std: Option<i32>,
    /// Dominant sequence at this position (e.g., "A" for mismatch, "ACGT" for insertion)
    pub sequence: Option<Vec<u8>>,
    /// Prevalence of the dominant sequence (percentage × 10, e.g., 853 = 85.3%)
    pub prevalence: Option<i16>,
    /// Codon category: 0=Synonymous, 1=Non-synonymous, 2=Intergenic
    pub codon_category: Option<u8>,
    /// Codon ID: 0-63 index into standard 64-codon table (255=none)
    pub codon_id: Option<u8>,
    /// Amino acid ID: 0-20 index into amino acid list (255=none)
    pub aa_id: Option<u8>,
}

/// Dense zoom level bin - stores only mean value.
#[derive(Clone, Debug)]
struct ZoomBinDense {
    mean: i32,
}

/// Sparse zoom level bin - stores only max value.
#[derive(Clone, Debug)]
struct ZoomBinSparse {
    max_value: i32,
}

// ============================================================================
// Encoding Utilities
// ============================================================================

/// Delta-encode a slice of i32 values.
/// Returns the deltas: [v[0], v[1]-v[0], v[2]-v[1], ...]
fn delta_encode(values: &[i32]) -> Vec<i32> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut deltas = Vec::with_capacity(values.len());
    deltas.push(values[0]);
    for i in 1..values.len() {
        deltas.push(values[i] - values[i - 1]);
    }
    deltas
}

/// Zigzag-encode signed integers to unsigned (maps negatives to odds, positives to evens).
/// This makes small-magnitude values (positive or negative) encode as small unsigned values.
#[inline]
fn zigzag_encode_val(val: i32) -> u32 {
    ((val << 1) ^ (val >> 31)) as u32
}

/// Encode unsigned integers as variable-length bytes (LEB128).
fn varint_encode(values: &[u32]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(values.len() * 2);
    for &val in values {
        let mut v = val;
        loop {
            let mut byte = (v & 0x7F) as u8;
            v >>= 7;
            if v != 0 {
                byte |= 0x80;
            }
            buf.push(byte);
            if v == 0 {
                break;
            }
        }
    }
    buf
}

/// Compress bytes with zstd.
fn zstd_compress(data: &[u8]) -> Vec<u8> {
    zstd::encode_all(data.as_ref(), ZSTD_LEVEL).unwrap_or_else(|_| data.to_vec())
}

/// Write a little-endian u16 to a buffer.
fn write_u16(buf: &mut Vec<u8>, val: u16) {
    buf.extend_from_slice(&val.to_le_bytes());
}

/// Write a little-endian u32 to a buffer.
fn write_u32(buf: &mut Vec<u8>, val: u32) {
    buf.extend_from_slice(&val.to_le_bytes());
}

/// Write a little-endian i32 to a buffer.
fn write_i32(buf: &mut Vec<u8>, val: i32) {
    buf.extend_from_slice(&val.to_le_bytes());
}

/// Patch a u32 value at a specific offset in the buffer (little-endian).
fn patch_u32(buf: &mut [u8], offset: usize, val: u32) {
    buf[offset..offset + 4].copy_from_slice(&val.to_le_bytes());
}

// ============================================================================
// Zoom Level Computation
// ============================================================================

/// Compute zoom levels for dense (per-position) data.
/// Stores only the mean value per bin.
fn compute_zoom_levels_dense(values: &[i32], contig_length: u32, bin_sizes: &[u32]) -> Vec<Vec<ZoomBinDense>> {
    let mut all_levels = Vec::with_capacity(bin_sizes.len());

    for &bin_size in bin_sizes {
        let num_bins = (contig_length + bin_size - 1) / bin_size;
        let mut bins = Vec::with_capacity(num_bins as usize);

        for bin_idx in 0..num_bins {
            let start = (bin_idx * bin_size) as usize;
            let end = ((bin_idx + 1) * bin_size) as usize;
            let end = end.min(values.len());

            if start >= values.len() {
                bins.push(ZoomBinDense { mean: 0 });
                continue;
            }

            let slice = &values[start..end];
            let count = slice.len() as i64;
            let sum: i64 = slice.iter().map(|&v| v as i64).sum();
            let mean = if count > 0 { (sum / count) as i32 } else { 0 };

            bins.push(ZoomBinDense { mean });
        }

        all_levels.push(bins);
    }

    all_levels
}

/// Compute zoom levels for sparse (event) data.
/// Stores only the max value per bin (0 if no events in bin).
fn compute_zoom_levels_sparse(
    positions: &[u32],
    values: &[i32],
    contig_length: u32,
    bin_sizes: &[u32],
) -> Vec<Vec<ZoomBinSparse>> {
    let mut all_levels = Vec::with_capacity(bin_sizes.len());

    for &bin_size in bin_sizes {
        let num_bins = (contig_length + bin_size - 1) / bin_size;
        let mut bins: Vec<ZoomBinSparse> = (0..num_bins)
            .map(|_| ZoomBinSparse { max_value: 0 })
            .collect();

        for (&pos, &val) in positions.iter().zip(values.iter()) {
            let bin_idx = (pos / bin_size) as usize;
            if bin_idx < bins.len() {
                bins[bin_idx].max_value = bins[bin_idx].max_value.max(val);
            }
        }

        all_levels.push(bins);
    }

    all_levels
}

/// Encode dense zoom levels to bytes.
/// Each bin stores only the mean value (i32).
fn encode_zoom_levels_dense(levels: &[Vec<ZoomBinDense>]) -> Vec<u8> {
    let mut buf = Vec::new();

    for level in levels {
        let mut level_buf = Vec::new();
        for bin in level {
            write_i32(&mut level_buf, bin.mean);
        }
        let compressed = zstd_compress(&level_buf);
        write_u32(&mut buf, compressed.len() as u32);
        buf.extend_from_slice(&compressed);
    }

    buf
}

/// Encode sparse zoom levels to bytes.
/// Only nonzero bins are stored as (bin_index, max_value) pairs using the same
/// encoding pipeline as sparse base data: delta+zigzag+varint for indices,
/// zigzag+varint for values, then zstd compression.
fn encode_zoom_levels_sparse(levels: &[Vec<ZoomBinSparse>]) -> Vec<u8> {
    let mut buf = Vec::new();

    for level in levels {
        // Collect nonzero bins
        let mut indices: Vec<u32> = Vec::new();
        let mut values: Vec<i32> = Vec::new();
        for (i, bin) in level.iter().enumerate() {
            if bin.max_value != 0 {
                indices.push(i as u32);
                values.push(bin.max_value);
            }
        }

        let mut level_buf = Vec::new();
        // Write nonzero count
        write_u32(&mut level_buf, indices.len() as u32);

        if !indices.is_empty() {
            // Bin indices: delta + zigzag + varint
            let idx_deltas = delta_encode(&indices.iter().map(|&i| i as i32).collect::<Vec<_>>());
            let idx_zigzagged: Vec<u32> = idx_deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
            let idx_varint = varint_encode(&idx_zigzagged);
            write_u32(&mut level_buf, idx_varint.len() as u32);
            level_buf.extend_from_slice(&idx_varint);

            // Max values: zigzag + varint (no delta)
            let val_zigzagged: Vec<u32> = values.iter().map(|&v| zigzag_encode_val(v)).collect();
            let val_varint = varint_encode(&val_zigzagged);
            write_u32(&mut level_buf, val_varint.len() as u32);
            level_buf.extend_from_slice(&val_varint);
        }

        let compressed = zstd_compress(&level_buf);
        write_u32(&mut buf, compressed.len() as u32);
        buf.extend_from_slice(&compressed);
    }

    buf
}

// ============================================================================
// Standalone Zoom BLOB
// ============================================================================

/// Build a standalone zoom BLOB that can be fetched independently of base-resolution data.
///
/// Format: ZOOM_MAGIC (4) + scale_code (1) + flags (1) + num_zoom_levels (1) + reserved (1)
///       + contig_length (4) + reserved (4) = 16-byte header + zoom level data
fn build_zoom_blob(scale: ValueScale, sparse: bool, num_zoom_levels: u8, contig_length: u32, zoom_bytes: &[u8]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(16 + zoom_bytes.len());
    buf.extend_from_slice(ZOOM_MAGIC);          // [0..4] magic
    buf.push(scale as u8);                       // [4] scale_code
    let flags: u8 = if sparse { 0x01 } else { 0x00 };
    buf.push(flags);                             // [5] flags (sparse bit)
    buf.push(num_zoom_levels);                   // [6] num_zoom_levels
    buf.push(0);                                 // [7] reserved
    write_u32(&mut buf, contig_length);          // [8..12] contig_length
    write_u32(&mut buf, 0);                      // [12..16] reserved
    buf.extend_from_slice(zoom_bytes);           // zoom level data
    buf
}

// ============================================================================
// Dense BLOB Encoding
// ============================================================================

/// Encode a dense (per-position) feature as a compressed BLOB.
///
/// Dense features have a value at every position (e.g., coverage, MAPQ).
/// The array index IS the position, so no position column is needed.
///
/// # Arguments
/// * `values` - Per-position values (length = contig_length)
/// * `scale` - Scale factor applied to values
/// * `contig_length` - Length of the contig in base pairs
pub fn encode_dense_blob(values: &[i32], scale: ValueScale, contig_length: u32) -> EncodedBlob {
    // Skip encoding when all values are zero — caller filters empty blobs
    if values.iter().all(|&v| v == 0) {
        return EncodedBlob { data: Vec::new(), zoom: Vec::new(), chunks: Vec::new() };
    }
    let mut blob = Vec::with_capacity(32 + values.len() * 2);

    // === Header (32 bytes) ===
    blob.extend_from_slice(MAGIC);                    // [0..4] magic
    blob.push(FORMAT_VERSION);                         // [4] version
    blob.push(MetadataFlags::default().to_byte());     // [5] flags (dense, no metadata)
    blob.push(scale as u8);                            // [6] scale_code
    blob.push(NUM_ZOOM_LEVELS);                        // [7] num_zoom_levels
    write_u32(&mut blob, contig_length);               // [8..12] contig_length
    // Placeholders for offsets (patched later)
    let base_block_offset_pos = blob.len();
    write_u32(&mut blob, 0);                           // [12..16] base_block_offset
    let base_block_size_pos = blob.len();
    write_u32(&mut blob, 0);                           // [16..20] base_block_compressed_size
    let zoom_index_offset_pos = blob.len();
    write_u32(&mut blob, 0);                           // [20..24] zoom_index_offset
    write_u32(&mut blob, 0);                           // [24..28] sparse_meta_offset (0=none)
    write_u32(&mut blob, 0);                           // [28..32] sparse_meta_compressed_size (0=none)

    // === Base Resolution (chunked) ===
    let base_block_start = blob.len() as u32;
    patch_u32(&mut blob, base_block_offset_pos, base_block_start);

    let num_chunks = (values.len() as u32 + CHUNK_SIZE - 1) / CHUNK_SIZE;
    write_u16(&mut blob, num_chunks as u16);
    write_u32(&mut blob, CHUNK_SIZE);

    let mut chunks = Vec::with_capacity(num_chunks as usize);
    for chunk_idx in 0..num_chunks {
        let start = (chunk_idx * CHUNK_SIZE) as usize;
        let end = ((chunk_idx + 1) * CHUNK_SIZE) as usize;
        let end = end.min(values.len());
        let chunk_values = &values[start..end];

        // Delta → zigzag → varint → zstd
        let deltas = delta_encode(chunk_values);
        let zigzagged: Vec<u32> = deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
        let varint_bytes = varint_encode(&zigzagged);
        let compressed = zstd_compress(&varint_bytes);

        write_u32(&mut blob, compressed.len() as u32);
        blob.extend_from_slice(&compressed);
        chunks.push(compressed);
    }

    let base_block_end = blob.len() as u32;
    patch_u32(&mut blob, base_block_size_pos, base_block_end - base_block_start);

    // === Zoom Levels ===
    let zoom_offset = blob.len() as u32;
    patch_u32(&mut blob, zoom_index_offset_pos, zoom_offset);

    let zoom_levels = compute_zoom_levels_dense(values, contig_length, ZOOM_BIN_SIZES);
    let zoom_bytes = encode_zoom_levels_dense(&zoom_levels);
    blob.extend_from_slice(&zoom_bytes);

    let zoom_blob = build_zoom_blob(scale, false, NUM_ZOOM_LEVELS, contig_length, &zoom_bytes);
    EncodedBlob { data: blob, zoom: zoom_blob, chunks }
}

// ============================================================================
// Sparse BLOB Encoding
// ============================================================================

/// Encode a sparse (event-based) feature as a compressed BLOB.
///
/// Sparse features have values only at certain positions (e.g., mismatches, clippings).
/// Positions are stored explicitly using delta + varint encoding.
///
/// # Arguments
/// * `positions` - 0-indexed positions where events occur
/// * `values` - Event values (parallel to positions)
/// * `metadata` - Optional per-event metadata
/// * `flags` - Which metadata fields are present
/// * `scale` - Scale factor applied to values
/// * `contig_length` - Length of the contig in base pairs
pub fn encode_sparse_blob(
    positions: &[u32],
    values: &[i32],
    metadata: Option<&[EventMeta]>,
    flags: MetadataFlags,
    scale: ValueScale,
    contig_length: u32,
) -> EncodedBlob {
    assert_eq!(positions.len(), values.len());

    // Skip encoding when no events — caller filters empty blobs
    if positions.is_empty() {
        return EncodedBlob { data: Vec::new(), zoom: Vec::new(), chunks: Vec::new() };
    }

    let mut blob = Vec::with_capacity(32 + positions.len() * 8);

    // === Header (32 bytes) ===
    blob.extend_from_slice(MAGIC);                         // [0..4] magic
    blob.push(FORMAT_VERSION);                              // [4] version
    let mut meta_flags = flags;
    meta_flags.sparse = true;
    blob.push(meta_flags.to_byte());                        // [5] flags
    blob.push(scale as u8);                                 // [6] scale_code
    blob.push(NUM_ZOOM_LEVELS);                             // [7] num_zoom_levels
    write_u32(&mut blob, contig_length);                    // [8..12] contig_length
    let base_block_offset_pos = blob.len();
    write_u32(&mut blob, 0);                                // [12..16] base_block_offset
    let base_block_size_pos = blob.len();
    write_u32(&mut blob, 0);                                // [16..20] base_block_compressed_size
    let zoom_index_offset_pos = blob.len();
    write_u32(&mut blob, 0);                                // [20..24] zoom_index_offset
    let sparse_meta_offset_pos = blob.len();
    write_u32(&mut blob, 0);                                // [24..28] sparse_meta_offset
    let sparse_meta_size_pos = blob.len();
    write_u32(&mut blob, 0);                                // [28..32] sparse_meta_compressed_size

    // === Sparse Base Resolution ===
    let base_block_start = blob.len() as u32;
    patch_u32(&mut blob, base_block_offset_pos, base_block_start);

    // Event count
    write_u32(&mut blob, positions.len() as u32);

    // Delta-encode positions, then varint compress, then zstd
    let pos_deltas = delta_encode(&positions.iter().map(|&p| p as i32).collect::<Vec<_>>());
    let pos_zigzagged: Vec<u32> = pos_deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
    let pos_varint = varint_encode(&pos_zigzagged);
    let pos_compressed = zstd_compress(&pos_varint);
    write_u32(&mut blob, pos_compressed.len() as u32);
    blob.extend_from_slice(&pos_compressed);

    // Values: zigzag + varint + zstd (no delta — sparse values don't correlate well)
    let val_zigzagged: Vec<u32> = values.iter().map(|&v| zigzag_encode_val(v)).collect();
    let val_varint = varint_encode(&val_zigzagged);
    let val_compressed = zstd_compress(&val_varint);
    write_u32(&mut blob, val_compressed.len() as u32);
    blob.extend_from_slice(&val_compressed);

    let base_block_end = blob.len() as u32;
    patch_u32(&mut blob, base_block_size_pos, base_block_end - base_block_start);

    // === Zoom Levels ===
    let zoom_offset = blob.len() as u32;
    patch_u32(&mut blob, zoom_index_offset_pos, zoom_offset);

    let zoom_levels = compute_zoom_levels_sparse(positions, values, contig_length, ZOOM_BIN_SIZES);
    let zoom_bytes = encode_zoom_levels_sparse(&zoom_levels);
    blob.extend_from_slice(&zoom_bytes);

    // === Sparse Metadata Block ===
    if let Some(meta) = metadata {
        if meta.len() == positions.len() && (flags.has_stats || flags.has_sequence || flags.has_codons) {
            let sparse_meta_start = blob.len() as u32;
            patch_u32(&mut blob, sparse_meta_offset_pos, sparse_meta_start);

            let meta_bytes = encode_sparse_metadata(meta, &flags);
            let meta_compressed = zstd_compress(&meta_bytes);
            write_u32(&mut blob, meta_compressed.len() as u32);
            blob.extend_from_slice(&meta_compressed);

            let meta_total_size = (blob.len() as u32) - sparse_meta_start;
            patch_u32(&mut blob, sparse_meta_size_pos, meta_total_size);
        }
    }

    // === Build sparse chunks: group events by CHUNK_SIZE bp position ranges ===
    let has_meta = metadata.is_some()
        && metadata.unwrap().len() == positions.len()
        && (flags.has_stats || flags.has_sequence || flags.has_codons);
    let num_chunks = (contig_length + CHUNK_SIZE - 1) / CHUNK_SIZE;
    let mut chunks = Vec::with_capacity(num_chunks as usize);
    for chunk_idx in 0..num_chunks {
        let range_start = chunk_idx * CHUNK_SIZE;
        let range_end = range_start + CHUNK_SIZE;
        // Binary search for events in [range_start, range_end)
        let lo = positions.partition_point(|&p| p < range_start);
        let hi = positions.partition_point(|&p| p < range_end);
        if lo >= hi {
            chunks.push(Vec::new()); // empty chunk
            continue;
        }
        let chunk_pos = &positions[lo..hi];
        let chunk_vals = &values[lo..hi];
        let n_events = (hi - lo) as u32;

        let mut cb = Vec::new();
        write_u32(&mut cb, n_events);
        // Positions: delta → zigzag → varint → zstd
        let p_deltas = delta_encode(&chunk_pos.iter().map(|&p| p as i32).collect::<Vec<_>>());
        let p_zig: Vec<u32> = p_deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
        let p_compressed = zstd_compress(&varint_encode(&p_zig));
        write_u32(&mut cb, p_compressed.len() as u32);
        cb.extend_from_slice(&p_compressed);
        // Values: zigzag → varint → zstd
        let v_zig: Vec<u32> = chunk_vals.iter().map(|&v| zigzag_encode_val(v)).collect();
        let v_compressed = zstd_compress(&varint_encode(&v_zig));
        write_u32(&mut cb, v_compressed.len() as u32);
        cb.extend_from_slice(&v_compressed);
        // Metadata slice — store flags byte so decoder knows which fields are present
        if has_meta {
            let chunk_meta = &metadata.unwrap()[lo..hi];
            let meta_bytes = encode_sparse_metadata(chunk_meta, &flags);
            let meta_compressed = zstd_compress(&meta_bytes);
            cb.push(flags.to_byte()); // flags byte (has_stats, has_sequence, has_codons)
            write_u32(&mut cb, meta_compressed.len() as u32);
            cb.extend_from_slice(&meta_compressed);
        } else {
            cb.push(0u8); // no metadata
        }
        chunks.push(cb);
    }

    let zoom_blob = build_zoom_blob(scale, true, NUM_ZOOM_LEVELS, contig_length, &zoom_bytes);
    EncodedBlob { data: blob, zoom: zoom_blob, chunks }
}

/// Encode sparse metadata in columnar format.
/// Each field type is packed together for better compression.
fn encode_sparse_metadata(metadata: &[EventMeta], flags: &MetadataFlags) -> Vec<u8> {
    let mut buf = Vec::new();

    // Stats: mean[], median[], std[] — each as i32 little-endian
    if flags.has_stats {
        for m in metadata {
            write_i32(&mut buf, m.mean.unwrap_or(0));
        }
        for m in metadata {
            write_i32(&mut buf, m.median.unwrap_or(0));
        }
        for m in metadata {
            write_i32(&mut buf, m.std.unwrap_or(0));
        }
    }

    // Sequences: [len: u8, bytes: u8[], prevalence: i16] per event
    if flags.has_sequence {
        for m in metadata {
            if let Some(ref seq) = m.sequence {
                let len = seq.len().min(255) as u8;
                buf.push(len);
                buf.extend_from_slice(&seq[..len as usize]);
            } else {
                buf.push(0);
            }
        }
        // Prevalence as i16 array (columnar, after all sequences)
        for m in metadata {
            let prev = m.prevalence.unwrap_or(0);
            buf.extend_from_slice(&prev.to_le_bytes());
        }
    }

    // Codons: [category: u8, codon_id: u8, aa_id: u8] per event — 3 bytes total
    if flags.has_codons {
        for m in metadata {
            buf.push(m.codon_category.unwrap_or(255));
        }
        for m in metadata {
            buf.push(m.codon_id.unwrap_or(255));
        }
        for m in metadata {
            buf.push(m.aa_id.unwrap_or(255));
        }
    }

    buf
}

// ============================================================================
// Contig_blob-Specific Encoding (only 10kbp zoom)
// ============================================================================

/// Encode a dense contig feature (GC content, GC skew) as a compressed BLOB.
/// Uses only 10kbp zoom level for contigs > 1 Mbp.
///
/// `window_size` is the genomic window each value represents (e.g. 500 for GC content).
/// Zoom bins are computed in window-index space so indices align with the values array.
pub fn encode_contig_dense_blob(values: &[i32], scale: ValueScale, contig_length: u32, window_size: u32) -> EncodedBlob {
    let mut blob = Vec::with_capacity(32 + values.len() * 2);

    // === Header (32 bytes) ===
    blob.extend_from_slice(MAGIC);                    // [0..4] magic
    blob.push(FORMAT_VERSION);                         // [4] version
    blob.push(MetadataFlags::default().to_byte());     // [5] flags (dense, no metadata)
    blob.push(scale as u8);                            // [6] scale_code
    blob.push(NUM_CONTIG_ZOOM_LEVELS);                 // [7] num_zoom_levels (1 for contig)
    write_u32(&mut blob, contig_length);               // [8..12] contig_length
    // Placeholders for offsets (patched later)
    let base_block_offset_pos = blob.len();
    write_u32(&mut blob, 0);                           // [12..16] base_block_offset
    let base_block_size_pos = blob.len();
    write_u32(&mut blob, 0);                           // [16..20] base_block_compressed_size
    let zoom_index_offset_pos = blob.len();
    write_u32(&mut blob, 0);                           // [20..24] zoom_index_offset
    write_u32(&mut blob, 0);                           // [24..28] sparse_meta_offset (0=none)
    write_u32(&mut blob, 0);                           // [28..32] sparse_meta_compressed_size (0=none)

    // === Base Resolution (chunked) ===
    let base_block_start = blob.len() as u32;
    patch_u32(&mut blob, base_block_offset_pos, base_block_start);

    let num_chunks = (values.len() as u32 + CHUNK_SIZE - 1) / CHUNK_SIZE;
    write_u16(&mut blob, num_chunks as u16);
    write_u32(&mut blob, CHUNK_SIZE);

    let mut chunks = Vec::with_capacity(num_chunks as usize);
    for chunk_idx in 0..num_chunks {
        let start = (chunk_idx * CHUNK_SIZE) as usize;
        let end = ((chunk_idx + 1) * CHUNK_SIZE) as usize;
        let end = end.min(values.len());
        let chunk_values = &values[start..end];

        // Delta → zigzag → varint → zstd
        let deltas = delta_encode(chunk_values);
        let zigzagged: Vec<u32> = deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
        let varint_bytes = varint_encode(&zigzagged);
        let compressed = zstd_compress(&varint_bytes);

        write_u32(&mut blob, compressed.len() as u32);
        blob.extend_from_slice(&compressed);
        chunks.push(compressed);
    }

    let base_block_end = blob.len() as u32;
    patch_u32(&mut blob, base_block_size_pos, base_block_end - base_block_start);

    // === Zoom Levels (only 10kbp for contig features) ===
    let zoom_offset = blob.len() as u32;
    patch_u32(&mut blob, zoom_index_offset_pos, zoom_offset);

    // Scale bin sizes from bp to window indices so zoom computation indexes correctly
    let scaled_bins: Vec<u32> = CONTIG_ZOOM_BIN_SIZES.iter()
        .map(|&bs| (bs / window_size).max(1))
        .collect();
    let zoom_levels = compute_zoom_levels_dense(values, values.len() as u32, &scaled_bins);
    let zoom_bytes = encode_zoom_levels_dense(&zoom_levels);
    blob.extend_from_slice(&zoom_bytes);

    let zoom_blob = build_zoom_blob(scale, false, NUM_CONTIG_ZOOM_LEVELS, contig_length, &zoom_bytes);
    EncodedBlob { data: blob, zoom: zoom_blob, chunks }
}

/// Encode a sparse contig feature (repeat features) as a compressed BLOB.
/// Uses 100bp, 1000bp, 10000bp zoom levels (same as sample-level blobs).
pub fn encode_contig_sparse_blob(
    positions: &[u32],
    values: &[i32],
    scale: ValueScale,
    contig_length: u32,
) -> EncodedBlob {
    assert_eq!(positions.len(), values.len());

    let mut blob = Vec::with_capacity(32 + positions.len() * 8);

    // === Header (32 bytes) ===
    blob.extend_from_slice(MAGIC);                         // [0..4] magic
    blob.push(FORMAT_VERSION);                              // [4] version
    let mut meta_flags = MetadataFlags::default();
    meta_flags.sparse = true;
    blob.push(meta_flags.to_byte());                        // [5] flags
    blob.push(scale as u8);                                 // [6] scale_code
    blob.push(NUM_ZOOM_LEVELS);                              // [7] num_zoom_levels (3: 100bp, 1000bp, 10000bp)
    write_u32(&mut blob, contig_length);                    // [8..12] contig_length
    let base_block_offset_pos = blob.len();
    write_u32(&mut blob, 0);                                // [12..16] base_block_offset
    let base_block_size_pos = blob.len();
    write_u32(&mut blob, 0);                                // [16..20] base_block_compressed_size
    let zoom_index_offset_pos = blob.len();
    write_u32(&mut blob, 0);                                // [20..24] zoom_index_offset
    write_u32(&mut blob, 0);                                // [24..28] sparse_meta_offset (0=none)
    write_u32(&mut blob, 0);                                // [28..32] sparse_meta_compressed_size (0=none)

    // === Sparse Base Resolution ===
    let base_block_start = blob.len() as u32;
    patch_u32(&mut blob, base_block_offset_pos, base_block_start);

    // Event count
    write_u32(&mut blob, positions.len() as u32);

    // Delta-encode positions, then varint compress, then zstd
    let pos_deltas = delta_encode(&positions.iter().map(|&p| p as i32).collect::<Vec<_>>());
    let pos_zigzagged: Vec<u32> = pos_deltas.iter().map(|&d| zigzag_encode_val(d)).collect();
    let pos_varint = varint_encode(&pos_zigzagged);
    let pos_compressed = zstd_compress(&pos_varint);
    write_u32(&mut blob, pos_compressed.len() as u32);
    blob.extend_from_slice(&pos_compressed);

    // Values: zigzag + varint + zstd
    let val_zigzagged: Vec<u32> = values.iter().map(|&v| zigzag_encode_val(v)).collect();
    let val_varint = varint_encode(&val_zigzagged);
    let val_compressed = zstd_compress(&val_varint);
    write_u32(&mut blob, val_compressed.len() as u32);
    blob.extend_from_slice(&val_compressed);

    let base_block_end = blob.len() as u32;
    patch_u32(&mut blob, base_block_size_pos, base_block_end - base_block_start);

    // === Zoom Levels (100bp, 1000bp, 10000bp — same as sample-level blobs) ===
    let zoom_offset = blob.len() as u32;
    patch_u32(&mut blob, zoom_index_offset_pos, zoom_offset);

    let zoom_levels = compute_zoom_levels_sparse(positions, values, contig_length, ZOOM_BIN_SIZES);
    let zoom_bytes = encode_zoom_levels_sparse(&zoom_levels);
    blob.extend_from_slice(&zoom_bytes);

    let zoom_blob = build_zoom_blob(scale, true, NUM_ZOOM_LEVELS, contig_length, &zoom_bytes);
    EncodedBlob { data: blob, zoom: zoom_blob, chunks: Vec::new() }
}

/// Chunk size in base pairs (must match CHUNK_SIZE constant used in encoding).
pub const BLOB_CHUNK_SIZE: u32 = CHUNK_SIZE;

// ============================================================================
// Feature ID Mapping
// ============================================================================

/// Standard codon table for encoding codon_id (0-63).
pub const CODON_TABLE: &[&str] = &[
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
];

/// Amino acid table for encoding aa_id (0-20).
pub const AMINO_ACID_TABLE: &[&str] = &[
    "A (Alanine)", "C (Cysteine)", "D (Aspartic acid)", "E (Glutamic acid)",
    "F (Phenylalanine)", "G (Glycine)", "H (Histidine)", "I (Isoleucine)",
    "K (Lysine)", "L (Leucine)", "M (Methionine)", "N (Asparagine)",
    "P (Proline)", "Q (Glutamine)", "R (Arginine)", "S (Serine)",
    "T (Threonine)", "V (Valine)", "W (Tryptophan)", "Y (Tyrosine)",
    "* (Stop)",
];

/// Codon category encoding: text → u8.
pub fn codon_category_to_id(category: &str) -> u8 {
    match category {
        "Synonymous" => 0,
        "Non-synonymous" => 1,
        "Intergenic" => 2,
        _ => 255,
    }
}

/// Codon string → u8 id (0-63).
pub fn codon_to_id(codon: &str) -> u8 {
    CODON_TABLE.iter().position(|&c| c == codon).map(|p| p as u8).unwrap_or(255)
}

/// Amino acid description string → u8 id (0-20).
pub fn aa_to_id(aa: &str) -> u8 {
    AMINO_ACID_TABLE.iter().position(|&a| a == aa).map(|p| p as u8).unwrap_or(255)
}

/// Map feature name to a feature_id for the Feature_blob table.
/// Returns the 1-based index into VARIABLES, or None if not found.
pub fn feature_id(name: &str) -> Option<u16> {
    use crate::types::VARIABLES;
    VARIABLES.iter().position(|v| v.name == name).map(|i| (i + 1) as u16)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_encode() {
        let values = vec![10, 12, 15, 13, 20];
        let deltas = delta_encode(&values);
        assert_eq!(deltas, vec![10, 2, 3, -2, 7]);
    }

    #[test]
    fn test_zigzag_encode() {
        assert_eq!(zigzag_encode_val(0), 0);
        assert_eq!(zigzag_encode_val(-1), 1);
        assert_eq!(zigzag_encode_val(1), 2);
        assert_eq!(zigzag_encode_val(-2), 3);
        assert_eq!(zigzag_encode_val(2), 4);
    }

    #[test]
    fn test_varint_encode() {
        let values = vec![0, 1, 127, 128, 300];
        let encoded = varint_encode(&values);
        // 0 → [0x00], 1 → [0x01], 127 → [0x7F], 128 → [0x80, 0x01], 300 → [0xAC, 0x02]
        assert_eq!(encoded, vec![0x00, 0x01, 0x7F, 0x80, 0x01, 0xAC, 0x02]);
    }

    #[test]
    fn test_dense_blob_roundtrip_header() {
        let values = vec![100, 102, 98, 200, 95];
        let blob = encode_dense_blob(&values, ValueScale::Raw, 5);

        // Check magic
        assert_eq!(&blob[0..4], MAGIC);
        // Check version
        assert_eq!(blob[4], FORMAT_VERSION);
        // Check flags (dense, no metadata)
        assert_eq!(blob[5], 0);
        // Check scale
        assert_eq!(blob[6], 0); // Raw
        // Check num_zoom_levels
        assert_eq!(blob[7], NUM_ZOOM_LEVELS);
        // Check contig_length
        assert_eq!(u32::from_le_bytes([blob[8], blob[9], blob[10], blob[11]]), 5);
    }

    #[test]
    fn test_sparse_blob_roundtrip_header() {
        let positions = vec![10, 50, 100];
        let values = vec![5, 3, 8];
        let flags = MetadataFlags { sparse: true, ..Default::default() };
        let blob = encode_sparse_blob(&positions, &values, None, flags, ValueScale::Times1000, 200);

        // Check magic
        assert_eq!(&blob[0..4], MAGIC);
        // Check flags (sparse=true)
        assert_eq!(blob[5] & 0x01, 1);
        // Check scale
        assert_eq!(blob[6], 2); // Times1000
    }

    #[test]
    fn test_feature_id_mapping() {
        assert!(feature_id("primary_reads").is_some());
        assert!(feature_id("mismatches").is_some());
        assert!(feature_id("nonexistent_feature").is_none());
    }

    #[test]
    fn test_codon_encoding() {
        assert_eq!(codon_to_id("ATG"), 14);
        assert_eq!(codon_to_id("TAA"), 48);
        assert_eq!(codon_to_id("XYZ"), 255);
        assert_eq!(codon_category_to_id("Synonymous"), 0);
        assert_eq!(codon_category_to_id("Non-synonymous"), 1);
        assert_eq!(codon_category_to_id("Intergenic"), 2);
    }

    #[test]
    fn test_encode_dense_blob_small() {
        // Small genome — should produce exactly 1 chunk
        let values: Vec<i32> = (0..1000).collect();
        let blob = encode_dense_blob(&values, ValueScale::Raw, 1000);

        // Verify chunk count
        let base_offset = u32::from_le_bytes([blob[12], blob[13], blob[14], blob[15]]) as usize;
        let num_chunks = u16::from_le_bytes([blob[base_offset], blob[base_offset + 1]]);
        assert_eq!(num_chunks, 1);
    }

    #[test]
    fn test_encode_dense_blob_large() {
        // Large genome — should produce multiple chunks
        let values: Vec<i32> = (0..100_000).collect();
        let blob = encode_dense_blob(&values, ValueScale::Raw, 100_000);

        let base_offset = u32::from_le_bytes([blob[12], blob[13], blob[14], blob[15]]) as usize;
        let num_chunks = u16::from_le_bytes([blob[base_offset], blob[base_offset + 1]]);
        assert_eq!(num_chunks, 2); // 100000 / 65536 = 2 chunks
    }

    #[test]
    fn test_sparse_with_metadata() {
        let positions = vec![10, 20, 30];
        let values = vec![5, 3, 8];
        let metadata = vec![
            EventMeta {
                mean: Some(100), median: Some(95), std: Some(10),
                sequence: Some(b"A".to_vec()), prevalence: Some(850),
                codon_category: Some(1), codon_id: Some(14), aa_id: Some(17),
                ..Default::default()
            },
            EventMeta {
                mean: Some(200), median: Some(190), std: Some(20),
                sequence: Some(b"ACGT".to_vec()), prevalence: Some(600),
                codon_category: Some(0), codon_id: Some(6), aa_id: Some(16),
                ..Default::default()
            },
            EventMeta {
                mean: Some(50), median: Some(48), std: Some(5),
                sequence: None, prevalence: None,
                codon_category: Some(2), codon_id: Some(255), aa_id: Some(255),
                ..Default::default()
            },
        ];
        let flags = MetadataFlags {
            sparse: true, has_stats: true, has_sequence: true, has_codons: true,
        };
        let blob = encode_sparse_blob(&positions, &values, Some(&metadata), flags, ValueScale::Times1000, 100);

        // Just verify it doesn't panic and produces reasonable output
        assert!(blob.len() > 32);
        // Check sparse_meta_offset is non-zero
        let meta_offset = u32::from_le_bytes([blob[24], blob[25], blob[26], blob[27]]);
        assert!(meta_offset > 0);
    }
}
