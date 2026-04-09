"""
Decode compressed BLOB feature data from Feature_blob table.

Binary format: Header (32 bytes) + BaseResolution + ZoomLevels + SparseMetadata?

This module mirrors the Rust blob encoder (src/blob.rs) and provides
numpy-based decoding for use in Bokeh visualization and CSV export.
"""

import struct
import numpy as np
from thebigbam_rs import decode_dense_chunk, decode_sparse_chunk

try:
    import zstandard as zstd
except ImportError:
    zstd = None
    import zlib  # fallback won't work, but gives clear error

# ============================================================================
# Constants
# ============================================================================

MAGIC = b"TBB\x01"
ZOOM_MAGIC = b"TBZ\x01"
CHUNK_SIZE = 65536
# Zoom level bin sizes (100bp, 1000bp, 10000bp)
ZOOM_BIN_SIZES = [100, 1000, 10000]

# Scale factor mapping
SCALE_DIVISORS = {0: 1, 1: 100, 2: 1000, 3: 10}

# Codon/AA tables (must match src/blob.rs)
CODON_TABLE = [
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
]

AMINO_ACID_TABLE = [
    "A (Alanine)", "C (Cysteine)", "D (Aspartic acid)", "E (Glutamic acid)",
    "F (Phenylalanine)", "G (Glycine)", "H (Histidine)", "I (Isoleucine)",
    "K (Lysine)", "L (Leucine)", "M (Methionine)", "N (Asparagine)",
    "P (Proline)", "Q (Glutamine)", "R (Arginine)", "S (Serine)",
    "T (Threonine)", "V (Valine)", "W (Tryptophan)", "Y (Tyrosine)",
    "* (Stop)",
]

CODON_CATEGORIES = {0: "Synonymous", 1: "Non-synonymous", 2: "Intergenic"}


# ============================================================================
# Low-level Decoding
# ============================================================================

def _zstd_decompress(data):
    """Decompress zstd data."""
    if zstd is None:
        raise ImportError("zstandard package required for BLOB decoding. Install with: pip install zstandard")
    dctx = zstd.ZstdDecompressor()
    return dctx.decompress(data, max_output_size=100 * 1024 * 1024)  # 100MB max


def _varint_decode(data):
    """Decode LEB128 variable-length integers from bytes."""
    values = []
    i = 0
    while i < len(data):
        val = 0
        shift = 0
        while True:
            if i >= len(data):
                break
            byte = data[i]
            i += 1
            val |= (byte & 0x7F) << shift
            shift += 7
            if byte & 0x80 == 0:
                break
        values.append(val)
    return values



# ============================================================================
# Header Parsing
# ============================================================================

def _parse_header(blob):
    """Parse the 32-byte BLOB header. Returns dict with header fields."""
    if len(blob) < 32:
        raise ValueError(f"BLOB too short ({len(blob)} bytes, need at least 32)")

    magic = blob[0:4]
    if magic != MAGIC:
        raise ValueError(f"Invalid BLOB magic: {magic!r} (expected {MAGIC!r})")

    version = blob[4]
    flags_byte = blob[5]
    scale_code = blob[6]
    num_zoom_levels = blob[7]
    contig_length = struct.unpack_from("<I", blob, 8)[0]
    base_block_offset = struct.unpack_from("<I", blob, 12)[0]
    base_block_compressed_size = struct.unpack_from("<I", blob, 16)[0]
    zoom_index_offset = struct.unpack_from("<I", blob, 20)[0]
    sparse_meta_offset = struct.unpack_from("<I", blob, 24)[0]
    sparse_meta_compressed_size = struct.unpack_from("<I", blob, 28)[0]

    return {
        "version": version,
        "sparse": bool(flags_byte & 0x01),
        "has_stats": bool(flags_byte & 0x02),
        "has_sequence": bool(flags_byte & 0x04),
        "has_codons": bool(flags_byte & 0x08),
        "scale_code": scale_code,
        "scale_divisor": SCALE_DIVISORS.get(scale_code, 1),
        "num_zoom_levels": num_zoom_levels,
        "contig_length": contig_length,
        "base_block_offset": base_block_offset,
        "base_block_compressed_size": base_block_compressed_size,
        "zoom_index_offset": zoom_index_offset,
        "sparse_meta_offset": sparse_meta_offset,
        "sparse_meta_compressed_size": sparse_meta_compressed_size,
    }


# ============================================================================
# Dense BLOB Decoding
# ============================================================================

def _decode_dense_base(blob, header):
    """Decode dense base-resolution data from BLOB."""
    offset = header["base_block_offset"]
    num_chunks = struct.unpack_from("<H", blob, offset)[0]
    chunk_size = struct.unpack_from("<I", blob, offset + 2)[0]
    offset += 6

    all_chunks = []
    for _ in range(num_chunks):
        compressed_size = struct.unpack_from("<I", blob, offset)[0]
        offset += 4
        compressed_data = blob[offset:offset + compressed_size]
        offset += compressed_size

        # Decompress: zstd → varint → vectorized zigzag+delta
        decompressed = _zstd_decompress(compressed_data)
        unsigned_vals = _varint_decode(decompressed)
        arr = np.array(unsigned_vals, dtype=np.int64)
        chunk_values = np.cumsum((arr >> 1) ^ -(arr & 1))
        all_chunks.append(chunk_values)

    # Trim to contig length (last chunk may have padding)
    contig_length = header["contig_length"]
    values = np.concatenate(all_chunks) if all_chunks else np.array([], dtype=np.int64)
    return values[:contig_length].astype(np.int32)


# ============================================================================
# Sparse BLOB Decoding
# ============================================================================

def _decode_sparse_base(blob, header):
    """Decode sparse base-resolution data from BLOB."""
    offset = header["base_block_offset"]
    event_count = struct.unpack_from("<I", blob, offset)[0]
    offset += 4

    if event_count == 0:
        return np.array([], dtype=np.uint32), np.array([], dtype=np.int32)

    # Positions: compressed block
    pos_compressed_size = struct.unpack_from("<I", blob, offset)[0]
    offset += 4
    pos_data = _zstd_decompress(blob[offset:offset + pos_compressed_size])
    offset += pos_compressed_size

    pos_arr = np.array(_varint_decode(pos_data), dtype=np.int64)
    positions = np.cumsum((pos_arr >> 1) ^ -(pos_arr & 1))  # vectorized zigzag+delta

    # Values: compressed block
    val_compressed_size = struct.unpack_from("<I", blob, offset)[0]
    offset += 4
    val_data = _zstd_decompress(blob[offset:offset + val_compressed_size])
    offset += val_compressed_size

    val_arr = np.array(_varint_decode(val_data), dtype=np.int64)
    values = (val_arr >> 1) ^ -(val_arr & 1)  # vectorized zigzag (no delta for values)

    return positions[:event_count].astype(np.uint32), values[:event_count].astype(np.int32)


def decode_raw_chunks(chunk_rows, scale_divisor, chunk_size=CHUNK_SIZE):
    """Decode raw chunk blobs fetched from Feature_blob_chunk / Contig_blob_chunk.

    Args:
        chunk_rows: list of (chunk_idx, raw_bytes) sorted by chunk_idx
        scale_divisor: value scale divisor (1, 100, 1000, or 10)
        chunk_size: positions per chunk (default 65536)

    Returns dict: {"x": ndarray, "y": ndarray, "sparse": False}
    """
    all_chunks = []
    first_chunk_idx = chunk_rows[0][0] if chunk_rows else 0
    for _chunk_idx, raw_bytes in chunk_rows:
        if isinstance(raw_bytes, memoryview):
            raw_bytes = bytes(raw_bytes)
        chunk_values = np.array(decode_dense_chunk(raw_bytes), dtype=np.int64)
        all_chunks.append(chunk_values)

    values = np.concatenate(all_chunks) if all_chunks else np.array([], dtype=np.int64)
    base_offset = first_chunk_idx * chunk_size
    n = len(values)
    x = np.arange(base_offset, base_offset + n, dtype=np.uint32)
    values = values.astype(np.float64)
    if scale_divisor != 1:
        values = values / scale_divisor
    return {"x": x, "y": values, "sparse": False}


def get_scale_from_zoom_blob(zoom_blob_bytes):
    """Extract scale divisor from a standalone zoom BLOB (TBZ format).

    Returns scale divisor (1, 100, 1000, or 10), or 1 if parsing fails.
    """
    if isinstance(zoom_blob_bytes, memoryview):
        zoom_blob_bytes = bytes(zoom_blob_bytes)
    if len(zoom_blob_bytes) < 16 or zoom_blob_bytes[:4] != ZOOM_MAGIC:
        return 1
    scale_code = zoom_blob_bytes[4]
    return SCALE_DIVISORS.get(scale_code, 1)


def is_sparse_zoom_blob(zoom_blob_bytes):
    """Check if a standalone zoom BLOB represents sparse data."""
    if isinstance(zoom_blob_bytes, memoryview):
        zoom_blob_bytes = bytes(zoom_blob_bytes)
    if len(zoom_blob_bytes) < 16 or zoom_blob_bytes[:4] != ZOOM_MAGIC:
        return False
    return bool(zoom_blob_bytes[5] & 0x01)



def _decode_chunk_metadata(meta_data, event_count, has_stats, has_sequence, has_codons):
    """Decode sparse metadata from raw decompressed bytes (same format as full BLOB metadata)."""
    result = {}
    pos = 0
    n = event_count

    if has_stats:
        result["mean"] = np.frombuffer(meta_data[pos:pos + n * 4], dtype="<i4") / 100.0
        pos += n * 4
        result["median"] = np.frombuffer(meta_data[pos:pos + n * 4], dtype="<i4") / 100.0
        pos += n * 4
        result["std"] = np.frombuffer(meta_data[pos:pos + n * 4], dtype="<i4") / 100.0
        pos += n * 4

    if has_sequence:
        sequences = []
        for _ in range(n):
            seq_len = meta_data[pos]
            pos += 1
            if seq_len > 0:
                sequences.append(meta_data[pos:pos + seq_len].decode("ascii", errors="replace"))
                pos += seq_len
            else:
                sequences.append(None)
        result["sequence"] = sequences
        result["sequence_prevalence"] = np.frombuffer(meta_data[pos:pos + n * 2], dtype="<i2").astype(np.float64) / 10.0
        pos += n * 2

    if has_codons:
        cat_bytes = meta_data[pos:pos + n]
        pos += n
        codon_id_bytes = meta_data[pos:pos + n]
        pos += n
        aa_id_bytes = meta_data[pos:pos + n]
        pos += n
        result["codon_category"] = [CODON_CATEGORIES.get(cat_bytes[i]) for i in range(n)]
        result["codon_change"] = [CODON_TABLE[codon_id_bytes[i]] if codon_id_bytes[i] < 64 else None for i in range(n)]
        result["aa_change"] = [AMINO_ACID_TABLE[aa_id_bytes[i]] if aa_id_bytes[i] < 21 else None for i in range(n)]

    return result


def decode_raw_sparse_chunks(chunk_rows, scale_divisor, chunk_size=CHUNK_SIZE):
    """Decode raw sparse chunk blobs from Feature_blob_chunk / Contig_blob_chunk.

    Each chunk blob format:
        [event_count: u32]
        [positions_compressed_size: u32][positions_compressed]
        [values_compressed_size: u32][values_compressed]
        [has_metadata: u8]
        [metadata_compressed_size: u32][metadata_compressed]  // if has_metadata

    Args:
        chunk_rows: list of (chunk_idx, raw_bytes) sorted by chunk_idx
        scale_divisor: value scale divisor (1, 100, 1000, or 10)

    Returns dict: {"x": ndarray, "y": ndarray, "sparse": True, ...metadata}
    """
    all_positions = []
    all_values = []
    all_meta = {}

    for _chunk_idx, raw_bytes in chunk_rows:
        if isinstance(raw_bytes, memoryview):
            raw_bytes = bytes(raw_bytes)
        if not raw_bytes:
            continue

        event_count = struct.unpack_from("<I", raw_bytes, 0)[0]
        if event_count == 0:
            continue

        # Decode positions + values in Rust (zstd + varint + zigzag + delta)
        pos_list, val_list = decode_sparse_chunk(raw_bytes)
        positions = np.array(pos_list, dtype=np.int64)
        values = np.array(val_list, dtype=np.int64)

        all_positions.append(positions[:event_count])
        all_values.append(values[:event_count])

        # Advance offset past positions and values to reach metadata
        off = 4
        pos_size = struct.unpack_from("<I", raw_bytes, off)[0]
        off += 4 + pos_size
        val_size = struct.unpack_from("<I", raw_bytes, off)[0]
        off += 4 + val_size

        # Metadata — flags byte encodes which fields are present (same as BLOB header flags)
        if off < len(raw_bytes):
            flags_byte = raw_bytes[off]
            off += 1
            if flags_byte:
                has_stats = bool(flags_byte & 0x02)
                has_seq = bool(flags_byte & 0x04)
                has_codons = bool(flags_byte & 0x08)
                meta_size = struct.unpack_from("<I", raw_bytes, off)[0]
                off += 4
                meta_data = _zstd_decompress(raw_bytes[off:off + meta_size])
                chunk_meta = _decode_chunk_metadata(meta_data, event_count, has_stats, has_seq, has_codons)
                for key, val in chunk_meta.items():
                    if key not in all_meta:
                        all_meta[key] = []
                    if isinstance(val, np.ndarray):
                        all_meta[key].append(val)
                    elif isinstance(val, list):
                        all_meta[key].extend(val)

    if not all_positions:
        return {"x": np.array([], dtype=np.uint32), "y": np.array([], dtype=np.float64), "sparse": True}

    x = np.concatenate(all_positions).astype(np.uint32)
    y = np.concatenate(all_values).astype(np.float64)
    if scale_divisor != 1:
        y = y / scale_divisor

    result = {"x": x, "y": y, "sparse": True}
    for key, val in all_meta.items():
        if isinstance(val, list) and val and isinstance(val[0], np.ndarray):
            result[key] = np.concatenate(val)
        else:
            result[key] = val
    return result


# ============================================================================
# Sparse Metadata Decoding
# ============================================================================

def _decode_sparse_metadata(blob, header, event_count):
    """Decode sparse metadata block."""
    meta_offset = header["sparse_meta_offset"]
    if meta_offset == 0 or event_count == 0:
        return {}

    meta_compressed_size = struct.unpack_from("<I", blob, meta_offset)[0]
    meta_offset += 4
    meta_data = _zstd_decompress(blob[meta_offset:meta_offset + meta_compressed_size])

    result = {}
    pos = 0
    n = event_count

    # Stats: mean[], median[], std[] as i32 arrays
    if header["has_stats"]:
        mean_bytes = meta_data[pos:pos + n * 4]
        result["mean"] = np.frombuffer(mean_bytes, dtype="<i4") / 100.0
        pos += n * 4

        median_bytes = meta_data[pos:pos + n * 4]
        result["median"] = np.frombuffer(median_bytes, dtype="<i4") / 100.0
        pos += n * 4

        std_bytes = meta_data[pos:pos + n * 4]
        result["std"] = np.frombuffer(std_bytes, dtype="<i4") / 100.0
        pos += n * 4

    # Sequences: [len: u8, bytes[len]] per event, then prevalence[n] as i16
    if header["has_sequence"]:
        sequences = []
        for _ in range(n):
            seq_len = meta_data[pos]
            pos += 1
            if seq_len > 0:
                seq = meta_data[pos:pos + seq_len].decode("ascii", errors="replace")
                pos += seq_len
                sequences.append(seq)
            else:
                sequences.append(None)
        result["sequence"] = sequences

        prev_bytes = meta_data[pos:pos + n * 2]
        result["sequence_prevalence"] = np.frombuffer(prev_bytes, dtype="<i2").astype(np.float64) / 10.0
        pos += n * 2

    # Codons: category[], codon_id[], aa_id[] each as u8 arrays
    if header["has_codons"]:
        cat_bytes = meta_data[pos:pos + n]
        pos += n
        codon_id_bytes = meta_data[pos:pos + n]
        pos += n
        aa_id_bytes = meta_data[pos:pos + n]
        pos += n

        categories = []
        codons = []
        amino_acids = []
        for i in range(n):
            cat = cat_bytes[i]
            categories.append(CODON_CATEGORIES.get(cat))
            cid = codon_id_bytes[i]
            codons.append(CODON_TABLE[cid] if cid < 64 else None)
            aid = aa_id_bytes[i]
            amino_acids.append(AMINO_ACID_TABLE[aid] if aid < 21 else None)

        result["codon_category"] = categories
        result["codon_change"] = codons
        result["aa_change"] = amino_acids

    return result


# ============================================================================
# Zoom Level Decoding
# ============================================================================

def _decode_zoom_levels(blob, header):
    """Decode all zoom levels from BLOB.

    Format: sparse features store only max_value(i32), dense features store only mean(i32).
    """
    offset = header["zoom_index_offset"]
    is_sparse = header["sparse"]
    levels = []

    for level_idx in range(header["num_zoom_levels"]):
        compressed_size = struct.unpack_from("<I", blob, offset)[0]
        offset += 4
        level_data = _zstd_decompress(blob[offset:offset + compressed_size])
        offset += compressed_size

        if header["num_zoom_levels"] == 1:
            zoom_bin_sizes = [10000]
        else:
            zoom_bin_sizes = ZOOM_BIN_SIZES
        bin_size = zoom_bin_sizes[level_idx] if level_idx < len(zoom_bin_sizes) else 10000
        contig_length = header["contig_length"]
        num_bins = (contig_length + bin_size - 1) // bin_size

        # Simplified format: sparse stores only nonzero bins, dense stores mean per bin
        if is_sparse:
            # Sparse: nonzero_count + delta-encoded bin_indices + zigzag-encoded max_values
            pos = 0
            nonzero_count = struct.unpack_from("<I", level_data, pos)[0]
            pos += 4
            bins = []
            if nonzero_count > 0:
                # Bin indices: length-prefixed varint block, delta+zigzag encoded
                idx_len = struct.unpack_from("<I", level_data, pos)[0]
                pos += 4
                idx_unsigned = _varint_decode(level_data[pos:pos + idx_len])
                pos += idx_len
                idx_arr = np.array(idx_unsigned, dtype=np.int64)
                bin_indices = np.cumsum((idx_arr >> 1) ^ -(idx_arr & 1))

                # Max values: length-prefixed varint block, zigzag encoded
                val_len = struct.unpack_from("<I", level_data, pos)[0]
                pos += 4
                val_unsigned = _varint_decode(level_data[pos:pos + val_len])
                val_arr = np.array(val_unsigned, dtype=np.int64)
                val_values = (val_arr >> 1) ^ -(val_arr & 1)

                for i in range(nonzero_count):
                    bins.append({"idx": bin_indices[i], "max": val_values[i], "mean": val_values[i]})
        else:
            # Dense: only mean(i32) = 4 bytes per bin
            bins = []
            pos = 0
            for _ in range(num_bins):
                if pos + 4 > len(level_data):
                    break
                mean_val = struct.unpack_from("<i", level_data, pos)[0]
                pos += 4
                bins.append({"mean": mean_val})

        levels.append({"bin_size": bin_size, "bins": bins})

    return levels


# ============================================================================
# Public API
# ============================================================================

def decode_blob(blob_bytes):
    """
    Decode a Feature_blob BLOB to arrays.

    Returns dict:
    - For dense: {"x": ndarray, "y": ndarray}
    - For sparse: {"x": ndarray, "y": ndarray, ...metadata fields}

    Values are automatically descaled (e.g., ÷100 for Times100, ÷1000 for Times1000).
    Positions in x are 0-indexed.
    """
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)

    header = _parse_header(blob_bytes)
    scale = header["scale_divisor"]

    if header["sparse"]:
        positions, values = _decode_sparse_base(blob_bytes, header)
        result = {
            "x": positions,
            "y": values.astype(np.float64) / scale if scale != 1 else values.astype(np.float64),
            "sparse": True,
        }
        # Decode metadata if present
        meta = _decode_sparse_metadata(blob_bytes, header, len(positions))
        result.update(meta)
        return result
    else:
        # decode_blob always decodes all chunks (no range)
        values = _decode_dense_base(blob_bytes, header)
        n = len(values)
        contig_length = header["contig_length"]
        if n > 0 and n < contig_length:
            # Windowed data (e.g. GC content 500bp, GC skew 1000bp): convert indices to window midpoints
            window_size = contig_length // n
            x = np.arange(n, dtype=np.uint32) * window_size + window_size // 2
        else:
            # Base-pair resolution: index = position
            x = np.arange(n, dtype=np.uint32)
        return {
            "x": x,
            "y": values.astype(np.float64) / scale if scale != 1 else values.astype(np.float64),
            "sparse": False,
        }


def decode_zoom_level(blob_bytes, level):
    """
    Decode pre-computed zoom summaries for a specific level index.

    Note: Prefer decode_zoom_by_bin_size() for more flexible access.

    Args:
        blob_bytes: Raw BLOB bytes
        level: Zoom level index (0=100bp, 1=1000bp, 2=10000bp)

    Returns dict:
        {"bin_start": ndarray, "bin_end": ndarray, "mean": ndarray, ...}
        For sparse features: includes "max", only nonzero bins
        For dense features: includes "mean", all bins
    """
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)

    header = _parse_header(blob_bytes)
    scale = header["scale_divisor"]
    levels = _decode_zoom_levels(blob_bytes, header)

    if level >= len(levels):
        raise ValueError(f"Zoom level {level} not available (have {len(levels)} levels)")

    zoom = levels[level]
    return _zoom_bins_to_dict(
        zoom["bins"], zoom["bin_size"], header["contig_length"], scale, header["sparse"]
    )


def _zoom_bins_to_dict(bins, bin_size, contig_length, scale, is_sparse):
    """Convert decoded zoom bins to a result dict with bin_start/bin_end/values."""
    n = len(bins)
    if n == 0:
        return None

    if is_sparse:
        bin_indices = np.array([b["idx"] for b in bins], dtype=np.uint32)
        bin_starts = bin_indices * bin_size
        bin_ends = bin_starts + bin_size - 1
        bin_ends[-1] = min(bin_ends[-1], contig_length - 1)
        return {
            "bin_start": bin_starts,
            "bin_end": bin_ends,
            "max": np.array([b["max"] for b in bins], dtype=np.float64) / scale,
            "mean": np.array([b["mean"] for b in bins], dtype=np.float64) / scale,
            "sparse": True,
        }
    else:
        bin_starts = np.arange(n, dtype=np.uint32) * bin_size
        bin_ends = bin_starts + bin_size - 1
        bin_ends[-1] = min(bin_ends[-1], contig_length - 1)
        return {
            "bin_start": bin_starts,
            "bin_end": bin_ends,
            "mean": np.array([b["mean"] for b in bins], dtype=np.float64) / scale,
            "sparse": False,
        }


def decode_zoom_by_bin_size(blob_bytes, target_bin_size):
    """
    Decode pre-computed zoom summaries by target bin size with fallback.

    Available bin sizes: 100bp, 1000bp, 10000bp

    If target_bin_size is not available, tries to use the next larger available bin size.
    For Contig_blob (1 zoom level), uses base resolution if 10kbp is not available.

    Args:
        blob_bytes: Raw BLOB bytes
        target_bin_size: Desired bin size (100, 1000, 10000, etc.)

    Returns dict or None if no suitable bin size is available.
    """
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)

    header = _parse_header(blob_bytes)
    scale = header["scale_divisor"]
    is_sparse = header["sparse"]
    contig_length = header["contig_length"]
    levels = _decode_zoom_levels(blob_bytes, header)

    # Try exact match first, then next larger bin size
    available_bins = [100, 1000, 10000]
    bin_sizes_to_try = [target_bin_size]
    fallback = [b for b in available_bins if b > target_bin_size]
    if fallback:
        bin_sizes_to_try.append(min(fallback))

    for try_size in bin_sizes_to_try:
        for zoom in levels:
            if zoom["bin_size"] == try_size:
                result = _zoom_bins_to_dict(zoom["bins"], try_size, contig_length, scale, is_sparse)
                return result

    return None


def decode_zoom_standalone(zoom_blob_bytes, target_bin_size):
    """
    Decode zoom levels from a standalone Zoom_data BLOB (TBZ format).

    This is the lightweight alternative to decode_zoom_by_bin_size() that works
    on the separate Zoom_data column, avoiding the need to fetch the full Data BLOB.

    Format: ZOOM_MAGIC (4) + scale_code (1) + flags (1) + num_zoom_levels (1)
          + reserved (1) + contig_length (4) + reserved (4) = 16-byte header + zoom data

    Returns dict or None if no suitable bin size is available.
    """
    if isinstance(zoom_blob_bytes, memoryview):
        zoom_blob_bytes = bytes(zoom_blob_bytes)

    if len(zoom_blob_bytes) < 16 or zoom_blob_bytes[:4] != ZOOM_MAGIC:
        return None

    scale_code = zoom_blob_bytes[4]
    flags_byte = zoom_blob_bytes[5]
    num_zoom_levels = zoom_blob_bytes[6]
    contig_length = struct.unpack_from("<I", zoom_blob_bytes, 8)[0]

    scale = SCALE_DIVISORS.get(scale_code, 1)
    is_sparse = (flags_byte & 0x01) != 0

    if num_zoom_levels == 1:
        zoom_bin_sizes = [10000]
    else:
        zoom_bin_sizes = ZOOM_BIN_SIZES

    # Decode zoom levels from offset 16
    offset = 16
    levels = []
    for level_idx in range(num_zoom_levels):
        if offset + 4 > len(zoom_blob_bytes):
            break
        compressed_size = struct.unpack_from("<I", zoom_blob_bytes, offset)[0]
        offset += 4
        if offset + compressed_size > len(zoom_blob_bytes):
            break
        level_data = _zstd_decompress(zoom_blob_bytes[offset:offset + compressed_size])
        offset += compressed_size

        bin_size = zoom_bin_sizes[level_idx] if level_idx < len(zoom_bin_sizes) else 10000
        num_bins = (contig_length + bin_size - 1) // bin_size

        if is_sparse:
            pos = 0
            nonzero_count = struct.unpack_from("<I", level_data, pos)[0]
            pos += 4
            bins = []
            if nonzero_count > 0:
                idx_len = struct.unpack_from("<I", level_data, pos)[0]
                pos += 4
                idx_unsigned = _varint_decode(level_data[pos:pos + idx_len])
                pos += idx_len
                idx_arr = np.array(idx_unsigned, dtype=np.int64)
                bin_indices = np.cumsum((idx_arr >> 1) ^ -(idx_arr & 1))

                val_len = struct.unpack_from("<I", level_data, pos)[0]
                pos += 4
                val_unsigned = _varint_decode(level_data[pos:pos + val_len])
                val_arr = np.array(val_unsigned, dtype=np.int64)
                val_values = (val_arr >> 1) ^ -(val_arr & 1)

                for i in range(nonzero_count):
                    bins.append({"idx": bin_indices[i], "max": val_values[i], "mean": val_values[i]})
        else:
            bins = []
            pos = 0
            for _ in range(num_bins):
                if pos + 4 > len(level_data):
                    break
                mean_val = struct.unpack_from("<i", level_data, pos)[0]
                pos += 4
                bins.append({"mean": mean_val})

        levels.append({"bin_size": bin_size, "bins": bins})

    # Try exact match first, then next larger bin size
    available_bins = [100, 1000, 10000]
    bin_sizes_to_try = [target_bin_size]
    fallback = [b for b in available_bins if b > target_bin_size]
    if fallback:
        bin_sizes_to_try.append(min(fallback))

    for try_size in bin_sizes_to_try:
        for zoom in levels:
            if zoom["bin_size"] == try_size:
                return _zoom_bins_to_dict(zoom["bins"], try_size, contig_length, scale, is_sparse)

    return None


def get_blob_header(blob_bytes):
    """Parse and return the BLOB header as a dict (for inspection/debugging)."""
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)
    return _parse_header(blob_bytes)


# ============================================================================
# Feature ID Mapping (read from DB Variable table)
# ============================================================================

# Contig_blob feature IDs (stored in Feature_id column)
# These are a separate ID space with fixed values assigned in Rust.
_CONTIG_BLOB_IDS = {
    1: "gc_content",
    2: "gc_skew",
    3: "direct_repeat_count",
    4: "inverted_repeat_count",
    5: "direct_repeat_identity",
    6: "inverted_repeat_identity",
}
_CONTIG_BLOB_NAMES = {v: k for k, v in _CONTIG_BLOB_IDS.items()}


def feature_name_to_id(name, conn):
    """Look up Feature_id from the Variable table (1-based).

    Uses the DB as source of truth so Python never goes out of sync
    with the Rust VARIABLES array.
    """
    row = conn.execute(
        "SELECT Variable_id FROM Variable WHERE Variable_name = ?", [name]
    ).fetchone()
    return row[0] if row else None


def contig_blob_name_to_id(name):
    """Convert contig feature name to Contig_blob feature_id."""
    return _CONTIG_BLOB_NAMES.get(name)


