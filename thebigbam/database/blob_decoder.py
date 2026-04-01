"""
Decode compressed BLOB feature data from Feature_blob table.

Binary format: Header (32 bytes) + BaseResolution + ZoomLevels + SparseMetadata?

This module mirrors the Rust blob encoder (src/blob.rs) and provides
numpy-based decoding for use in Bokeh visualization and CSV export.
"""

import struct
import numpy as np

try:
    import zstandard as zstd
except ImportError:
    zstd = None
    import zlib  # fallback won't work, but gives clear error

# ============================================================================
# Constants
# ============================================================================

MAGIC = b"TBB\x01"
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


def _zigzag_decode(values):
    """Decode zigzag-encoded unsigned integers back to signed."""
    return [(v >> 1) ^ -(v & 1) for v in values]


def _delta_decode(deltas):
    """Reconstruct original values from deltas using prefix sum."""
    if not deltas:
        return []
    values = [deltas[0]]
    for i in range(1, len(deltas)):
        values.append(values[-1] + deltas[i])
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

    all_values = []
    for _ in range(num_chunks):
        compressed_size = struct.unpack_from("<I", blob, offset)[0]
        offset += 4
        compressed_data = blob[offset:offset + compressed_size]
        offset += compressed_size

        # Decompress: zstd → varint → zigzag → delta
        decompressed = _zstd_decompress(compressed_data)
        unsigned_vals = _varint_decode(decompressed)
        signed_deltas = _zigzag_decode(unsigned_vals)
        chunk_values = _delta_decode(signed_deltas)
        all_values.extend(chunk_values)

    # Trim to contig length (last chunk may have padding)
    contig_length = header["contig_length"]
    all_values = all_values[:contig_length]

    return np.array(all_values, dtype=np.int32)


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

    pos_unsigned = _varint_decode(pos_data)
    pos_signed = _zigzag_decode(pos_unsigned)
    positions = _delta_decode(pos_signed)

    # Values: compressed block
    val_compressed_size = struct.unpack_from("<I", blob, offset)[0]
    offset += 4
    val_data = _zstd_decompress(blob[offset:offset + val_compressed_size])
    offset += val_compressed_size

    val_unsigned = _varint_decode(val_data)
    values = _zigzag_decode(val_unsigned)

    return np.array(positions[:event_count], dtype=np.uint32), np.array(values[:event_count], dtype=np.int32)


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
                idx_signed = _zigzag_decode(idx_unsigned)
                bin_indices = _delta_decode(idx_signed)

                # Max values: length-prefixed varint block, zigzag encoded
                val_len = struct.unpack_from("<I", level_data, pos)[0]
                pos += 4
                val_unsigned = _varint_decode(level_data[pos:pos + val_len])
                val_values = _zigzag_decode(val_unsigned)

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
        }
        # Decode metadata if present
        meta = _decode_sparse_metadata(blob_bytes, header, len(positions))
        result.update(meta)
        return result
    else:
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
        }
    else:
        bin_starts = np.arange(n, dtype=np.uint32) * bin_size
        bin_ends = bin_starts + bin_size - 1
        bin_ends[-1] = min(bin_ends[-1], contig_length - 1)
        return {
            "bin_start": bin_starts,
            "bin_end": bin_ends,
            "mean": np.array([b["mean"] for b in bins], dtype=np.float64) / scale,
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
    if header["num_zoom_levels"] == 1:
        available_bins = [10000]
    else:
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


def decode_dense_chunk(blob_bytes, chunk_idx):
    """
    Decode a single 64K chunk for random access on large genomes.

    Args:
        blob_bytes: Raw BLOB bytes
        chunk_idx: 0-based chunk index

    Returns: numpy array of i32 values for this chunk
    """
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)

    header = _parse_header(blob_bytes)
    if header["sparse"]:
        raise ValueError("Cannot decode chunks from sparse BLOB")

    offset = header["base_block_offset"]
    num_chunks = struct.unpack_from("<H", blob_bytes, offset)[0]
    offset += 6  # skip num_chunks(2) + chunk_size(4)

    if chunk_idx >= num_chunks:
        raise ValueError(f"Chunk {chunk_idx} out of range (have {num_chunks} chunks)")

    # Skip to the requested chunk
    for i in range(chunk_idx):
        compressed_size = struct.unpack_from("<I", blob_bytes, offset)[0]
        offset += 4 + compressed_size

    # Decode the target chunk
    compressed_size = struct.unpack_from("<I", blob_bytes, offset)[0]
    offset += 4
    compressed_data = blob_bytes[offset:offset + compressed_size]

    decompressed = _zstd_decompress(compressed_data)
    unsigned_vals = _varint_decode(decompressed)
    signed_deltas = _zigzag_decode(unsigned_vals)
    values = _delta_decode(signed_deltas)

    return np.array(values, dtype=np.int32)


def get_blob_header(blob_bytes):
    """Parse and return the BLOB header as a dict (for inspection/debugging)."""
    if isinstance(blob_bytes, memoryview):
        blob_bytes = bytes(blob_bytes)
    return _parse_header(blob_bytes)


# ============================================================================
# Feature ID Mapping (must match VARIABLES in src/types.rs)
# ============================================================================

# This mapping is derived from the VARIABLES array order in types.rs
# Feature_id = 1-based index into VARIABLES
_FEATURE_NAMES = [
    "direct_repeat_count", "inverted_repeat_count",
    "direct_repeat_identity", "inverted_repeat_identity",
    "gc_content", "gc_skew",
    "primary_reads", "primary_reads_plus_only", "primary_reads_minus_only",
    "secondary_reads", "supplementary_reads", "mapq",
    "left_clippings", "right_clippings", "insertions", "deletions", "mismatches",
    "read_lengths",
    "insert_sizes", "non_inward_pairs", "mate_not_mapped", "mate_on_another_contig",
    "coverage_reduced", "reads_starts", "reads_ends",
]

_NAME_TO_ID = {name: i + 1 for i, name in enumerate(_FEATURE_NAMES)}
_ID_TO_NAME = {i + 1: name for i, name in enumerate(_FEATURE_NAMES)}

# Contig_blob feature IDs (stored in Feature_id column)
_CONTIG_BLOB_IDS = {
    1: "gc_content",
    2: "gc_skew",
    3: "direct_repeat_count",
    4: "inverted_repeat_count",
    5: "direct_repeat_identity",
    6: "inverted_repeat_identity",
}
_CONTIG_BLOB_NAMES = {v: k for k, v in _CONTIG_BLOB_IDS.items()}


def feature_name_to_id(name):
    """Convert feature name to feature_id (1-based) for Feature_blob."""
    return _NAME_TO_ID.get(name)


def feature_id_to_name(feature_id):
    """Convert feature_id (1-based) to feature name for Feature_blob."""
    return _ID_TO_NAME.get(feature_id)


def contig_blob_name_to_id(name):
    """Convert contig feature name to Contig_blob feature_id."""
    return _CONTIG_BLOB_NAMES.get(name)


def contig_blob_id_to_name(feature_id):
    """Convert Contig_blob feature_id to feature name."""
    return _CONTIG_BLOB_IDS.get(feature_id)
