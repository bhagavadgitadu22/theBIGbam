"""
CLI command: thebigbam inspect

Decode and display Feature_blob or Contig_blob contents as TSV.
Supports comma-separated lists for --feature, --contig, and --sample.
"""

import sys

import numpy as np
import duckdb

from thebigbam.database.blob_decoder import (
    decode_zoom_standalone,
    decode_raw_chunks, decode_raw_sparse_chunks,
    get_scale_from_zoom_blob, is_sparse_zoom_blob,
    feature_name_to_id, contig_blob_name_to_id,
    CHUNK_SIZE, SCALE_DIVISORS, ZOOM_MAGIC, CODON_CATEGORIES,
)


def add_inspect_args(parser):
    parser.add_argument('-d', '--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--feature', required=True, help='Feature name(s), comma-separated (e.g. mismatches,coverage,gc_content)')
    parser.add_argument('--contig', required=True, help='Contig name(s), comma-separated')
    parser.add_argument('--sample', default=None, help='Sample name(s), comma-separated (required for sample-level features, omit for contig-level)')
    parser.add_argument('--region', default=None, help='Genomic region to display (e.g. 1000-2000)')
    parser.add_argument('--zoom', type=int, default=None, choices=[0, 1, 2],
                        help='Show zoom level instead of base resolution (0=100bp, 1=1000bp, 2=10000bp)')
    parser.add_argument('--header-only', action='store_true', help='Only show BLOB header info')


def run_inspect(args):
    if args.header_only and args.zoom is not None:
        sys.exit("ERROR: --header-only and --zoom are mutually exclusive.")

    conn = duckdb.connect(args.db, read_only=True)

    features = [f.strip() for f in args.feature.split(',')]
    contigs = [c.strip() for c in args.contig.split(',')]
    samples = [s.strip() for s in args.sample.split(',')] if args.sample else []

    # Parse region filter
    region_start, region_end = None, None
    if args.region:
        parts = args.region.split('-')
        if len(parts) != 2:
            sys.exit("ERROR: --region must be in format START-END (e.g. 1000-2000)")
        region_start, region_end = int(parts[0]), int(parts[1])

    zoom_bin_size = {0: 100, 1: 1000, 2: 10000}.get(args.zoom)

    for contig_name in contigs:
        # Resolve contig_id
        row = conn.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name]).fetchone()
        if row is None:
            print(f"# ERROR: Contig '{contig_name}' not found, skipping.", file=sys.stderr)
            continue
        contig_id = row[0]

        for feature_name in features:
            contig_fid = contig_blob_name_to_id(feature_name)
            sample_fid = feature_name_to_id(feature_name, conn)

            if contig_fid is not None:
                # Contig-level feature
                _inspect_contig_feature(
                    conn, contig_id, contig_name, feature_name, contig_fid,
                    region_start, region_end, zoom_bin_size, args.header_only,
                )
            elif sample_fid is not None:
                # Sample-level feature
                if not samples:
                    print(f"# ERROR: '{feature_name}' is sample-level — --sample required, skipping.", file=sys.stderr)
                    continue
                for sample_name in samples:
                    row = conn.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", [sample_name]).fetchone()
                    if row is None:
                        print(f"# ERROR: Sample '{sample_name}' not found, skipping.", file=sys.stderr)
                        continue
                    sample_id = row[0]
                    _inspect_sample_feature(
                        conn, contig_id, contig_name, sample_name, sample_id,
                        feature_name, sample_fid,
                        region_start, region_end, zoom_bin_size, args.header_only,
                    )
            else:
                print(f"# ERROR: Unknown feature '{feature_name}', skipping.", file=sys.stderr)

    conn.close()
    return 0


def _inspect_contig_feature(conn, contig_id, contig_name, feature_name, fid,
                            region_start, region_end, zoom_bin_size, header_only):
    """Inspect a contig-level feature from Contig_blob / Contig_blob_chunk."""
    # Fetch zoom blob
    row = conn.execute(
        "SELECT Zoom_data FROM Contig_blob WHERE Contig_id = ? AND Feature_id = ?",
        [contig_id, fid],
    ).fetchone()
    if row is None:
        print(f"# No data for contig='{contig_name}', feature='{feature_name}'", file=sys.stderr)
        return
    zoom_blob = bytes(row[0])

    print(f"# contig={contig_name} feature={feature_name}")

    if header_only:
        _print_zoom_header(zoom_blob, contig_name, feature_name)
        return

    if zoom_bin_size is not None:
        _print_zoom_data(zoom_blob, zoom_bin_size, region_start, region_end)
        return

    # Base resolution from chunks
    scale_div = get_scale_from_zoom_blob(zoom_blob)
    chunk_rows = _fetch_contig_chunks(conn, contig_id, fid, region_start, region_end, feature_name)
    if not chunk_rows:
        print("# No chunk data found", file=sys.stderr)
        return

    data = decode_raw_chunks(chunk_rows, scale_div)
    x_arr = data["x"]
    y_arr = data["y"]

    # GC content/skew: convert window indices to genomic positions
    gc_window = 500 if feature_name == "gc_content" else (1000 if feature_name == "gc_skew" else 1)
    if gc_window > 1:
        x_arr = x_arr * gc_window

    _apply_region_and_print(x_arr, y_arr, data, region_start, region_end)


def _inspect_sample_feature(conn, contig_id, contig_name, sample_name, sample_id,
                            feature_name, fid,
                            region_start, region_end, zoom_bin_size, header_only):
    """Inspect a sample-level feature from Feature_blob / Feature_blob_chunk."""
    # Fetch zoom blob
    row = conn.execute(
        "SELECT Zoom_data FROM Feature_blob WHERE Contig_id = ? AND Sample_id = ? AND Feature_id = ?",
        [contig_id, sample_id, fid],
    ).fetchone()
    if row is None:
        print(f"# No data for contig='{contig_name}', sample='{sample_name}', feature='{feature_name}'", file=sys.stderr)
        return
    zoom_blob = bytes(row[0])

    print(f"# contig={contig_name} sample={sample_name} feature={feature_name}")

    if header_only:
        _print_zoom_header(zoom_blob, contig_name, feature_name, sample_name)
        return

    if zoom_bin_size is not None:
        _print_zoom_data(zoom_blob, zoom_bin_size, region_start, region_end)
        return

    # Base resolution from chunks
    scale_div = get_scale_from_zoom_blob(zoom_blob)
    sparse = is_sparse_zoom_blob(zoom_blob)
    chunk_rows = _fetch_feature_chunks(conn, contig_id, sample_id, fid, region_start, region_end)
    if not chunk_rows:
        print("# No chunk data found", file=sys.stderr)
        return

    if sparse:
        data = decode_raw_sparse_chunks(chunk_rows, scale_div)
    else:
        data = decode_raw_chunks(chunk_rows, scale_div)

    x_arr = data["x"]
    y_arr = data["y"]
    _apply_region_and_print(x_arr, y_arr, data, region_start, region_end)


def _fetch_contig_chunks(conn, contig_id, fid, region_start, region_end, feature_name):
    """Fetch Contig_blob_chunk rows, optionally limited to region."""
    gc_window = 500 if feature_name == "gc_content" else (1000 if feature_name == "gc_skew" else 1)
    if region_start is not None:
        # Convert genomic positions to chunk indices
        first_idx = max(0, (region_start - 1) // gc_window // CHUNK_SIZE)
        last_idx = (region_end - 1) // gc_window // CHUNK_SIZE
        rows = conn.execute(
            "SELECT Chunk_idx, Data FROM Contig_blob_chunk "
            "WHERE Contig_id=? AND Feature_id=? AND Chunk_idx BETWEEN ? AND ? ORDER BY Chunk_idx",
            [contig_id, fid, first_idx, last_idx],
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT Chunk_idx, Data FROM Contig_blob_chunk "
            "WHERE Contig_id=? AND Feature_id=? ORDER BY Chunk_idx",
            [contig_id, fid],
        ).fetchall()
    return [(r[0], r[1]) for r in rows] if rows else []


def _fetch_feature_chunks(conn, contig_id, sample_id, fid, region_start, region_end):
    """Fetch Feature_blob_chunk rows, optionally limited to region."""
    if region_start is not None:
        first_idx = max(0, (region_start - 1) // CHUNK_SIZE)
        last_idx = (region_end - 1) // CHUNK_SIZE
        rows = conn.execute(
            "SELECT Chunk_idx, Data FROM Feature_blob_chunk "
            "WHERE Contig_id=? AND Sample_id=? AND Feature_id=? AND Chunk_idx BETWEEN ? AND ? ORDER BY Chunk_idx",
            [contig_id, sample_id, fid, first_idx, last_idx],
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT Chunk_idx, Data FROM Feature_blob_chunk "
            "WHERE Contig_id=? AND Sample_id=? AND Feature_id=? ORDER BY Chunk_idx",
            [contig_id, sample_id, fid],
        ).fetchall()
    return [(r[0], r[1]) for r in rows] if rows else []


def _print_zoom_header(zoom_blob, contig_name, feature_name, sample_name=None):
    """Print zoom blob header info."""
    if len(zoom_blob) < 16 or zoom_blob[:4] != ZOOM_MAGIC:
        print("  Invalid zoom blob", file=sys.stderr)
        return
    scale_code = zoom_blob[4]
    flags = zoom_blob[5]
    num_levels = zoom_blob[6]
    import struct
    contig_length = struct.unpack_from("<I", zoom_blob, 8)[0]
    sparse = bool(flags & 0x01)

    parts = [f"contig={contig_name}", f"feature={feature_name}"]
    if sample_name:
        parts.append(f"sample={sample_name}")
    print(f"  {'  '.join(parts)}")
    print(f"  Zoom blob: {len(zoom_blob):,} bytes")
    print(f"  Contig length: {contig_length:,}")
    print(f"  Sparse: {sparse}, Scale: ÷{SCALE_DIVISORS.get(scale_code, 1)}")
    print(f"  Zoom levels: {num_levels}")


def _print_zoom_data(zoom_blob, bin_size, region_start, region_end):
    """Print zoom summary bins as TSV."""
    data = decode_zoom_standalone(zoom_blob, bin_size)
    if data is None:
        print(f"# No zoom data for bin_size={bin_size}", file=sys.stderr)
        return

    bin_starts = data["bin_start"]
    bin_ends = data["bin_end"]

    if region_start is not None:
        mask = (bin_starts >= (region_start - 1)) & (bin_starts <= (region_end - 1))
        bin_starts = bin_starts[mask]
        bin_ends = bin_ends[mask]
        for k in data:
            if k not in ("bin_start", "bin_end") and hasattr(data[k], '__len__') and len(data[k]) == len(mask):
                data[k] = data[k][mask]

    if len(bin_starts) == 0:
        return

    has_max = "max" in data
    if has_max:
        print("bin_start\tbin_end\tmax\tmean")
        for i in range(len(bin_starts)):
            print(f"{int(bin_starts[i]) + 1}\t{int(bin_ends[i]) + 1}\t{data['max'][i]:.4f}\t{data['mean'][i]:.4f}")
    else:
        print("bin_start\tbin_end\tmean")
        for i in range(len(bin_starts)):
            print(f"{int(bin_starts[i]) + 1}\t{int(bin_ends[i]) + 1}\t{data['mean'][i]:.4f}")


def _apply_region_and_print(x_arr, y_arr, data, region_start, region_end):
    """Apply region filter and print base-resolution TSV."""
    if region_start is not None:
        mask = (x_arr >= (region_start - 1)) & (x_arr <= (region_end - 1))
        orig_len = len(x_arr)
        x_arr = x_arr[mask]
        y_arr = y_arr[mask]
        for k in data:
            if k in ("x", "y", "sparse"):
                continue
            v = data[k]
            if hasattr(v, '__len__') and len(v) == orig_len:
                data[k] = v[mask] if isinstance(v, np.ndarray) else [v[i] for i, m in enumerate(mask) if m]

    n = len(x_arr)
    if n == 0:
        return

    has_stats = "mean" in data
    has_seq = "sequence" in data
    has_prevalence = "sequence_prevalence" in data
    has_codons = "codon_category" in data

    cols = ["position", "value"]
    if has_stats:
        cols.extend(["mean", "median", "std"])
    if has_seq:
        cols.append("sequence")
    if has_prevalence:
        cols.append("sequence_prevalence")
    if has_codons:
        cols.extend(["codon_category", "codon", "amino_acid"])

    print("\t".join(cols))

    for i in range(n):
        parts = [str(int(x_arr[i]) + 1), f"{y_arr[i]:.4f}"]
        if has_stats:
            parts.extend([
                f"{data['mean'][i]:.2f}",
                f"{data['median'][i]:.2f}",
                f"{data['std'][i]:.2f}",
            ])
        if has_seq:
            seq = data["sequence"][i]
            parts.append(seq if seq else "")
        if has_prevalence:
            parts.append(f"{data['sequence_prevalence'][i]:.1f}")
        if has_codons:
            cat = data["codon_category"][i]
            cat_str = CODON_CATEGORIES.get(cat, str(cat)) if isinstance(cat, (int, float)) else str(cat)
            parts.append(cat_str)
            parts.append(data["codon_change"][i] if "codon_change" in data and data["codon_change"][i] else "")
            parts.append(data["aa_change"][i] if "aa_change" in data and data["aa_change"][i] else "")
        print("\t".join(parts))
