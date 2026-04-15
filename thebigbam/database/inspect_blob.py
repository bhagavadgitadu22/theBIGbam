"""
CLI command: thebigbam inspect

Decode and display Feature_blob or Contig_blob contents as TSV.
Supports comma-separated lists for --feature, --contig, and --sample.

Output format: contig\tsample\tfeature\tposition_start\tposition_end\tvalue
Consecutive positions with the same value are collapsed into one row (RLE).
"""

import sys
import io as _io

import numpy as np
import duckdb

from thebigbam.database.blob_decoder import (
    decode_zoom_standalone,
    decode_raw_chunks, decode_raw_sparse_chunks,
    get_scale_from_zoom_blob, is_sparse_zoom_blob,
    feature_name_to_id, is_contig_blob_feature,
    CHUNK_SIZE, SCALE_DIVISORS, ZOOM_MAGIC,
)


def add_inspect_args(parser):
    parser.add_argument('-d', '--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--feature', required=True, help='Feature name(s), comma-separated (e.g. mismatches,coverage,gc_content)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--contig', help='Contig name(s), comma-separated')
    group.add_argument('--mag', help='MAG name(s), comma-separated (expands to all member contigs; adds MAG column to output)')
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
    samples = [s.strip() for s in args.sample.split(',')] if args.sample else []

    # Parse region filter
    region_start, region_end = None, None
    if args.region:
        parts = args.region.split('-')
        if len(parts) != 2:
            sys.exit("ERROR: --region must be in format START-END (e.g. 1000-2000)")
        region_start, region_end = int(parts[0]), int(parts[1])

    zoom_bin_size = {0: 100, 1: 1000, 2: 10000}.get(args.zoom)

    # Resolve contigs and MAG association
    mag_of_contig = {}   # contig_name → mag_name (empty when using --contig)
    contig_offsets = {}  # contig_name → Offset_in_MAG (0 when using --contig)
    mags = [m.strip() for m in args.mag.split(',')] if args.mag else []
    if args.mag:
        # Expand MAG names → member contigs, recording the MAG association and offset
        contigs = []
        for mag_name in mags:
            rows = conn.execute(
                "SELECT c.Contig_name, mca.Offset_in_MAG FROM Contig c "
                "JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id "
                "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
                "WHERE mg.MAG_name = ?",
                [mag_name],
            ).fetchall()
            if not rows:
                print(f"# ERROR: MAG '{mag_name}' not found or has no member contigs, skipping.", file=sys.stderr)
                continue
            for cn, offset in rows:
                contigs.append(cn)
                mag_of_contig[cn] = mag_name
                contig_offsets[cn] = int(offset)
    else:
        contigs = [c.strip() for c in args.contig.split(',')]

    # Print header once (use buffer.write to avoid mixed text/binary buffering)
    if not args.header_only:
        if mag_of_contig:
            sys.stdout.buffer.write(b"mag\tcontig\tsample\tfeature\tposition_start\tposition_end\tvalue\n")
        else:
            sys.stdout.buffer.write(b"contig\tsample\tfeature\tposition_start\tposition_end\tvalue\n")

    # When using --mag with --zoom, query MAG-level zoom blobs directly.
    # They are pre-computed in MAG-global coordinates — no per-contig offset needed.
    if args.mag and zoom_bin_size is not None:
        from thebigbam.database.database_getters import (
            get_mag_id, get_mag_feature_zoom, get_mag_contig_zoom,
        )
        for mag_name in mags:
            mag_id = get_mag_id(conn, mag_name)
            if mag_id is None:
                print(f"# ERROR: MAG '{mag_name}' not found, skipping.", file=sys.stderr)
                continue
            for feature_name in features:
                if is_contig_blob_feature(feature_name):
                    zoom_blob = get_mag_contig_zoom(conn, mag_id, feature_name)
                    if zoom_blob is None:
                        continue
                    _print_zoom_rows(zoom_blob, zoom_bin_size, "-", "", feature_name,
                                     region_start, region_end, mag_name=mag_name)
                else:
                    if not samples:
                        print(f"# ERROR: '{feature_name}' is sample-level — --sample required, skipping.", file=sys.stderr)
                        continue
                    for sample_name in samples:
                        s_row = conn.execute(
                            "SELECT Sample_id FROM Sample WHERE Sample_name = ?", [sample_name]
                        ).fetchone()
                        if s_row is None:
                            print(f"# ERROR: Sample '{sample_name}' not found, skipping.", file=sys.stderr)
                            continue
                        zoom_blob = get_mag_feature_zoom(conn, mag_id, s_row[0], feature_name)
                        if zoom_blob is None:
                            continue
                        _print_zoom_rows(zoom_blob, zoom_bin_size, "-", sample_name, feature_name,
                                         region_start, region_end, mag_name=mag_name)
        conn.close()
        return 0

    for contig_name in contigs:
        row = conn.execute(
            "SELECT Contig_id, Contig_length FROM Contig WHERE Contig_name = ?", [contig_name]
        ).fetchone()
        if row is None:
            print(f"# ERROR: Contig '{contig_name}' not found, skipping.", file=sys.stderr)
            continue
        contig_id, contig_length = int(row[0]), int(row[1])
        contig_mag = mag_of_contig.get(contig_name, "")
        mag_offset = contig_offsets.get(contig_name, 0)

        # Translate MAG-global region to contig-local coordinates
        if region_start is not None and mag_offset > 0:
            local_rs = region_start - mag_offset
            local_re = region_end - mag_offset
            # Skip this contig entirely if the region doesn't overlap it
            if local_re < 1 or local_rs > contig_length:
                continue
            local_rs = max(1, local_rs)
            local_re = min(contig_length, local_re)
        else:
            local_rs = region_start
            local_re = region_end

        for feature_name in features:
            fid = feature_name_to_id(feature_name, conn)

            if is_contig_blob_feature(feature_name):
                if fid is None:
                    print(f"# ERROR: Unknown feature '{feature_name}', skipping.", file=sys.stderr)
                    continue
                _inspect_contig_feature(
                    conn, contig_id, contig_name, feature_name, fid,
                    local_rs, local_re, zoom_bin_size, args.header_only,
                    mag_name=contig_mag, mag_offset=mag_offset,
                )
            elif fid is not None:
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
                        feature_name, fid,
                        local_rs, local_re, zoom_bin_size, args.header_only,
                        mag_name=contig_mag, mag_offset=mag_offset,
                    )
            else:
                print(f"# ERROR: Unknown feature '{feature_name}', skipping.", file=sys.stderr)

    conn.close()
    return 0


def _inspect_contig_feature(conn, contig_id, contig_name, feature_name, fid,
                            region_start, region_end, zoom_bin_size, header_only,
                            mag_name="", mag_offset=0):
    """Inspect a contig-level feature from Contig_blob / Contig_blob_chunk."""
    row = conn.execute(
        "SELECT Zoom_data FROM Contig_blob WHERE Contig_id = ? AND Feature_id = ?",
        [contig_id, fid],
    ).fetchone()
    if row is None:
        return
    zoom_blob = bytes(row[0])

    if header_only:
        _print_zoom_header(zoom_blob, contig_name, feature_name)
        return

    if zoom_bin_size is not None:
        _print_zoom_rows(zoom_blob, zoom_bin_size, contig_name, "", feature_name, region_start, region_end,
                         mag_name=mag_name, mag_offset=mag_offset)
        return

    scale_div = get_scale_from_zoom_blob(zoom_blob)
    chunk_rows = _fetch_contig_chunks(conn, contig_id, fid, region_start, region_end, feature_name)
    if not chunk_rows:
        return

    data = decode_raw_chunks(chunk_rows, scale_div)
    x_arr = data["x"]
    y_arr = data["y"]

    gc_window = 500 if feature_name == "gc_content" else (1000 if feature_name == "gc_skew" else 1)
    if gc_window > 1:
        x_arr = x_arr * gc_window

    _print_rows(contig_name, "", feature_name, x_arr, y_arr, region_start, region_end,
                window_size=gc_window, mag_name=mag_name, mag_offset=mag_offset)


def _inspect_sample_feature(conn, contig_id, contig_name, sample_name, sample_id,
                            feature_name, fid,
                            region_start, region_end, zoom_bin_size, header_only,
                            mag_name="", mag_offset=0):
    """Inspect a sample-level feature from Feature_blob / Feature_blob_chunk."""
    row = conn.execute(
        "SELECT Zoom_data FROM Feature_blob WHERE Contig_id = ? AND Sample_id = ? AND Feature_id = ?",
        [contig_id, sample_id, fid],
    ).fetchone()
    if row is None:
        return
    zoom_blob = bytes(row[0])

    if header_only:
        _print_zoom_header(zoom_blob, contig_name, feature_name, sample_name)
        return

    if zoom_bin_size is not None:
        _print_zoom_rows(zoom_blob, zoom_bin_size, contig_name, sample_name, feature_name, region_start, region_end,
                         mag_name=mag_name, mag_offset=mag_offset)
        return

    scale_div = get_scale_from_zoom_blob(zoom_blob)
    sparse = is_sparse_zoom_blob(zoom_blob)
    chunk_rows = _fetch_feature_chunks(conn, contig_id, sample_id, fid, region_start, region_end)
    if not chunk_rows:
        return

    if sparse:
        data = decode_raw_sparse_chunks(chunk_rows, scale_div)
    else:
        data = decode_raw_chunks(chunk_rows, scale_div)

    _print_rows(contig_name, sample_name, feature_name, data["x"], data["y"], region_start, region_end,
                mag_name=mag_name, mag_offset=mag_offset)


def _fetch_contig_chunks(conn, contig_id, fid, region_start, region_end, feature_name):
    """Fetch Contig_blob_chunk rows, optionally limited to region."""
    gc_window = 500 if feature_name == "gc_content" else (1000 if feature_name == "gc_skew" else 1)
    if region_start is not None:
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
    """Print zoom blob header info (diagnostic, to stderr)."""
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


def _print_zoom_rows(zoom_blob, bin_size, contig_name, sample_name, feature_name,
                     region_start, region_end, mag_name="", mag_offset=0):
    """Print zoom summary bins as RLE TSV rows."""
    data = decode_zoom_standalone(zoom_blob, bin_size)
    if data is None:
        return

    bin_starts = data["bin_start"]
    bin_ends = data["bin_end"]
    if region_start is not None:
        mask = (bin_starts >= (region_start - 1)) & (bin_starts <= (region_end - 1))
        bin_starts = bin_starts[mask]
        bin_ends = bin_ends[mask]
        if "mean" in data:
            data["mean"] = data["mean"][mask]

    if len(bin_starts) == 0:
        return

    means = np.round(data["mean"], 4)
    pos_starts = (bin_starts + 1 + mag_offset).astype(np.int64)
    pos_ends = (bin_ends + 1 + mag_offset).astype(np.int64)

    # RLE: merge consecutive bins with same value
    _write_rle(contig_name, sample_name, feature_name, pos_starts, pos_ends, means, mag_name=mag_name)


def _print_rows(contig_name, sample_name, feature_name, x_arr, y_arr,
                region_start, region_end, window_size=1, mag_name="", mag_offset=0):
    """Apply region filter and print RLE TSV rows.

    Args:
        window_size: For windowed features (gc_content=500, gc_skew=1000),
                     each entry spans window_size bp. position_end = position_start + window_size - 1.
        mag_name: MAG the contig belongs to (empty string when not using --mag mode).
        mag_offset: Offset of this contig within the MAG (0 when not in MAG mode).
                    Added to all output positions so they are in MAG-global coordinates.
    """
    if region_start is not None:
        mask = (x_arr >= (region_start - 1)) & (x_arr <= (region_end - 1))
        x_arr = x_arr[mask]
        y_arr = y_arr[mask]

    n = len(x_arr)
    if n == 0:
        return

    y_rounded = np.round(y_arr, 4)
    pos_starts = (x_arr + 1 + mag_offset).astype(np.int64)  # 0-based → 1-based + MAG offset
    pos_ends = (x_arr + window_size + mag_offset).astype(np.int64)  # end of each window/position

    # Find runs: same value AND consecutive entries (next start == prev end + 1)
    if n == 1:
        _write_rle(contig_name, sample_name, feature_name,
                   pos_starts, pos_ends, y_rounded, mag_name=mag_name)
        return

    value_changes = np.diff(y_rounded) != 0
    # Entries are adjacent if the next start immediately follows the previous end
    gaps = pos_starts[1:] != (pos_ends[:-1] + 1)
    breaks = value_changes | gaps

    change_idx = np.where(breaks)[0] + 1
    run_starts = np.concatenate([[0], change_idx])
    run_ends = np.concatenate([change_idx - 1, [n - 1]])

    rle_pos_starts = pos_starts[run_starts]
    rle_pos_ends = pos_ends[run_ends]
    run_values = y_rounded[run_starts]

    _write_rle(contig_name, sample_name, feature_name, rle_pos_starts, rle_pos_ends, run_values,
               mag_name=mag_name)


def _write_rle(contig_name, sample_name, feature_name, pos_starts, pos_ends, values, mag_name=""):
    """Write RLE rows to stdout in bulk."""
    n = len(pos_starts)
    if n == 0:
        return

    if mag_name:
        prefix = f"{mag_name}\t{contig_name}\t{sample_name}\t{feature_name}\t"
    else:
        prefix = f"{contig_name}\t{sample_name}\t{feature_name}\t"
    prefix_b = prefix.encode()

    # Format numeric columns with numpy, then prepend prefix
    buf = _io.BytesIO()
    np.savetxt(buf, np.column_stack([pos_starts, pos_ends, values]),
               fmt=['%d', '%d', '%.4f'], delimiter='\t')
    raw_lines = buf.getvalue().rstrip(b'\n').split(b'\n')
    sys.stdout.buffer.write(b'\n'.join(prefix_b + line for line in raw_lines) + b'\n')
