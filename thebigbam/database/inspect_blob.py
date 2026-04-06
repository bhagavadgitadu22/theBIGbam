"""
CLI command: thebigbam inspect

Decode and display Feature_blob or Contig_blob contents as TSV for debugging and exploration.
"""

import sys

import duckdb

from thebigbam.database.blob_decoder import (
    decode_blob, decode_zoom_level, get_blob_header,
    feature_name_to_id, contig_blob_name_to_id,
    CODON_CATEGORIES,
)


def add_inspect_args(parser):
    parser.add_argument('-d', '--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--feature', required=True, help='Feature name (e.g. mismatches, primary_reads, gc_content)')
    parser.add_argument('--contig', required=True, help='Contig name')
    parser.add_argument('--sample', default=None, help='Sample name (required for sample-level features, omit for contig-level)')
    parser.add_argument('--region', default=None, help='Genomic region to display (e.g. 1000-2000)')
    parser.add_argument('--zoom', type=int, default=None, choices=[0, 1, 2],
                        help='Show zoom level instead of base resolution (0=100bp, 1=1000bp, 2=10000bp)')
    parser.add_argument('--header-only', action='store_true', help='Only show BLOB header info')


def run_inspect(args):
    if args.header_only and args.zoom is not None:
        sys.exit("ERROR: --header-only and --zoom are mutually exclusive.")

    conn = duckdb.connect(args.db, read_only=True)

    # Resolve contig_id
    row = conn.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", [args.contig]).fetchone()
    if row is None:
        sys.exit(f"ERROR: Contig '{args.contig}' not found in database.")
    contig_id = row[0]

    # Detect contig-level vs sample-level feature
    contig_feature_id = contig_blob_name_to_id(args.feature)
    sample_feature_id = feature_name_to_id(args.feature, conn)

    if contig_feature_id is not None:
        # Contig-level feature (Contig_blob)
        if args.sample is not None:
            sys.exit(f"ERROR: '{args.feature}' is a contig-level feature — --sample is not applicable.")

        row = conn.execute(
            "SELECT Data FROM Contig_blob WHERE Contig_id = ? AND Feature_id = ?",
            [contig_id, contig_feature_id],
        ).fetchone()
        conn.close()

        if row is None:
            sys.exit(f"ERROR: No Contig_blob found for contig='{args.contig}', feature='{args.feature}'.")

        blob_bytes = bytes(row[0])
        header = get_blob_header(blob_bytes)
        sample_id = None

    elif sample_feature_id is not None:
        # Sample-level feature (Feature_blob)
        if args.sample is None:
            sys.exit(f"ERROR: '{args.feature}' is a sample-level feature — --sample is required.")

        row = conn.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", [args.sample]).fetchone()
        if row is None:
            sys.exit(f"ERROR: Sample '{args.sample}' not found in database.")
        sample_id = row[0]

        row = conn.execute(
            "SELECT Data FROM Feature_blob WHERE Contig_id = ? AND Sample_id = ? AND Feature_id = ?",
            [contig_id, sample_id, sample_feature_id],
        ).fetchone()
        conn.close()

        if row is None:
            sys.exit(f"ERROR: No BLOB found for contig='{args.contig}', sample='{args.sample}', feature='{args.feature}'.")

        blob_bytes = bytes(row[0])
        header = get_blob_header(blob_bytes)

    else:
        conn.close()
        sys.exit(f"ERROR: Unknown feature '{args.feature}'.")

    if args.header_only:
        _print_header(args, blob_bytes, header, contig_id, sample_id)
        return 0

    # Parse region filter
    region_start, region_end = None, None
    if args.region:
        parts = args.region.split('-')
        if len(parts) != 2:
            sys.exit("ERROR: --region must be in format START-END (e.g. 1000-2000)")
        region_start, region_end = int(parts[0]), int(parts[1])

    if args.zoom is not None:
        _print_zoom(blob_bytes, header, args.zoom, region_start, region_end)
    else:
        _print_base_resolution(blob_bytes, header, region_start, region_end)

    return 0


def _print_header(args, blob_bytes, header, contig_id, sample_id):
    """Print BLOB header info."""
    print(f"Feature: {args.feature}")
    if sample_id is not None:
        print(f"Contig: {args.contig} (id={contig_id}), Sample: {args.sample} (id={sample_id})")
    else:
        print(f"Contig: {args.contig} (id={contig_id})")
    print(f"BLOB size: {len(blob_bytes):,} bytes")
    print(f"Contig length: {header['contig_length']:,}")
    print(f"Sparse: {header['sparse']}, Scale: {header['scale_code']}")
    print(f"Flags: has_stats={header['has_stats']}, has_sequence={header['has_sequence']}, has_codons={header['has_codons']}")
    print(f"Zoom levels: {header['num_zoom_levels']}")


def _print_zoom(blob_bytes, header, level, region_start, region_end):
    """Print zoom level summary bins as TSV."""
    data = decode_zoom_level(blob_bytes, level)
    if data is None:
        sys.exit(f"ERROR: Zoom level {level} not available in this BLOB.")

    bin_starts = data["bin_start"]
    bin_ends = data["bin_end"]

    # Apply region filter
    if region_start is not None:
        mask = (bin_starts >= (region_start - 1)) & (bin_starts <= (region_end - 1))
        bin_starts = bin_starts[mask]
        bin_ends = bin_ends[mask]
        data = {k: v[mask] if hasattr(v, '__len__') and len(v) == len(mask) else v for k, v in data.items()}

    n = len(bin_starts)
    if n == 0:
        return

    is_sparse = header["sparse"]
    if is_sparse:
        print("bin_start\tbin_end\tmax\tmean")
        for i in range(n):
            print(f"{int(bin_starts[i]) + 1}\t{int(bin_ends[i]) + 1}\t{data['max'][i]:.4f}\t{data['mean'][i]:.4f}")
    else:
        print("bin_start\tbin_end\tmean")
        for i in range(n):
            print(f"{int(bin_starts[i]) + 1}\t{int(bin_ends[i]) + 1}\t{data['mean'][i]:.4f}")


def _print_base_resolution(blob_bytes, header, region_start, region_end):
    """Print base-resolution data as TSV."""
    data = decode_blob(blob_bytes)

    x = data.get("x")
    y = data.get("y")
    if x is None or y is None:
        print("No base-resolution data found.", file=sys.stderr)
        return

    # Apply region filter
    if region_start is not None:
        mask = (x >= (region_start - 1)) & (x <= (region_end - 1))
        orig_len = len(x)
        x = x[mask]
        y = y[mask]
        data = {k: v[mask] if hasattr(v, '__len__') and len(v) == orig_len else v for k, v in data.items()}

    n = len(x)
    if n == 0:
        return

    has_stats = "mean" in data
    has_seq = "sequence" in data
    has_prevalence = "sequence_prevalence" in data
    has_codons = "codon_category" in data

    # Build header
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
        parts = [str(int(x[i]) + 1), f"{y[i]:.4f}"]
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
