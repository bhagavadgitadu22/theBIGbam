#!/usr/bin/env python3
"""Export all annotation info per CDS from a theBIGbam database.

Produces a TSV with fixed columns (contig_name, gene_name, start, end,
gene_length, strand, main_isoform, contig_gc_content, gene_gc_content)
followed by one column per distinct qualifier Key found for CDS features
in Annotation_qualifier.

Gene names follow the <contig_name>_tbb_<N> convention, numbered
sequentially per contig by Start position — consistent with
mapping_patterns_per_CDS_per_contig_per_sample.py.
"""

import argparse
import csv
import sys
from collections import defaultdict

import numpy as np
import duckdb

from thebigbam.database.blob_decoder import (
    feature_name_to_id,
    decode_raw_chunks,
    get_scale_from_zoom_blob,
    GC_CONTENT_WINDOW_SIZE,
)


FIXED_COLUMNS = [
    "contig_name",
    "gene_name",
    "start",
    "end",
    "gene_length",
    "strand",
    "main_isoform",
    "contig_gc_content",
    "gene_gc_content",
]


def _load_contig_gc(conn, contig_id, contig_length):
    """Decode GC content from Contig_blob for one contig. Returns full-length array."""
    fid = feature_name_to_id("gc_content", conn)
    zoom_row = conn.execute(
        "SELECT Zoom_data FROM Contig_blob WHERE Contig_id=? AND Feature_id=?",
        [contig_id, fid]
    ).fetchone()
    if zoom_row is None:
        return None

    chunk_rows = conn.execute(
        "SELECT Chunk_idx, Data FROM Contig_blob_chunk "
        "WHERE Contig_id=? AND Feature_id=? ORDER BY Chunk_idx",
        [contig_id, fid]
    ).fetchall()
    if not chunk_rows:
        return None

    zoom_blob = bytes(zoom_row[0])
    scale_div = get_scale_from_zoom_blob(zoom_blob)
    chunk_rows = [(r[0], r[1]) for r in chunk_rows]
    data = decode_raw_chunks(chunk_rows, scale_div)

    window_size = GC_CONTENT_WINDOW_SIZE
    arr = np.full(contig_length, np.nan, dtype=np.float64)
    for i in range(len(data["x"])):
        start = int(data["x"][i]) * window_size
        end = min(start + window_size, contig_length)
        arr[start:end] = data["y"][i]
    return arr


DESCRIPTION = "Export all annotation info per CDS from a theBIGbam database."


def add_args(parser):
    """Register arguments on an existing subparser."""
    parser.add_argument("--db", required=True, help="Path to the theBIGbam DuckDB database")
    parser.add_argument("--output", required=True, help="Path to the output TSV file")


def run(args):
    """Entry point called by CLI dispatcher."""
    conn = duckdb.connect(args.db, read_only=True)

    # Detect MAG mode and build contig→MAG map
    from thebigbam.database.database_getters import is_mag_mode
    mag_mode = is_mag_mode(conn)
    contig_to_mag = {}
    if mag_mode:
        rows = conn.execute(
            "SELECT c.Contig_name, mg.MAG_name "
            "FROM Contig c "
            "JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id "
            "JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
        ).fetchall()
        contig_to_mag = {cname: mname for cname, mname in rows}

    # Discover all qualifier keys used for CDS features
    qualifier_keys = [
        r[0] for r in conn.execute(
            """
            SELECT DISTINCT aq."Key"
            FROM Annotation_qualifier aq
            JOIN Contig_annotation_core ca ON aq.Annotation_id = ca.Annotation_id
            WHERE ca."Type" = 'CDS'
            ORDER BY aq."Key"
            """
        ).fetchall()
    ]

    # Get all CDS with core fields + contig info
    cds_rows = conn.execute(
        """
        SELECT ca.Annotation_id, c.Contig_name, ca.Contig_id,
               ca."Start", ca."End", ca.Strand, ca.Main_isoform,
               c.Contig_length, c.GC_mean
        FROM Contig_annotation_core ca
        JOIN Contig c ON ca.Contig_id = c.Contig_id
        WHERE ca."Type" = 'CDS'
        ORDER BY ca.Contig_id, ca."Start"
        """
    ).fetchall()

    if not cds_rows:
        print("No CDS features found in database.", file=sys.stderr)
        conn.close()
        sys.exit(0)

    # Get all qualifiers for CDS annotations in one query
    quals = defaultdict(dict)
    for ann_id, key, value in conn.execute(
        """
        SELECT aq.Annotation_id, aq."Key", aq."Value"
        FROM Annotation_qualifier aq
        JOIN Contig_annotation_core ca ON aq.Annotation_id = ca.Annotation_id
        WHERE ca."Type" = 'CDS'
        """
    ).fetchall():
        quals[ann_id][key] = value

    # Precompute GC arrays per contig
    contig_lengths = {}
    for row in cds_rows:
        contig_lengths[row[2]] = row[7]  # contig_id → contig_length

    gc_arrays = {}
    for contig_id, contig_length in contig_lengths.items():
        gc_arr = _load_contig_gc(conn, contig_id, contig_length)
        if gc_arr is not None:
            gc_arrays[contig_id] = gc_arr

    conn.close()

    # Build gene names: <contig_name>_tbb_<N> per contig
    counter = defaultdict(int)
    gene_names = []
    for row in cds_rows:
        contig_id = row[2]
        contig_name = row[1]
        counter[contig_id] += 1
        gene_names.append(f"{contig_name}_tbb_{counter[contig_id]}")

    # Write TSV
    columns = (["mag_name"] + FIXED_COLUMNS + qualifier_keys
               if mag_mode else FIXED_COLUMNS + qualifier_keys)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)

        for i, row in enumerate(cds_rows):
            ann_id = row[0]
            contig_name = row[1]
            contig_id = row[2]
            start = row[3]
            end = row[4]
            strand = row[5]
            main_isoform = row[6]
            contig_gc = row[8]  # GC_mean from Contig table

            gene_length = end - start + 1

            # Gene GC from blob
            gene_gc = None
            gc_arr = gc_arrays.get(contig_id)
            if gc_arr is not None:
                s0 = start - 1
                e0 = end
                gene_slice = gc_arr[s0:e0]
                valid = ~np.isnan(gene_slice)
                if np.any(valid):
                    gene_gc = round(float(np.mean(gene_slice[valid])), 1)

            fixed = [contig_name, gene_names[i], start, end, gene_length,
                     strand, main_isoform, contig_gc, gene_gc]
            if mag_mode:
                fixed = [contig_to_mag.get(contig_name, "")] + fixed
            dynamic = [quals[ann_id].get(key, "") for key in qualifier_keys]

            writer.writerow(fixed + dynamic)

    print(f"Wrote {len(cds_rows)} CDS to {args.output} ({len(qualifier_keys)} qualifier columns)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
    )
    add_args(parser)
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
