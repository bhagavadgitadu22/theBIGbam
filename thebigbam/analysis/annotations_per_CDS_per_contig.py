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
    get_blob_scale, get_chunk_size, get_gc_window_size,
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

CONTIG_BATCH_SIZE = 10_000


def _prefetch_batch_gc(conn, contig_id_batch, gc_feature_id):
    """Prefetch GC blob data for a batch of contigs in 2 SQL queries."""
    if gc_feature_id is None:
        return {}, {}

    placeholders = ",".join("?" * len(contig_id_batch))
    params = list(contig_id_batch) + [gc_feature_id]

    zoom_by_contig = {}
    for row in conn.execute(
        f"SELECT Contig_id, Zoom_data FROM Contig_blob "
        f"WHERE Contig_id IN ({placeholders}) AND Feature_id=?",
        params,
    ).fetchall():
        zoom_by_contig[row[0]] = bytes(row[1])

    chunks_by_contig = defaultdict(list)
    for row in conn.execute(
        f"SELECT Contig_id, Chunk_idx, Data FROM Contig_blob_chunk "
        f"WHERE Contig_id IN ({placeholders}) AND Feature_id=? "
        f"ORDER BY Contig_id, Chunk_idx",
        params,
    ).fetchall():
        chunks_by_contig[row[0]].append((row[1], row[2]))

    return zoom_by_contig, chunks_by_contig


def _decode_gc_array(contig_id, contig_length, zoom_by_contig, chunks_by_contig,
                     scale_div, chunk_sz, window_size):
    """Decode GC content array for one contig from prefetched data."""
    zoom_blob = zoom_by_contig.get(contig_id)
    chunk_rows = chunks_by_contig.get(contig_id)
    if zoom_blob is None or not chunk_rows:
        return None

    data = decode_raw_chunks(chunk_rows, scale_div, chunk_sz)

    arr = np.full(contig_length, np.nan, dtype=np.float64)
    for i in range(len(data["x"])):
        start = int(data["x"][i]) * window_size
        end = min(start + window_size, contig_length)
        arr[start:end] = data["y"][i]
    return arr


DESCRIPTION = """\
Export all annotation info per CDS from a theBIGbam database.

Output columns (one row per CDS):
- mag_name: (MAG mode only) MAG the contig belongs to
- contig_name: contig the CDS belongs to
- gene_name: stable identifier <contig_name>_tbb_<N>, numbered per contig by start position
- start: 1-based start position of the CDS on the contig
- end: 1-based end position of the CDS on the contig
- gene_length: length of the CDS in base pairs (end - start + 1)
- strand: coding strand (+ or -)
- main_isoform: whether this CDS is the main isoform at its locus
- contig_gc_content: mean GC content of the contig (percentage)
- gene_gc_content: mean GC content across the CDS (percentage)
- <qualifier columns>: one column per distinct qualifier Key found for CDS features \
(e.g. product, gene, locus_tag)
"""


def add_args(parser):
    """Register arguments on an existing subparser."""
    import argparse

    parser.formatter_class = argparse.RawDescriptionHelpFormatter
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

    # Get all contigs that have CDS features
    contigs_with_cds = conn.execute(
        """
        SELECT DISTINCT c.Contig_id, c.Contig_name, c.Contig_length, c.GC_mean
        FROM Contig_annotation_core ca
        JOIN Contig c ON ca.Contig_id = c.Contig_id
        WHERE ca."Type" = 'CDS'
        ORDER BY c.Contig_id
        """
    ).fetchall()

    if not contigs_with_cds:
        print("No CDS features found in database.", file=sys.stderr, flush=True)
        conn.close()
        sys.exit(0)

    total_contigs = len(contigs_with_cds)
    print(f"Processing {total_contigs} contigs ({len(qualifier_keys)} qualifier columns)...", file=sys.stderr, flush=True)

    columns = (["mag_name"] + FIXED_COLUMNS + qualifier_keys
               if mag_mode else FIXED_COLUMNS + qualifier_keys)

    total_cds = 0
    gc_feature_id = feature_name_to_id("gc_content", conn)
    gc_scale_div = get_blob_scale(conn, "gc_content")
    gc_chunk_sz = get_chunk_size(conn)
    gc_window_sz = get_gc_window_size(conn, "gc_content")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)

        for batch_start in range(0, total_contigs, CONTIG_BATCH_SIZE):
            batch = contigs_with_cds[batch_start:batch_start + CONTIG_BATCH_SIZE]
            batch_ids = [r[0] for r in batch]
            placeholders = ",".join("?" * len(batch_ids))

            batch_cds = conn.execute(
                f"""
                SELECT ca.Contig_id, ca.Annotation_id, ca."Start", ca."End",
                       ca.Strand, ca.Main_isoform
                FROM Contig_annotation_core ca
                WHERE ca.Contig_id IN ({placeholders}) AND ca."Type" = 'CDS'
                ORDER BY ca.Contig_id, ca."Start"
                """,
                batch_ids,
            ).fetchall()

            ann_ids = [r[1] for r in batch_cds]
            quals = defaultdict(dict)
            if ann_ids:
                ann_ph = ",".join("?" * len(ann_ids))
                for ann_id, key, value in conn.execute(
                    f"""
                    SELECT aq.Annotation_id, aq."Key", aq."Value"
                    FROM Annotation_qualifier aq
                    WHERE aq.Annotation_id IN ({ann_ph})
                    """,
                    ann_ids,
                ).fetchall():
                    quals[ann_id][key] = value

            gc_zoom, gc_chunks = _prefetch_batch_gc(conn, batch_ids, gc_feature_id)

            cds_by_contig = defaultdict(list)
            for row in batch_cds:
                cds_by_contig[row[0]].append(row[1:])

            for contig_id, contig_name, contig_length, contig_gc in batch:
                contig_cds = cds_by_contig.get(contig_id)
                if not contig_cds:
                    continue

                gc_arr = _decode_gc_array(contig_id, contig_length, gc_zoom, gc_chunks,
                                         gc_scale_div, gc_chunk_sz, gc_window_sz)
                mag_name = contig_to_mag.get(contig_name, "") if mag_mode else None

                for gene_idx, (ann_id, start, end, strand, main_isoform) in enumerate(contig_cds, start=1):
                    gene_length = end - start + 1
                    gene_name = f"{contig_name}_tbb_{gene_idx}"

                    gene_gc = None
                    if gc_arr is not None:
                        s0 = start - 1
                        e0 = end
                        gene_slice = gc_arr[s0:e0]
                        valid = ~np.isnan(gene_slice)
                        if np.any(valid):
                            gene_gc = round(float(np.mean(gene_slice[valid])), 1)

                    fixed = [contig_name, gene_name, start, end, gene_length,
                             strand, main_isoform, contig_gc, gene_gc]
                    if mag_mode:
                        fixed = [mag_name] + fixed
                    dynamic = [quals[ann_id].get(key, "") for key in qualifier_keys]

                    writer.writerow(fixed + dynamic)

                total_cds += len(contig_cds)

            processed = min(batch_start + CONTIG_BATCH_SIZE, total_contigs)
            print(f"  Processed {processed}/{total_contigs} contigs ({total_cds} CDS)...", file=sys.stderr, flush=True)

    conn.close()

    print(f"Wrote {total_cds} CDS to {args.output} ({len(qualifier_keys)} qualifier columns)", file=sys.stderr, flush=True)


def main():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
    )
    add_args(parser)
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
