#!/usr/bin/env python3
"""Export per-gene feature signals from a theBIGbam database.

Produces a TSV with one row per (sample, gene) combination, containing
coverage, GC, mismatch, and structural variant metrics.

Sparse feature values (mismatches, insertions, deletions, clippings) are
stored as prevalence per mille in blobs (ValueScale::Times1000). After
decoding, values are fractional prevalences. For each event type we export
both the number of affected positions and the sum of prevalences across
those positions.

All feature data is decoded from Feature_blob / Feature_blob_chunk tables
using Rust-accelerated decompression. SQL queries are batched: 2 per
(contig, sample) instead of per-feature.
"""

import argparse
import csv
import sys
from collections import defaultdict

import numpy as np
import duckdb

from thebigbam.database.blob_decoder import (
    feature_name_to_id, contig_blob_name_to_id,
    decode_raw_chunks, decode_raw_sparse_chunks,
    get_scale_from_zoom_blob, is_sparse_zoom_blob,
)


COLUMNS = [
    "sample_name",
    "contig_name",
    "contig_length",
    "gene_start",
    "gene_end",
    "gene_length",
    "gene_product",
    "gene_function",
    "phrog_id",
    "contig_gc_content",
    "gene_gc_content",
    "contig_aligned_fraction",
    "contig_coverage_median",
    "gene_coverage_median",
    "coverage_ratio",
    "gene_secondary_coverage_median",
    "mismatches_positions",
    "mismatches_prevalence",
    "synonymous_positions",
    "synonymous_prevalence",
    "nonsynonymous_positions",
    "nonsynonymous_prevalence",
    "s_sites",
    "n_sites",
    "dnds_ratio",
    "insertions_positions",
    "insertions_prevalence",
    "deletions_positions",
    "deletions_prevalence",
    "left_clippings_positions",
    "left_clippings_prevalence",
    "right_clippings_positions",
    "right_clippings_prevalence",
]

# Features to decode per (contig, sample)
_DENSE_FEATURES = ["primary_reads", "secondary_reads"]
_SPARSE_FEATURES = ["mismatches", "insertions", "deletions", "left_clippings", "right_clippings"]

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_BASES = "ACGT"


def _compute_sn_sites(nuc_seq):
    """Compute (S_sites, N_sites) from a nucleotide sequence."""
    seq = nuc_seq.upper()
    syn = 0
    nonsyn = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) != 3 or not all(b in _BASES for b in codon):
            continue
        orig_aa = _CODON_TABLE.get(codon)
        if orig_aa is None or orig_aa == "*":
            continue
        for pos in range(3):
            for alt in _BASES:
                if alt == codon[pos]:
                    continue
                mutant = codon[:pos] + alt + codon[pos + 1:]
                mut_aa = _CODON_TABLE.get(mutant)
                if mut_aa is not None:
                    if mut_aa == orig_aa:
                        syn += 1
                    else:
                        nonsyn += 1
    return (syn, nonsyn)


def _column_exists(conn, table_name, column_name):
    """Check if a column exists in a table."""
    try:
        conn.execute(f'SELECT "{column_name}" FROM {table_name} LIMIT 0')
        return True
    except (duckdb.CatalogException, duckdb.BinderException):
        return False


def _table_exists(conn, table_name):
    """Check if a table exists in the database."""
    return (
        conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = ?",
            [table_name],
        ).fetchone()
        is not None
    )


# ============================================================================
# Blob decoding — batched per (contig, sample)
# ============================================================================

def _resolve_feature_ids(conn):
    """Build {feature_id: feature_name} map for all needed features."""
    id_to_name = {}
    name_to_id = {}
    for name in _DENSE_FEATURES + _SPARSE_FEATURES:
        fid = feature_name_to_id(name, conn)
        if fid is not None:
            id_to_name[fid] = name
            name_to_id[name] = fid
    return id_to_name, name_to_id


def _load_all_features(conn, contig_id, sample_id, contig_length, id_to_name, name_to_id):
    """Decode all sample-level features for one (contig, sample) in 2 SQL queries.

    Returns dict: {feature_name: ndarray (dense) or sparse_dict (sparse)}
    Missing features get zeros/empty arrays.
    """
    # Query 1: all zoom headers
    zoom_map = {}
    for row in conn.execute(
        "SELECT Feature_id, Zoom_data FROM Feature_blob "
        "WHERE Contig_id=? AND Sample_id=?", [contig_id, sample_id]
    ).fetchall():
        fid = row[0]
        if fid in id_to_name:
            zoom_map[fid] = bytes(row[1])

    # Query 2: all chunks
    chunks_by_fid = defaultdict(list)
    for row in conn.execute(
        "SELECT Feature_id, Chunk_idx, Data FROM Feature_blob_chunk "
        "WHERE Contig_id=? AND Sample_id=? ORDER BY Feature_id, Chunk_idx",
        [contig_id, sample_id]
    ).fetchall():
        fid = row[0]
        if fid in id_to_name:
            chunks_by_fid[fid].append((row[1], row[2]))

    result = {}

    # Decode dense features → full-length numpy arrays
    for name in _DENSE_FEATURES:
        fid = name_to_id.get(name)
        if fid is None or fid not in zoom_map or fid not in chunks_by_fid:
            result[name] = np.zeros(contig_length, dtype=np.float64)
            continue
        scale_div = get_scale_from_zoom_blob(zoom_map[fid])
        data = decode_raw_chunks(chunks_by_fid[fid], scale_div)
        arr = np.zeros(contig_length, dtype=np.float64)
        x, y = data["x"], data["y"]
        valid = x < contig_length
        arr[x[valid].astype(np.intp)] = y[valid]
        result[name] = arr

    # Decode sparse features → position/value/metadata dicts
    _empty_sparse = {"x": np.array([], dtype=np.uint32), "y": np.array([], dtype=np.float64)}
    for name in _SPARSE_FEATURES:
        fid = name_to_id.get(name)
        if fid is None or fid not in zoom_map or fid not in chunks_by_fid:
            result[name] = _empty_sparse
            continue
        scale_div = get_scale_from_zoom_blob(zoom_map[fid])
        result[name] = decode_raw_sparse_chunks(chunks_by_fid[fid], scale_div)

    return result


def _load_contig_gc(conn, contig_id, contig_length):
    """Decode GC content from Contig_blob for one contig. Returns full-length array."""
    fid = contig_blob_name_to_id("gc_content")
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

    # x = window indices (0-based), y = GC values; 500bp windows
    window_size = 500
    arr = np.full(contig_length, np.nan, dtype=np.float64)
    for i in range(len(data["x"])):
        start = int(data["x"][i]) * window_size
        end = min(start + window_size, contig_length)
        arr[start:end] = data["y"][i]
    return arr


# ============================================================================
# Sparse feature helpers
# ============================================================================

def _sparse_stats(sparse_dict, start_0, end_0):
    """Count positions and sum prevalence for a sparse feature within [start_0, end_0)."""
    x = sparse_dict["x"]
    if len(x) == 0:
        return None, None
    mask = (x >= start_0) & (x < end_0)
    count = int(np.sum(mask))
    if count == 0:
        return None, None
    return count, round(float(np.sum(sparse_dict["y"][mask])), 4)


def _mismatch_stats(mm_dict, start_0, end_0):
    """Mismatch positions/prevalence with synonymous/nonsynonymous split.

    Returns (total_pos, total_prev, syn_pos, syn_prev, nonsyn_pos, nonsyn_prev).
    """
    x = mm_dict["x"]
    if len(x) == 0:
        return None, None, None, None, None, None

    mask = (x >= start_0) & (x < end_0)
    total_pos = int(np.sum(mask))
    if total_pos == 0:
        return None, None, None, None, None, None

    total_prev = round(float(np.sum(mm_dict["y"][mask])), 4)

    cats = mm_dict.get("codon_category")
    if cats is None:
        return total_pos, total_prev, None, None, None, None

    cat_arr = np.array(cats) if not isinstance(cats, np.ndarray) else cats
    syn_mask = mask & (cat_arr == "Synonymous")
    nonsyn_mask = mask & (cat_arr == "Non-synonymous")
    syn_pos = int(np.sum(syn_mask))
    syn_prev = round(float(np.sum(mm_dict["y"][syn_mask])), 4) if syn_pos else None
    nonsyn_pos = int(np.sum(nonsyn_mask))
    nonsyn_prev = round(float(np.sum(mm_dict["y"][nonsyn_mask])), 4) if nonsyn_pos else None

    return total_pos, total_prev, syn_pos, syn_prev, nonsyn_pos, nonsyn_prev


# ============================================================================
# Main processing
# ============================================================================

def process_sample(conn, sample_id, sample_name, contig_info, cov_map, af_map,
                   gc_arrays, id_to_name, name_to_id, has_sn_sites):
    """Process all genes for one sample. Yields row tuples."""
    contig_ids = list(cov_map.keys())
    if not contig_ids:
        return

    id_list = ",".join(str(c) for c in contig_ids)

    # Fetch genes with Phrog from qualifier table
    if has_sn_sites:
        genes = conn.execute(
            f"""
            SELECT ca.Contig_id, ca."Start", ca."End", ca."End" - ca."Start" + 1,
                   ca.Product, ca."Function",
                   aq_phrog."Value",
                   ca.S_sites, ca.N_sites
            FROM Contig_annotation ca
            LEFT JOIN Annotation_qualifier aq_phrog
              ON ca.Annotation_id = aq_phrog.Annotation_id AND aq_phrog."Key" = 'phrog'
            WHERE ca.Contig_id IN ({id_list})
            ORDER BY ca.Contig_id, ca."Start"
            """
        ).fetchall()
    else:
        raw = conn.execute(
            f"""
            SELECT ca.Contig_id, ca."Start", ca."End", ca."End" - ca."Start" + 1,
                   ca.Product, ca."Function",
                   aq_phrog."Value",
                   ca.Nucleotide_sequence
            FROM Contig_annotation ca
            LEFT JOIN Annotation_qualifier aq_phrog
              ON ca.Annotation_id = aq_phrog.Annotation_id AND aq_phrog."Key" = 'phrog'
            WHERE ca.Contig_id IN ({id_list})
            ORDER BY ca.Contig_id, ca."Start"
            """
        ).fetchall()
        genes = []
        for r in raw:
            nuc_seq = r[7]
            if nuc_seq:
                s, n = _compute_sn_sites(nuc_seq)
            else:
                s, n = None, None
            genes.append((*r[:7], s, n))

    if not genes:
        return

    # Group genes by contig
    genes_by_contig = defaultdict(list)
    for g in genes:
        genes_by_contig[g[0]].append(g)

    # Process each contig: batch-decode all features, then iterate genes
    for contig_id in contig_ids:
        contig_genes = genes_by_contig.get(contig_id)
        if not contig_genes:
            continue

        contig_name, contig_length, contig_gc = contig_info[contig_id]
        contig_cov_med = cov_map[contig_id]
        contig_af = af_map[contig_id]

        # 2 SQL queries → all features decoded
        features = _load_all_features(
            conn, contig_id, sample_id, contig_length, id_to_name, name_to_id
        )

        gc_arr = gc_arrays.get(contig_id)

        for g in contig_genes:
            gene_start, gene_end, gene_length = g[1], g[2], g[3]
            s_sites, n_sites = g[7], g[8]

            # 0-based half-open slice
            s0 = gene_start - 1
            e0 = gene_end

            # Dense: gene median
            gene_cov_med = float(np.median(features["primary_reads"][s0:e0]))
            gene_sec_med = float(np.median(features["secondary_reads"][s0:e0]))

            coverage_ratio = None
            if contig_cov_med is not None and contig_cov_med > 0:
                coverage_ratio = round(gene_cov_med / contig_cov_med, 4)

            # GC per gene
            gene_gc = None
            if gc_arr is not None:
                gene_slice = gc_arr[s0:e0]
                valid = ~np.isnan(gene_slice)
                if np.any(valid):
                    gene_gc = round(float(np.mean(gene_slice[valid])), 1)

            # Mismatches with syn/nonsyn split
            mm = _mismatch_stats(features["mismatches"], s0, e0)
            total_pos, total_prev = mm[0], mm[1]
            syn_pos, syn_prev = mm[2], mm[3]
            nonsyn_pos, nonsyn_prev = mm[4], mm[5]

            # dN/dS
            dnds = None
            if (syn_prev is not None and nonsyn_prev is not None
                    and syn_prev > 0 and n_sites and s_sites):
                dnds = round((nonsyn_prev * s_sites) / (syn_prev * n_sites), 4)

            # Sparse event features
            ins_pos, ins_prev = _sparse_stats(features["insertions"], s0, e0)
            del_pos, del_prev = _sparse_stats(features["deletions"], s0, e0)
            lc_pos, lc_prev = _sparse_stats(features["left_clippings"], s0, e0)
            rc_pos, rc_prev = _sparse_stats(features["right_clippings"], s0, e0)

            yield (
                sample_name,
                contig_name,
                contig_length,
                gene_start,
                gene_end,
                gene_length,
                g[4],  # product
                g[5],  # function
                g[6],  # phrog
                contig_gc,
                gene_gc,
                contig_af,
                contig_cov_med,
                gene_cov_med,
                coverage_ratio,
                gene_sec_med,
                total_pos,
                total_prev,
                syn_pos,
                syn_prev,
                nonsyn_pos,
                nonsyn_prev,
                s_sites,
                n_sites,
                dnds,
                ins_pos,
                ins_prev,
                del_pos,
                del_prev,
                lc_pos,
                lc_prev,
                rc_pos,
                rc_prev,
            )


def main():
    parser = argparse.ArgumentParser(
        description="Export per-gene feature signals from a theBIGbam database."
    )
    parser.add_argument("db", help="Path to the theBIGbam DuckDB database")
    parser.add_argument("output", help="Path to the output TSV file")
    parser.add_argument(
        "--min_aligned_fraction",
        type=float,
        default=50,
        help="Minimum aligned fraction percentage (0-100) to include a contig-sample pair (default: 50)",
    )
    parser.add_argument(
        "--min_coverage_depth",
        type=float,
        default=10,
        help="Minimum median coverage depth to include a contig-sample pair (default: 10)",
    )
    args = parser.parse_args()

    if not 0 <= args.min_aligned_fraction <= 100:
        print("Error: --min_aligned_fraction must be between 0 and 100.", file=sys.stderr)
        sys.exit(1)

    conn = duckdb.connect(args.db, read_only=True)

    if not _table_exists(conn, "Contig_annotation"):
        print("Error: No Contig_annotation table. Run pharokka annotation first.", file=sys.stderr)
        sys.exit(1)

    if not _table_exists(conn, "Sample"):
        print("Error: No Sample table. Database has no BAM data.", file=sys.stderr)
        sys.exit(1)

    samples = conn.execute(
        "SELECT Sample_id, Sample_name FROM Sample ORDER BY Sample_id"
    ).fetchall()
    if not samples:
        print("Error: No samples found in database.", file=sys.stderr)
        sys.exit(1)

    # Contig info (sample-independent)
    contigs = conn.execute(
        "SELECT Contig_id, Contig_name, Contig_length, GC_mean FROM Contig ORDER BY Contig_id"
    ).fetchall()
    contig_info = {r[0]: (r[1], r[2], r[3]) for r in contigs}

    # Resolve feature IDs once
    id_to_name, name_to_id = _resolve_feature_ids(conn)

    has_sn_sites = _column_exists(conn, "Contig_annotation", "S_sites")
    has_coverage = _table_exists(conn, "Coverage")

    # Precompute GC arrays per contig (sample-independent)
    gc_arrays = {}
    for contig_id, (_, contig_length, _) in contig_info.items():
        gc_arr = _load_contig_gc(conn, contig_id, contig_length)
        if gc_arr is not None:
            gc_arrays[contig_id] = gc_arr

    print(
        f"Processing {len(samples)} sample(s), "
        f"{len(name_to_id)} feature(s) available...",
        file=sys.stderr,
    )
    print(
        f"Filters: aligned_fraction >= {args.min_aligned_fraction}%, "
        f"coverage_depth >= {args.min_coverage_depth}",
        file=sys.stderr,
    )

    # Thresholds in stored units (×10)
    af_threshold = args.min_aligned_fraction * 10
    cov_threshold = args.min_coverage_depth * 10

    total_rows = 0

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(COLUMNS)

        for sample_id, sample_name in samples:
            if has_coverage:
                qualifying = conn.execute(
                    """
                    SELECT Contig_id, Aligned_fraction_percentage / 10.0,
                           Coverage_median / 10.0
                    FROM Coverage
                    WHERE Sample_id = ?
                      AND Aligned_fraction_percentage >= ?
                      AND Coverage_median >= ?
                    """,
                    [sample_id, af_threshold, cov_threshold],
                ).fetchall()
            else:
                qualifying = []

            af_map = {r[0]: r[1] for r in qualifying}
            cov_map = {r[0]: r[2] for r in qualifying}

            sample_rows = 0
            for row in process_sample(
                conn, sample_id, sample_name, contig_info, cov_map, af_map,
                gc_arrays, id_to_name, name_to_id, has_sn_sites,
            ):
                writer.writerow(row)
                sample_rows += 1

            total_rows += sample_rows
            print(f"  Done: {sample_name} ({sample_rows} genes)", file=sys.stderr)

    conn.close()
    print(f"Wrote {total_rows} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
