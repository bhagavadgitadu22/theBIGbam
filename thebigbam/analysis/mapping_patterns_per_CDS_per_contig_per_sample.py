#!/usr/bin/env python3
"""Export per-CDS mapping signals from a theBIGbam database.

Produces a TSV with one row per (sample, CDS) combination, containing
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
    feature_name_to_id,
    decode_raw_chunks, decode_raw_sparse_chunks,
    get_scale_from_zoom_blob,
)


_BASE_COLUMNS = [
    "sample_name",
    "contig_name",
    "gene_name",
    "contig_aligned_fraction",
    "gene_aligned_fraction",
    "contig_coverage_trimmed_mean",
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

_MAG_EXTRA_COLUMNS = {
    "contig_aligned_fraction": "mag_aligned_fraction",
    "contig_coverage_trimmed_mean": "mag_coverage_trimmed_mean",
}


def _build_columns(mag_mode):
    cols = []
    for c in _BASE_COLUMNS:
        if mag_mode and c in _MAG_EXTRA_COLUMNS:
            cols.append(_MAG_EXTRA_COLUMNS[c])
        cols.append(c)
    return cols

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
    try:
        conn.execute(f'SELECT "{column_name}" FROM {table_name} LIMIT 0')
        return True
    except (duckdb.CatalogException, duckdb.BinderException):
        return False


def _table_exists(conn, table_name):
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
    id_to_name = {}
    name_to_id = {}
    for name in _DENSE_FEATURES + _SPARSE_FEATURES:
        fid = feature_name_to_id(name, conn)
        if fid is not None:
            id_to_name[fid] = name
            name_to_id[name] = fid
    return id_to_name, name_to_id


def _load_all_features(conn, contig_id, sample_id, contig_length, id_to_name, name_to_id):
    """Decode all sample-level features for one (contig, sample) in 2 SQL queries."""
    zoom_map = {}
    for row in conn.execute(
        "SELECT Feature_id, Zoom_data FROM Feature_blob "
        "WHERE Contig_id=? AND Sample_id=?", [contig_id, sample_id]
    ).fetchall():
        fid = row[0]
        if fid in id_to_name:
            zoom_map[fid] = bytes(row[1])

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

    _empty_sparse = {"x": np.array([], dtype=np.uint32), "y": np.array([], dtype=np.float64)}
    for name in _SPARSE_FEATURES:
        fid = name_to_id.get(name)
        if fid is None or fid not in zoom_map or fid not in chunks_by_fid:
            result[name] = _empty_sparse
            continue
        scale_div = get_scale_from_zoom_blob(zoom_map[fid])
        result[name] = decode_raw_sparse_chunks(chunks_by_fid[fid], scale_div)

    return result


# ============================================================================
# Sparse feature helpers
# ============================================================================

def _sparse_stats(sparse_dict, start_0, end_0):
    x = sparse_dict["x"]
    if len(x) == 0:
        return 0, 0
    mask = (x >= start_0) & (x < end_0)
    count = int(np.sum(mask))
    if count == 0:
        return 0, 0
    return count, round(float(np.sum(sparse_dict["y"][mask])), 4)


def _mismatch_stats(mm_dict, start_0, end_0):
    x = mm_dict["x"]
    if len(x) == 0:
        return 0, 0, 0, 0, 0, 0

    mask = (x >= start_0) & (x < end_0)
    total_pos = int(np.sum(mask))
    if total_pos == 0:
        return 0, 0, 0, 0, 0, 0

    total_prev = round(float(np.sum(mm_dict["y"][mask])), 4)

    cats = mm_dict.get("codon_category")
    if cats is None:
        return total_pos, total_prev, 0, 0, 0, 0

    cat_arr = np.array(cats) if not isinstance(cats, np.ndarray) else cats
    syn_mask = mask & (cat_arr == "Synonymous")
    nonsyn_mask = mask & (cat_arr == "Non-synonymous")
    syn_pos = int(np.sum(syn_mask))
    syn_prev = round(float(np.sum(mm_dict["y"][syn_mask])), 4) if syn_pos else 0
    nonsyn_pos = int(np.sum(nonsyn_mask))
    nonsyn_prev = round(float(np.sum(mm_dict["y"][nonsyn_mask])), 4) if nonsyn_pos else 0

    return total_pos, total_prev, syn_pos, syn_prev, nonsyn_pos, nonsyn_prev


# ============================================================================
# Gene naming
# ============================================================================

def _build_gene_names(genes, contig_info):
    """Assign <contig_name>_tbb_<N> names to CDS features, numbered per contig by Start."""
    counter = defaultdict(int)
    names = []
    for g in genes:
        contig_id = g[0]
        contig_name = contig_info[contig_id][0]
        counter[contig_id] += 1
        names.append(f"{contig_name}_tbb_{counter[contig_id]}")
    return names


# ============================================================================
# Main processing
# ============================================================================

def process_sample(conn, sample_id, sample_name, contig_info, cov_map, af_map,
                   id_to_name, name_to_id, genes_by_contig, gene_names_by_contig,
                   mag_mode=False, contig_to_mag=None, mag_cov_map=None):
    """Process all CDS for one sample. Yields row tuples."""
    contig_ids = list(cov_map.keys())
    if not contig_ids:
        return

    for contig_id in contig_ids:
        contig_genes = genes_by_contig.get(contig_id)
        if not contig_genes:
            continue

        contig_name, contig_length, _ = contig_info[contig_id]
        contig_cov_tmean = cov_map[contig_id]
        contig_af = af_map[contig_id]

        mag_af = None
        mag_cov_tmean = None
        if mag_mode and contig_to_mag and mag_cov_map:
            mag_id = contig_to_mag.get(contig_id)
            if mag_id is not None:
                mag_af, mag_cov_tmean = mag_cov_map.get((mag_id, sample_id), (None, None))

        features = _load_all_features(
            conn, contig_id, sample_id, contig_length, id_to_name, name_to_id
        )

        gene_names = gene_names_by_contig.get(contig_id, [])

        for idx, g in enumerate(contig_genes):
            gene_start, gene_end = g[1], g[2]
            s_sites, n_sites = g[4], g[5]
            gene_name = gene_names[idx] if idx < len(gene_names) else ""

            s0 = gene_start - 1
            e0 = gene_end

            gene_slice = features["primary_reads"][s0:e0]
            gene_cov_med = float(np.median(gene_slice))
            gene_af = round(float(np.count_nonzero(gene_slice) / len(gene_slice) * 100), 2) if len(gene_slice) > 0 else 0
            gene_sec_med = float(np.median(features["secondary_reads"][s0:e0]))

            coverage_ratio = 0
            if contig_cov_tmean is not None and contig_cov_tmean > 0:
                coverage_ratio = round(gene_cov_med / contig_cov_tmean, 4)

            mm = _mismatch_stats(features["mismatches"], s0, e0)
            total_pos, total_prev = mm[0], mm[1]
            syn_pos, syn_prev = mm[2], mm[3]
            nonsyn_pos, nonsyn_prev = mm[4], mm[5]

            dnds = 0
            if (syn_prev and nonsyn_prev
                    and syn_prev > 0 and n_sites and s_sites):
                dnds = round((nonsyn_prev * s_sites) / (syn_prev * n_sites), 4)

            ins_pos, ins_prev = _sparse_stats(features["insertions"], s0, e0)
            del_pos, del_prev = _sparse_stats(features["deletions"], s0, e0)
            lc_pos, lc_prev = _sparse_stats(features["left_clippings"], s0, e0)
            rc_pos, rc_prev = _sparse_stats(features["right_clippings"], s0, e0)

            row = [
                sample_name,
                contig_name,
                gene_name,
            ]
            if mag_mode:
                row.append(mag_af)
            row.append(contig_af)
            row.append(gene_af)
            if mag_mode:
                row.append(mag_cov_tmean)
            row.extend([
                contig_cov_tmean,
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
            ])
            yield tuple(row)


DESCRIPTION = """\
Export per-CDS mapping signals from a theBIGbam database.

Output columns (one row per sample x CDS):
- sample_name: BAM sample name
- contig_name: contig the CDS belongs to
- gene_name: stable identifier <contig_name>_tbb_<N>, numbered per contig by start position
- mag_aligned_fraction: (MAG mode only) percentage of the MAG covered by aligned reads (0-100)
- contig_aligned_fraction: percentage of the contig covered by aligned reads (0-100)
- gene_aligned_fraction: percentage of positions in the CDS covered by at least one primary read (0-100)
- mag_coverage_trimmed_mean: (MAG mode only) trimmed mean primary read depth across the whole MAG
- contig_coverage_trimmed_mean: trimmed mean primary read depth across the whole contig
- gene_coverage_median: median primary read depth across the CDS
- coverage_ratio: gene_coverage_median / contig_coverage_trimmed_mean
- gene_secondary_coverage_median: median secondary (non-primary) read depth across the CDS
- mismatches_positions: number of positions in the CDS with at least one mismatch
- mismatches_prevalence: sum of per-position mismatch prevalences (count/coverage) across the CDS
- synonymous_positions: number of mismatch positions classified as synonymous
- synonymous_prevalence: sum of per-position prevalences for synonymous mismatches
- nonsynonymous_positions: number of mismatch positions classified as non-synonymous
- nonsynonymous_prevalence: sum of per-position prevalences for non-synonymous mismatches
- s_sites: number of synonymous sites in the CDS (from codon degeneracy)
- n_sites: number of non-synonymous sites in the CDS
- dnds_ratio: (nonsynonymous_prevalence x s_sites) / (synonymous_prevalence x n_sites)
- insertions_positions: number of positions in the CDS with at least one insertion
- insertions_prevalence: sum of per-position insertion prevalences across the CDS
- deletions_positions: number of positions in the CDS with at least one deletion
- deletions_prevalence: sum of per-position deletion prevalences across the CDS
- left_clippings_positions: number of positions with left soft-clipping events
- left_clippings_prevalence: sum of per-position left-clipping prevalences
- right_clippings_positions: number of positions with right soft-clipping events
- right_clippings_prevalence: sum of per-position right-clipping prevalences
"""


def add_args(parser):
    """Register arguments on an existing subparser."""
    parser.add_argument("--db", required=True, help="Path to the theBIGbam DuckDB database")
    parser.add_argument("--output", required=True, help="Path to the output TSV file")
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
        help="Minimum trimmed mean coverage depth to include a contig-sample pair (default: 10)",
    )


def run(args):
    """Entry point called by CLI dispatcher."""

    if not 0 <= args.min_aligned_fraction <= 100:
        print("Error: --min_aligned_fraction must be between 0 and 100.", file=sys.stderr)
        sys.exit(1)

    conn = duckdb.connect(args.db, read_only=True)

    from thebigbam.database.database_getters import is_mag_mode
    mag_mode = is_mag_mode(conn)
    contig_to_mag = {}
    mag_cov_map = {}
    if mag_mode:
        contig_to_mag = {
            r[0]: r[1] for r in conn.execute(
                "SELECT Contig_id, MAG_id FROM MAG_contigs_association"
            ).fetchall()
        }
        if _table_exists(conn, "MAG_coverage"):
            mag_cov_map = {
                (r[0], r[1]): (r[2] / 10.0, r[3] / 10.0)
                for r in conn.execute(
                    "SELECT MAG_id, Sample_id, Aligned_fraction_percentage, "
                    "Coverage_trimmed_mean FROM MAG_coverage"
                ).fetchall()
            }

    if not _table_exists(conn, "Contig_annotation"):
        print("Error: No Contig_annotation table. Run annotation first.", file=sys.stderr)
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

    contigs = conn.execute(
        "SELECT Contig_id, Contig_name, Contig_length, GC_mean FROM Contig ORDER BY Contig_id"
    ).fetchall()
    contig_info = {r[0]: (r[1], r[2], r[3]) for r in contigs}

    id_to_name, name_to_id = _resolve_feature_ids(conn)
    has_sn_sites = _column_exists(conn, "Contig_annotation", "S_sites")
    has_coverage = _table_exists(conn, "Coverage")

    # Fetch all CDS ordered by contig, start
    all_contig_ids = ",".join(str(c) for c in contig_info.keys())
    if has_sn_sites:
        all_genes = conn.execute(
            f"""
            SELECT ca.Contig_id, ca."Start", ca."End", ca."End" - ca."Start" + 1,
                   ca.S_sites, ca.N_sites
            FROM Contig_annotation ca
            WHERE ca.Contig_id IN ({all_contig_ids}) AND ca."Type" = 'CDS'
            ORDER BY ca.Contig_id, ca."Start"
            """
        ).fetchall()
    else:
        raw = conn.execute(
            f"""
            SELECT ca.Contig_id, ca."Start", ca."End", ca."End" - ca."Start" + 1,
                   ca.Nucleotide_sequence
            FROM Contig_annotation ca
            WHERE ca.Contig_id IN ({all_contig_ids}) AND ca."Type" = 'CDS'
            ORDER BY ca.Contig_id, ca."Start"
            """
        ).fetchall()
        all_genes = []
        for r in raw:
            nuc_seq = r[4]
            if nuc_seq:
                s, n = _compute_sn_sites(nuc_seq)
            else:
                s, n = None, None
            all_genes.append((r[0], r[1], r[2], r[3], s, n))

    # Build gene names and group by contig
    gene_names = _build_gene_names(all_genes, contig_info)
    genes_by_contig = defaultdict(list)
    gene_names_by_contig = defaultdict(list)
    for i, g in enumerate(all_genes):
        genes_by_contig[g[0]].append(g)
        gene_names_by_contig[g[0]].append(gene_names[i])

    print(
        f"Processing {len(samples)} sample(s), "
        f"{len(name_to_id)} feature(s), "
        f"{len(all_genes)} CDS...",
        file=sys.stderr,
    )
    print(
        f"Filters: aligned_fraction >= {args.min_aligned_fraction}%, "
        f"coverage_depth >= {args.min_coverage_depth}",
        file=sys.stderr,
    )

    af_threshold = args.min_aligned_fraction * 10
    cov_threshold = args.min_coverage_depth * 10

    total_rows = 0

    columns = _build_columns(mag_mode)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)

        for sample_id, sample_name in samples:
            if has_coverage:
                qualifying = conn.execute(
                    """
                    SELECT Contig_id, Aligned_fraction_percentage / 10.0,
                           Coverage_trimmed_mean / 10.0
                    FROM Coverage
                    WHERE Sample_id = ?
                      AND Aligned_fraction_percentage >= ?
                      AND Coverage_trimmed_mean >= ?
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
                id_to_name, name_to_id, genes_by_contig, gene_names_by_contig,
                mag_mode=mag_mode, contig_to_mag=contig_to_mag, mag_cov_map=mag_cov_map,
            ):
                writer.writerow(row)
                sample_rows += 1

            total_rows += sample_rows
            print(f"  Done: {sample_name} ({sample_rows} CDS)", file=sys.stderr)

    conn.close()
    print(f"Wrote {total_rows} rows to {args.output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_args(parser)
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
