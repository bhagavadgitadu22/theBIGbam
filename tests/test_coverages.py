"""
Coverage validation tests for theBIGbam pipeline.

Tests that coverage values in the database match expected values from CSV files.
All coverage features are now stored as compressed BLOBs in the Feature_blob table.
"""

import os
import pytest
import duckdb
import numpy as np


def get_contig_length(conn, contig_id):
    """Get contig length from database."""
    return conn.execute(
        "SELECT Contig_length FROM Contig WHERE Contig_id = ?", (contig_id,)
    ).fetchone()[0]


def blob_to_int_array(conn, contig_id, sample_id, feature_name, length):
    """Decode a Feature_blob entry to a per-position integer array.

    Returns a list of integers (length = contig length), 0-filled for positions
    without data. Dense blobs cover every position; sparse blobs only have events.
    """
    from thebigbam.database.blob_decoder import (
        feature_name_to_id, decode_raw_chunks, decode_raw_sparse_chunks,
        get_scale_from_zoom_blob, is_sparse_zoom_blob,
    )

    fid = feature_name_to_id(feature_name, conn)

    # Get scale/sparse info from zoom blob
    row = conn.execute(
        "SELECT Zoom_data FROM Feature_blob WHERE Contig_id=? AND Sample_id=? AND Feature_id=?",
        (contig_id, sample_id, fid)
    ).fetchone()
    if row is None:
        return [0] * length
    zoom_blob = bytes(row[0])
    scale_div = get_scale_from_zoom_blob(zoom_blob)
    sparse = is_sparse_zoom_blob(zoom_blob)

    # Decode from chunks
    chunk_rows = conn.execute(
        "SELECT Chunk_idx, Data FROM Feature_blob_chunk "
        "WHERE Contig_id=? AND Sample_id=? AND Feature_id=? ORDER BY Chunk_idx",
        (contig_id, sample_id, fid)
    ).fetchall()
    if not chunk_rows:
        return [0] * length
    chunk_rows = [(r[0], r[1]) for r in chunk_rows]

    data = decode_raw_sparse_chunks(chunk_rows, scale_div) if sparse else decode_raw_chunks(chunk_rows, scale_div)
    x = data["x"]
    y = data["y"]

    arr = [0] * length
    for i in range(len(x)):
        pos = int(x[i])
        if 0 <= pos < length:
            arr[pos] = int(round(y[i]))
    return arr


def fold_coverage(arr, target_length):
    """Fold a doubled-contig array back to original length by summing halves."""
    if len(arr) <= target_length:
        return arr
    folded = arr[:target_length]
    for i in range(target_length, len(arr)):
        folded[i - target_length] += arr[i]
    return folded


def get_coverage_arrays(conn, contig_id, sample_id):
    """Get primary, secondary, supplementary coverage arrays from Feature_blob."""
    length = get_contig_length(conn, contig_id)
    primary = blob_to_int_array(conn, contig_id, sample_id, "primary_reads", length)
    secondary = blob_to_int_array(conn, contig_id, sample_id, "secondary_reads", length)
    supplementary = blob_to_int_array(conn, contig_id, sample_id, "supplementary_reads", length)
    return primary, secondary, supplementary


def get_coverage_reduced(conn, contig_id, sample_id):
    """Get coverage_reduced array from Feature_blob."""
    length = get_contig_length(conn, contig_id)
    return blob_to_int_array(conn, contig_id, sample_id, "coverage_reduced", length)


@pytest.mark.parametrize("sample_name,coverage_file", [
    ("50_read_pairs", "coverage_50_read_pairs.csv"),
    ("50_read_pairs_inverted", "coverage_50_read_pairs_inverted.csv"),
    ("5000_read_pairs_concat", "coverage_5000_read_pairs_concat.csv"),
    ("1000_long_reads", "coverage_1000_long_reads.csv"),
    ("100_long_reads_concat", "coverage_100_long_reads_concat.csv"),
])
def test_coverage_linear_database(linear_db, tests_dir, sample_name, coverage_file):
    """Test that file_coverage = primary + secondary + supplementary for linear DB."""
    coverage_path = os.path.join(tests_dir, "coverages", coverage_file)

    with open(coverage_path) as f:
        expected = [int(line.strip()) for line in f]

    conn = duckdb.connect(linear_db, read_only=True)
    cur = conn.cursor()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    primary, secondary, supplementary = get_coverage_arrays(conn, contig_id, sample_id)
    conn.close()

    computed = [p + s + sup for p, s, sup in zip(primary, secondary, supplementary)]
    assert computed == expected, f"Coverage mismatch for {sample_name}"


@pytest.mark.parametrize("sample_name,coverage_file", [
    ("50_read_pairs_circular", "coverage_50_read_pairs_circular.csv"),
    ("50_read_pairs_inverted_circular", "coverage_50_read_pairs_inverted_circular.csv"),
    ("1000_long_reads_circular", "coverage_1000_long_reads_circular.csv"),
    ("100_long_reads_concat_circular", "coverage_100_long_reads_concat_circular.csv"),
    ("5000_read_pairs_concat_circular", "coverage_5000_read_pairs_concat_circular.csv"),
])
def test_coverage_circular_database(circular_db, tests_dir, sample_name, coverage_file):
    """Test that file_coverage = primary + secondary + supplementary for circular DB."""
    coverage_path = os.path.join(tests_dir, "coverages", coverage_file)

    with open(coverage_path) as f:
        expected = [int(line.strip()) for line in f]

    conn = duckdb.connect(circular_db, read_only=True)
    cur = conn.cursor()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    primary, secondary, supplementary = get_coverage_arrays(conn, contig_id, sample_id)
    conn.close()

    computed = [p + s + sup for p, s, sup in zip(primary, secondary, supplementary)]
    computed = fold_coverage(computed, len(expected))
    assert computed == expected, f"Coverage mismatch for {sample_name}"


@pytest.mark.parametrize("linear_sample,circular_sample", [
    ("50_read_pairs", "50_read_pairs_circular"),
    ("50_read_pairs_inverted", "50_read_pairs_inverted_circular"),
    ("1000_long_reads", "1000_long_reads_circular"),
])
def test_validate_circular_coverages(linear_db, circular_db, linear_sample, circular_sample):
    """Test that coverage is the same between linear and circular mapping modes."""
    conn_linear = duckdb.connect(linear_db, read_only=True)
    cur = conn_linear.cursor()
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (linear_sample,))
    linear_sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    linear_contig_id = cur.fetchone()[0]
    linear_primary, linear_secondary, linear_supplementary = get_coverage_arrays(
        conn_linear, linear_contig_id, linear_sample_id
    )
    conn_linear.close()
    linear_total = [p + s + sup for p, s, sup in zip(linear_primary, linear_secondary, linear_supplementary)]
    linear_total = fold_coverage(linear_total, 10000)

    conn_circular = duckdb.connect(circular_db, read_only=True)
    cur = conn_circular.cursor()
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (circular_sample,))
    circular_sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    circular_contig_id = cur.fetchone()[0]
    circular_primary, circular_secondary, circular_supplementary = get_coverage_arrays(
        conn_circular, circular_contig_id, circular_sample_id
    )
    conn_circular.close()
    circular_total = [p + s + sup for p, s, sup in zip(circular_primary, circular_secondary, circular_supplementary)]
    circular_total = fold_coverage(circular_total, 10000)

    assert linear_total == circular_total, f"Coverage mismatch between {linear_sample} and {circular_sample}"


@pytest.mark.parametrize("db_type,sample_name", [
    ("linear", "50_read_pairs"),
    ("circular", "50_read_pairs_circular"),
    ("linear", "50_read_pairs_inverted"),
    ("circular", "50_read_pairs_inverted_circular"),
    ("linear", "1000_long_reads"),
    ("circular", "1000_long_reads_circular"),
])
def test_validate_coverage_reduced(linear_db, circular_db, db_type, sample_name):
    """Test that primary_reads equals coverage_reduced for simple samples."""
    db_path = linear_db if db_type == "linear" else circular_db

    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    primary, _, _ = get_coverage_arrays(conn, contig_id, sample_id)
    coverage_reduced = get_coverage_reduced(conn, contig_id, sample_id)

    conn.close()

    assert primary == coverage_reduced, f"primary_reads != coverage_reduced for {sample_name}"


@pytest.mark.parametrize("db_type,sample_name", [
    ("linear", "1000_long_reads"),
    ("circular", "1000_long_reads_circular"),
    ("circular", "5000_read_pairs_concat_circular"),
    ("circular", "100_long_reads_concat_circular"),
])
def test_validate_simple_coverage(linear_db_simple, circular_db_simple, db_type, sample_name):
    """Test that primary_reads is uniform (single value) without percentage options."""
    db_path = linear_db_simple if db_type == "linear" else circular_db_simple

    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]
    length = get_contig_length(conn, contig_id)

    primary = blob_to_int_array(conn, contig_id, sample_id, "primary_reads", length)
    conn.close()

    unique_values = set(primary)
    assert len(unique_values) == 1, f"Expected uniform coverage for {sample_name}, got {len(unique_values)} distinct values"
