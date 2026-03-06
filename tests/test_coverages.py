"""
Coverage validation tests for theBIGbam pipeline.

Tests that coverage values in the database match expected values from CSV files.
"""

import os
import pytest
import duckdb


def expand_rle_to_array(rows, length):
    """Convert RLE rows (first_pos, last_pos, value) to per-position array.

    Database positions are 1-indexed (1 to length), convert to 0-indexed.
    """
    arr = [0] * length
    for first_pos, last_pos, value in rows:
        start = first_pos - 1
        end = last_pos
        arr[start:end] = [value] * (end - start)
    return arr


def get_contig_length(conn, contig_id):
    """Get contig length from database."""
    return conn.execute(
        "SELECT Contig_length FROM Contig WHERE Contig_id = ?", (contig_id,)
    ).fetchone()[0]


def fold_coverage(arr, target_length):
    """Fold a doubled-contig array back to original length by summing halves."""
    if len(arr) <= target_length:
        return arr
    folded = arr[:target_length]
    for i in range(target_length, len(arr)):
        folded[i - target_length] += arr[i]
    return folded


def get_coverage_arrays(conn, contig_id, sample_id):
    """Get primary, secondary, supplementary coverage arrays from database."""
    length = get_contig_length(conn, contig_id)
    cur = conn.cursor()

    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_primary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    primary = expand_rle_to_array(cur.fetchall(), length)

    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_secondary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    secondary = expand_rle_to_array(cur.fetchall(), length)

    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_supplementary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    supplementary = expand_rle_to_array(cur.fetchall(), length)

    return primary, secondary, supplementary


def get_coverage_reduced(conn, contig_id, sample_id):
    """Get coverage_reduced array from database."""
    length = get_contig_length(conn, contig_id)
    cur = conn.cursor()
    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_coverage_reduced
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    return expand_rle_to_array(cur.fetchall(), length)


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
def test_validate_simple_RLE(linear_db_simple, circular_db_simple, db_type, sample_name):
    """Test that primary_reads contains only one RLE row (uniform coverage) without percentage options."""
    db_path = linear_db_simple if db_type == "linear" else circular_db_simple

    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    cur.execute("""
        SELECT COUNT(*)
        FROM Feature_primary_reads
        WHERE Contig_id = ? AND Sample_id = ?
    """, (contig_id, sample_id))
    row_count = cur.fetchone()[0]

    conn.close()

    assert row_count == 1, f"Expected 1 RLE row for {sample_name}, got {row_count}"
