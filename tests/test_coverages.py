"""
Integration tests for MGFeatureViewer pipeline.

Tests the calculate command with generated BAM files.
"""

import os
import subprocess
import pytest
import duckdb

# Output database paths
TEST_DB_LINEAR = os.path.join(os.path.dirname(__file__), "test_10kbp_linear.db")
TEST_DB_CIRCULAR = os.path.join(os.path.dirname(__file__), "test_10kbp_circular.db")

def cleanup_dbs():
    """Remove test databases before test run."""
    for db in [TEST_DB_LINEAR, TEST_DB_CIRCULAR]:
        if os.path.exists(db):
            os.remove(db)

def verify_database(db_path: str):
    """Verify database structure and content."""
    assert os.path.exists(db_path), f"Database was not created: {db_path}"

    conn = duckdb.connect(db_path, read_only=True)

    # Check that main tables exist
    tables = conn.execute("SHOW TABLES").fetchall()
    table_names = [t[0] for t in tables]

    assert "Contig" in table_names, "Missing 'Contig' table"
    assert "Sample" in table_names, "Missing 'Sample' table"
    assert "Variable" in table_names, "Missing 'Variable' table"

    # Check that samples were processed
    sample_count = conn.execute("SELECT COUNT(*) FROM Sample").fetchone()[0]
    assert sample_count > 0, "No samples in database"
    print(f"  Processed {sample_count} samples")

    # Check that contigs were processed
    contig_count = conn.execute("SELECT COUNT(*) FROM Contig").fetchone()[0]
    assert contig_count > 0, "No contigs in database"
    print(f"  Found {contig_count} contigs")

    # Check that variables were created
    var_count = conn.execute("SELECT COUNT(*) FROM Variable").fetchone()[0]
    assert var_count > 0, "No variables in database"
    print(f"  Created {var_count} variables")

    # List variables
    variables = conn.execute("SELECT Feature_table_name FROM Variable").fetchall()
    print(f"  Variables: {[v[0] for v in variables]}")

    # Check for expected coverage variable
    var_names = [v[0] for v in variables]
    assert "Feature_primary_reads" in var_names, f"Missing primary_reads variable. Found: {var_names}"

    conn.close()

def run_calculate(bam_dir: str, output_db: str, circular: bool):
    """Run calculate command on BAM files."""
    cmd = [
        "mgfeatureviewer", "calculate",
        "-b", bam_dir,
        "-m", "coverage,phagetermini,assemblycheck",
        "-o", output_db,
        "-t", "4",
        "--coverage_percentage", "0",
        "--variation_percentage", "0",
    ]
    if circular:
        cmd.append("--circular")

    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    return result

def test_calculate_linear_bams(test_bams, tests_dir):
    """
    Test running calculate on linear (non-circular) BAM files.
    """
    # Clean up before test
    if os.path.exists(TEST_DB_LINEAR):
        os.remove(TEST_DB_LINEAR)

    linear_bams = test_bams["linear"]
    assert len(linear_bams) > 0, "No linear BAM files available"
    print(f"\nTesting with {len(linear_bams)} linear BAMs")

    linear_dir = os.path.join(tests_dir, "linear_bams")
    result = run_calculate(linear_dir, TEST_DB_LINEAR, circular=False)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"

    print(f"\nVerifying linear database: {TEST_DB_LINEAR}")
    verify_database(TEST_DB_LINEAR)

def test_calculate_circular_bams(test_bams, tests_dir):
    """
    Test running calculate on circular BAM files.
    """
    # Clean up before test
    if os.path.exists(TEST_DB_CIRCULAR):
        os.remove(TEST_DB_CIRCULAR)

    circular_bams = test_bams["circular"]
    assert len(circular_bams) > 0, "No circular BAM files available"
    print(f"\nTesting with {len(circular_bams)} circular BAMs")

    circular_dir = os.path.join(tests_dir, "circular_bams")
    result = run_calculate(circular_dir, TEST_DB_CIRCULAR, circular=True)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"

    print(f"\nVerifying circular database: {TEST_DB_CIRCULAR}")
    verify_database(TEST_DB_CIRCULAR)


# =============================================================================
# Coverage Validation Tests
# =============================================================================

def expand_rle_to_array(rows, length=10000):
    """Convert RLE rows (first_pos, last_pos, value) to per-position array.

    Database positions are 1-indexed (1 to 10000), convert to 0-indexed (0 to 9999).
    """
    arr = [0] * length
    for first_pos, last_pos, value in rows:
        # Convert from 1-indexed (database) to 0-indexed (Python)
        # Use slice assignment for efficiency
        start = first_pos - 1
        end = last_pos
        arr[start:end] = [value] * (end - start)
    return arr


def get_coverage_arrays(conn, contig_id, sample_id):
    """Get primary, secondary, supplementary coverage arrays from database."""
    cur = conn.cursor()

    # Primary reads
    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_primary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    primary = expand_rle_to_array(cur.fetchall())

    # Secondary reads
    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_secondary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    secondary = expand_rle_to_array(cur.fetchall())

    # Supplementary reads
    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_supplementary_reads
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    supplementary = expand_rle_to_array(cur.fetchall())

    return primary, secondary, supplementary


@pytest.mark.parametrize("sample_name,coverage_file", [
    ("50_read_pairs", "coverage_50_read_pairs.csv"),
    ("50_read_pairs_inverted", "coverage_50_read_pairs_inverted.csv"),
    ("5000_read_pairs_concat", "coverage_5000_read_pairs_concat.csv"),
    ("1000_long_reads", "coverage_1000_long_reads.csv"),
    ("100_long_reads_concat", "coverage_100_long_reads_concat.csv"),
])
def test_coverage_linear_database(tests_dir, sample_name, coverage_file):
    """Test that file_coverage = primary + secondary + supplementary for linear DB."""
    db_path = os.path.join(tests_dir, "test_10kbp_linear.db")
    coverage_path = os.path.join(tests_dir, "coverages", coverage_file)

    # Load expected coverage from file
    with open(coverage_path) as f:
        expected = [int(line.strip()) for line in f]

    # Query database
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # Get IDs
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    primary, secondary, supplementary = get_coverage_arrays(conn, contig_id, sample_id)
    conn.close()

    # Compare: file_coverage == primary + secondary + supplementary
    computed = [p + s + sup for p, s, sup in zip(primary, secondary, supplementary)]
    assert computed == expected, f"Coverage mismatch for {sample_name}"


@pytest.mark.parametrize("sample_name,coverage_file", [
    ("50_read_pairs_circular", "coverage_50_read_pairs_circular.csv"),
    ("50_read_pairs_inverted_circular", "coverage_50_read_pairs_inverted_circular.csv"),
    ("1000_long_reads_circular", "coverage_1000_long_reads_circular.csv"),
    ("100_long_reads_concat_circular", "coverage_100_long_reads_concat_circular.csv"),
    ("5000_read_pairs_concat_circular", "coverage_5000_read_pairs_concat_circular.csv"),
])
def test_coverage_circular_database(tests_dir, sample_name, coverage_file):
    """Test that file_coverage = primary + secondary + supplementary for circular DB."""
    db_path = os.path.join(tests_dir, "test_10kbp_circular.db")
    coverage_path = os.path.join(tests_dir, "coverages", coverage_file)

    # Load expected coverage from file
    with open(coverage_path) as f:
        expected = [int(line.strip()) for line in f]

    # Query database
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # Get IDs
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    primary, secondary, supplementary = get_coverage_arrays(conn, contig_id, sample_id)
    conn.close()

    # Compare: file_coverage == primary + secondary + supplementary
    computed = [p + s + sup for p, s, sup in zip(primary, secondary, supplementary)]
    assert computed == expected, f"Coverage mismatch for {sample_name}"


def get_coverage_reduced(conn, contig_id, sample_id):
    """Get coverage_reduced array from database."""
    cur = conn.cursor()
    cur.execute("""
        SELECT First_position, Last_position, Value
        FROM Feature_coverage_reduced
        WHERE Contig_id = ? AND Sample_id = ?
        ORDER BY First_position
    """, (contig_id, sample_id))
    return expand_rle_to_array(cur.fetchall())


@pytest.mark.parametrize("linear_sample,circular_sample", [
    ("50_read_pairs", "50_read_pairs_circular"),
    ("50_read_pairs_inverted", "50_read_pairs_inverted_circular"),
    ("1000_long_reads", "1000_long_reads_circular"),
])
def test_validate_circular_coverages(tests_dir, linear_sample, circular_sample):
    """Test that coverage is the same between linear and circular mapping modes."""
    linear_db = os.path.join(tests_dir, "test_10kbp_linear.db")
    circular_db = os.path.join(tests_dir, "test_10kbp_circular.db")

    # Get linear coverage
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

    # Get circular coverage
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

    # Compare
    assert linear_total == circular_total, f"Coverage mismatch between {linear_sample} and {circular_sample}"


@pytest.mark.parametrize("db_type,sample_name", [
    ("linear", "50_read_pairs"),
    ("circular", "50_read_pairs_circular"),
    ("linear", "50_read_pairs_inverted"),
    ("circular", "50_read_pairs_inverted_circular"),
    ("linear", "1000_long_reads"),
    ("circular", "1000_long_reads_circular"),
])
def test_validate_coverage_reduced(tests_dir, db_type, sample_name):
    """Test that primary_reads equals coverage_reduced for simple samples."""
    if db_type == "linear":
        db_path = os.path.join(tests_dir, "test_10kbp_linear.db")
    else:
        db_path = os.path.join(tests_dir, "test_10kbp_circular.db")

    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # Get IDs
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    # Get primary reads
    primary, _, _ = get_coverage_arrays(conn, contig_id, sample_id)

    # Get coverage_reduced
    coverage_reduced = get_coverage_reduced(conn, contig_id, sample_id)

    conn.close()

    # Compare
    assert primary == coverage_reduced, f"primary_reads != coverage_reduced for {sample_name}"


# =============================================================================
# Simple RLE Validation Tests (without coverage_percentage/variation_percentage)
# =============================================================================

TEST_DB_LINEAR_SIMPLE = os.path.join(os.path.dirname(__file__), "test_10kbp_linear_simple.db")
TEST_DB_CIRCULAR_SIMPLE = os.path.join(os.path.dirname(__file__), "test_10kbp_circular_simple.db")


def run_calculate_simple(bam_dir: str, output_db: str, circular: bool):
    """Run calculate command without coverage_percentage/variation_percentage options."""
    cmd = [
        "mgfeatureviewer", "calculate",
        "-b", bam_dir,
        "-m", "coverage,phagetermini,assemblycheck",
        "-o", output_db,
        "-t", "4",
    ]
    if circular:
        cmd.append("--circular")

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


def test_calculate_linear_simple(test_bams, tests_dir):
    """Create linear database without coverage_percentage/variation_percentage options."""
    if os.path.exists(TEST_DB_LINEAR_SIMPLE):
        os.remove(TEST_DB_LINEAR_SIMPLE)

    linear_dir = os.path.join(tests_dir, "linear_bams")
    result = run_calculate_simple(linear_dir, TEST_DB_LINEAR_SIMPLE, circular=False)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"


def test_calculate_circular_simple(test_bams, tests_dir):
    """Create circular database without coverage_percentage/variation_percentage options."""
    if os.path.exists(TEST_DB_CIRCULAR_SIMPLE):
        os.remove(TEST_DB_CIRCULAR_SIMPLE)

    circular_dir = os.path.join(tests_dir, "circular_bams")
    result = run_calculate_simple(circular_dir, TEST_DB_CIRCULAR_SIMPLE, circular=True)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"


@pytest.mark.parametrize("db_type,sample_name", [
    ("linear", "1000_long_reads"),
    ("circular", "1000_long_reads_circular"),
    ("circular", "5000_read_pairs_concat_circular"),
    ("circular", "100_long_reads_concat_circular"),
])
def test_validate_simple_RLE(tests_dir, db_type, sample_name):
    """Test that primary_reads contains only one RLE row (uniform coverage) without percentage options."""
    if db_type == "linear":
        db_path = TEST_DB_LINEAR_SIMPLE
    else:
        db_path = TEST_DB_CIRCULAR_SIMPLE

    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # Get IDs
    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name = ?", (sample_name,))
    sample_id = cur.fetchone()[0]
    cur.execute("SELECT Contig_id FROM Contig LIMIT 1")
    contig_id = cur.fetchone()[0]

    # Count RLE rows for primary_reads
    cur.execute("""
        SELECT COUNT(*)
        FROM Feature_primary_reads
        WHERE Contig_id = ? AND Sample_id = ?
    """, (contig_id, sample_id))
    row_count = cur.fetchone()[0]

    conn.close()

    print(f"Sample {sample_name} has {row_count} RLE rows in primary_reads")

    assert row_count == 1, f"Expected 1 RLE row for {sample_name}, got {row_count}"
