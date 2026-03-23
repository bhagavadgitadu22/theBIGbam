"""
Tests for BAM generation and database creation verification.
"""

import os
import duckdb


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

    # Check that contigs were processed
    contig_count = conn.execute("SELECT COUNT(*) FROM Contig").fetchone()[0]
    assert contig_count > 0, "No contigs in database"

    # Check that variables were created
    var_count = conn.execute("SELECT COUNT(*) FROM Variable").fetchone()[0]
    assert var_count > 0, "No variables in database"

    # Check that Feature_blob table exists and has data
    tables = conn.execute("SHOW TABLES").fetchall()
    table_names = [t[0] for t in tables]
    assert "Feature_blob" in table_names, f"Missing 'Feature_blob' table. Found: {table_names}"

    blob_count = conn.execute("SELECT COUNT(*) FROM Feature_blob").fetchone()[0]
    assert blob_count > 0, "No feature blobs in database"

    conn.close()


def test_bam_generation(test_bams):
    """Verify BAMs were created for both linear and circular modes."""
    assert len(test_bams["linear"]) > 0, "No linear BAM files generated"
    assert len(test_bams["circular"]) > 0, "No circular BAM files generated"


def test_linear_database(linear_db):
    """Verify linear database structure."""
    verify_database(linear_db)


def test_circular_database(circular_db):
    """Verify circular database structure."""
    verify_database(circular_db)


def test_linear_database_simple(linear_db_simple):
    """Verify simple linear database structure."""
    verify_database(linear_db_simple)


def test_circular_database_simple(circular_db_simple):
    """Verify simple circular database structure."""
    verify_database(circular_db_simple)
