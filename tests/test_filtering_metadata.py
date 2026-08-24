import duckdb

from thebigbam.database.database_getters import get_filtering_metadata


def test_filtering_metadata_excludes_relation_identifiers(tmp_path):
    db_path = tmp_path / "metadata.duckdb"
    connection = duckdb.connect(str(db_path))
    connection.execute("CREATE TABLE Database_metadata (Key VARCHAR, Value VARCHAR)")
    connection.execute("INSERT INTO Database_metadata VALUES ('View_mode', 'mag')")
    connection.execute(
        "CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR, Contig_length INTEGER)"
    )
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR, Group_name VARCHAR)")
    connection.execute("CREATE TABLE MAG (MAG_id INTEGER, MAG_name VARCHAR, Completeness DOUBLE)")
    connection.execute(
        "CREATE TABLE Explicit_coverage ("
        "Contig_id INTEGER, Contig_name VARCHAR, Sample_id INTEGER, Sample_name VARCHAR, "
        "Mean_coverage DOUBLE, Coverage_class VARCHAR)"
    )
    connection.execute(
        "CREATE TABLE Explicit_coverage_per_MAG ("
        "MAG_id INTEGER, MAG_name VARCHAR, Sample_id INTEGER, Sample_name VARCHAR, "
        "Mean_coverage DOUBLE, Coverage_class VARCHAR)"
    )
    connection.close()

    metadata = get_filtering_metadata(str(db_path))

    assert set(metadata["Coverage"]["columns"]) == {"Mean_coverage", "Coverage_class"}
    assert set(metadata["MAG coverage"]["columns"]) == {"Mean_coverage", "Coverage_class"}
