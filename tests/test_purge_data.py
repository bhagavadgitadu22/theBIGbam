import duckdb
import pytest

from thebigbam.database.purge_data import purge_mapping_data


@pytest.mark.parametrize("mode", ["contig", "mag", "legacy"])
def test_purge_mapping_data_handles_sample_count_schemas(tmp_path, monkeypatch, mode):
    monkeypatch.setattr("importlib.metadata.version", lambda _name: "test-version")
    source = tmp_path / "source.db"
    output = tmp_path / "purged.db"
    with duckdb.connect(str(source)) as conn:
        conn.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR)")
        conn.execute("INSERT INTO Contig VALUES (1, 'contig-1')")
        if mode == "contig":
            conn.execute("ALTER TABLE Contig ADD COLUMN Number_of_samples INTEGER DEFAULT 1")
        if mode == "mag":
            conn.execute("CREATE TABLE MAG (MAG_id INTEGER, MAG_name VARCHAR, Number_of_samples INTEGER)")
            conn.execute("INSERT INTO MAG VALUES (1, 'bin-1', 1)")
            conn.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
            conn.execute("INSERT INTO MAG_contigs_association VALUES (1, 1)")
            conn.execute("CREATE TABLE MAG_coverage (MAG_id INTEGER, Sample_id INTEGER)")
            conn.execute("INSERT INTO MAG_coverage VALUES (1, 1)")
        conn.execute("CREATE TABLE Sample (Sample_id INTEGER)")
        conn.execute("INSERT INTO Sample VALUES (1)")
        conn.execute("CREATE TABLE Coverage (Contig_id INTEGER, Sample_id INTEGER)")
        conn.execute("INSERT INTO Coverage VALUES (1, 1)")
        conn.execute("CREATE TABLE Database_metadata (Key VARCHAR, Value VARCHAR)")
        conn.execute("INSERT INTO Database_metadata VALUES ('View_mode', ?), ('Min_coverage_depth', '1')", [mode])

    purge_mapping_data(str(source), str(output))

    with duckdb.connect(str(output), read_only=True) as conn:
        assert conn.execute("SELECT Contig_id, Contig_name FROM Contig").fetchall() == [(1, "contig-1")]
        for table in ("Sample", "Coverage"):
            assert conn.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone() == (0,)
        if mode in ("contig", "mag"):
            table = "MAG" if mode == "mag" else "Contig"
            assert conn.execute(f'SELECT Number_of_samples FROM "{table}"').fetchall() == [(None,)]
        if mode == "mag":
            assert conn.execute("SELECT MAG_id, MAG_name FROM MAG").fetchall() == [(1, "bin-1")]
            assert conn.execute("SELECT * FROM MAG_contigs_association").fetchall() == [(1, 1)]
            assert conn.execute("SELECT COUNT(*) FROM MAG_coverage").fetchone() == (0,)
        assert conn.execute("SELECT Value FROM Database_metadata WHERE Key = 'View_mode'").fetchone() == (mode,)
        assert conn.execute("SELECT Value FROM Database_metadata WHERE Key = 'Min_coverage_depth'").fetchone() is None

    with duckdb.connect(str(source), read_only=True) as conn:
        assert conn.execute("SELECT COUNT(*) FROM Coverage").fetchone() == (1,)
        if mode in ("contig", "mag"):
            table = "MAG" if mode == "mag" else "Contig"
            assert conn.execute(f'SELECT Number_of_samples FROM "{table}"').fetchone() == (1,)
