import duckdb

from thebigbam.plotting.repositories.filtering import FilteringRepository


def test_unbounded_filtered_samples_respects_contig_and_mag_scope():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    connection.execute("CREATE TABLE MAG (MAG_id INTEGER, MAG_name VARCHAR)")
    connection.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    connection.execute("INSERT INTO Sample VALUES (1, 'only_match'), (2, 'other_mag'), (3, 'other_contig')")
    connection.execute("INSERT INTO MAG VALUES (10, 'selected'), (20, 'other')")
    connection.execute("INSERT INTO MAG_contigs_association VALUES (10, 100), (20, 100), (10, 200)")

    repository = FilteringRepository(connection)
    sql = "SELECT * FROM (VALUES (100, 1), (100, 2), (200, 3)) AS f(Contig_id, Sample_id)"
    assert repository.filtered_samples(sql, [], limit=None) == ("only_match", "other_contig", "other_mag")

    samples = repository.filtered_samples(
        sql,
        [],
        contig_id=100,
        mag_name="selected",
        limit=None,
    )

    assert samples == ("only_match", "other_mag")


def test_unbounded_filtered_samples_is_not_limited_to_autocomplete_payload():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Sample AS SELECT i AS Sample_id, 'sample_' || i AS Sample_name FROM range(150) t(i)")
    repository = FilteringRepository(connection)
    sql = "SELECT 1 AS Contig_id, i AS Sample_id FROM range(150) t(i)"

    assert len(repository.filtered_samples(sql, [], contig_id=1)) == 100
    assert len(repository.filtered_samples(sql, [], contig_id=1, limit=None)) == 150
