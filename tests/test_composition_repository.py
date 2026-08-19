from types import SimpleNamespace

import duckdb
import pytest

from thebigbam.plotting.models.composition import MagCompositionRequest, MagOrdering
from thebigbam.plotting.repositories.composition import CompositionRepository, MagContext
from thebigbam.plotting.services.composition import CompositionDataService


def test_repository_validates_contigs_and_optional_samples():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR, Contig_length INTEGER)")
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    connection.execute("INSERT INTO Contig VALUES (2, 'c2', 500)")
    connection.execute("INSERT INTO Sample VALUES (3, 's3')")
    repository = CompositionRepository(connection)

    assert repository.contig("c2").length == 500
    assert repository.sample("s3") == (3, "s3")
    assert repository.sample(None) == (None, None)
    with pytest.raises(ValueError, match="Contig not found"):
        repository.contig("missing")
    with pytest.raises(ValueError, match="Sample not found"):
        repository.sample("missing")
    connection.close()


def test_repository_filters_orders_and_limits_mag_samples():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    connection.execute("CREATE TABLE Feature_blob (Sample_id INTEGER, Contig_id INTEGER)")
    connection.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    connection.execute("INSERT INTO Sample VALUES (1, 'zeta'), (2, 'alpha'), (3, 'beta')")
    connection.execute("INSERT INTO Feature_blob VALUES (1, 10), (2, 10), (3, 10)")
    connection.execute("INSERT INTO MAG_contigs_association VALUES (7, 10)")
    request = SimpleNamespace(
        allowed_samples=frozenset({"zeta", "alpha"}),
        sample_ordering=SimpleNamespace(source=None, metric=None, ascending=True),
        max_samples=1,
        mag_name="m7",
    )

    repository = CompositionRepository(connection)
    rows = CompositionDataService(object(), repository).ordered_mag_samples(request, 7)

    assert rows == [(2, "alpha")]
    connection.close()


def test_mag_context_projects_named_genome_members_without_database_lookup():
    context = MagContext(
        mag_id=7,
        members=(("c2", 20, 10), ("c1", 10, 0)),
        feature_members=((22, 20, 10), (11, 10, 0)),
        total_length=30,
        xstart=1,
        xend=30,
    )

    assert context.genome_members == ((22, "c2", 20, 10), (11, "c1", 10, 0))


def test_sorted_mag_context_reuses_one_full_member_query(monkeypatch):
    calls = []

    monkeypatch.setattr("thebigbam.plotting.repositories.composition.get_mag_id", lambda _conn, _name: 7)

    def load_sorted(*args):
        calls.append(args)
        return [(11, "c1", 10, 0), (22, "c2", 20, 10)]

    monkeypatch.setattr("thebigbam.plotting.repositories.composition.get_mag_members_full_sorted", load_sorted)
    request = MagCompositionRequest(
        mag_name="m7",
        sample_name="s1",
        ordering=MagOrdering(source="Coverage", metric="Depth", ascending=False, sample_name="s1"),
    )

    context = CompositionRepository(object()).mag(request)

    assert len(calls) == 1
    assert context.members == (("c1", 10, 0), ("c2", 20, 10))
    assert context.feature_members == ((11, 10, 0), (22, 20, 10))
