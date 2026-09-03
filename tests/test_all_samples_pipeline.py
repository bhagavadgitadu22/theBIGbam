import sys
import types

import duckdb
import numpy as np
import pytest

# Unit tests mock the compiled Rust chunk primitives; integration environments
# continue to use the real extension through blob_decoder.
if "thebigbam_rs" not in sys.modules:
    rust_stub = types.ModuleType("thebigbam_rs")
    rust_stub.decode_dense_chunk = lambda data: []
    rust_stub.decode_sparse_chunk = lambda data: ([], [])
    sys.modules["thebigbam_rs"] = rust_stub

from thebigbam.plotting.models.plots import (
    AllSamplesPlotData,
    AllSamplesPlotRequest,
    DisplayOptions,
    FeatureSeries,
    FeatureStyle,
    GenomicRegion,
    SampleOrdering,
    SampleTrack,
)
from thebigbam.plotting.renderers.all_samples import AllSamplesRenderer
from thebigbam.plotting.repositories.all_samples import AllSamplesRepository
from thebigbam.plotting.services.all_samples import AllSamplesDataService


def test_genomic_region_allows_display_window_before_first_base():
    region = GenomicRegion(-100, 100)

    assert region.start == -100
    assert region.length == 200


@pytest.fixture
def repository_db(tmp_path):
    path = tmp_path / "repository.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR, Contig_length INTEGER)")
    conn.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR, rank INTEGER)")
    conn.execute("CREATE TABLE Coverage (Contig_id INTEGER, Sample_id INTEGER)")
    conn.execute("CREATE TABLE Explicit_coverage (Contig_id INTEGER, Sample_id INTEGER, Coverage_class VARCHAR)")
    conn.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    conn.execute("CREATE TABLE Explicit_coverage_per_MAG (MAG_id INTEGER, Sample_id INTEGER, Coverage_mean REAL)")
    conn.execute(
        'CREATE TABLE Variable (Variable_id INTEGER, Variable_name VARCHAR, Subplot VARCHAR, "Type" VARCHAR, '
        'Encoding VARCHAR, Color VARCHAR, Alpha REAL, Fill_alpha REAL, "Size" REAL, Title VARCHAR, '
        "Feature_table_name VARCHAR)"
    )
    conn.execute(
        "CREATE TABLE Feature_blob_chunk (Contig_id INTEGER, Sample_id INTEGER, Feature_id INTEGER, "
        "Chunk_idx INTEGER, Data BLOB)"
    )
    conn.execute("INSERT INTO Contig VALUES (1, 'c1', 100000)")
    conn.executemany(
        "INSERT INTO Sample VALUES (?, ?, ?)",
        [(index, f"sample_{index:03d}", 151 - index) for index in range(1, 151)],
    )
    conn.executemany("INSERT INTO Coverage VALUES (1, ?)", [(index,) for index in range(1, 151)])
    conn.executemany(
        "INSERT INTO Explicit_coverage VALUES (1, ?, ?)",
        [(index, "high" if index % 2 else "low") for index in range(1, 151)],
    )
    conn.execute("INSERT INTO MAG_contigs_association VALUES (10, 1)")
    conn.executemany(
        "INSERT INTO Explicit_coverage_per_MAG VALUES (10, ?, ?)",
        [(1, 100.0), (2, 200.0)],
    )
    conn.executemany(
        "INSERT INTO Feature_blob_chunk VALUES (1, ?, 1, 0, ?)",
        [(index, bytes([index % 256])) for index in range(1, 151)],
    )
    return conn


def test_repository_does_not_apply_autocomplete_limit(repository_db):
    repository = AllSamplesRepository(repository_db)
    samples = repository.available_samples(1)
    assert len(samples) == 150


def test_service_orders_then_limits_repository_samples(repository_db):
    repository = AllSamplesRepository(repository_db)
    request = AllSamplesPlotRequest(
        "c1",
        None,
        GenomicRegion(1, 10),
        ordering=SampleOrdering(source="Sample", column="rank", ascending=True, max_samples=3),
    )
    samples = AllSamplesDataService(repository)._resolve_samples(1, request)
    assert [name for _, name in samples] == ["sample_150", "sample_149", "sample_148"]


def test_service_orders_samples_alphabetically_by_text_metric(repository_db):
    repository = AllSamplesRepository(repository_db)
    ordering = SampleOrdering(source="Explicit_coverage", column="Coverage_class", ascending=True, max_samples=None)

    samples = repository.order_samples(1, [(2, "sample_002"), (1, "sample_001")], ordering)

    assert [name for _, name in samples] == ["sample_001", "sample_002"]


def test_contig_view_can_order_samples_by_parent_mag_metric(repository_db):
    repository = AllSamplesRepository(repository_db)
    ordering = SampleOrdering(source="Explicit_coverage_per_MAG", column="Coverage_mean", ascending=False)

    samples = repository.order_samples(1, [(1, "sample_001"), (2, "sample_002")], ordering)

    assert [name for _, name in samples] == ["sample_002", "sample_001"]


def test_repository_fetches_all_sample_chunks_in_one_query(repository_db):
    repository = AllSamplesRepository(repository_db)
    before = repository.query_count
    chunks = repository.sample_chunks(1, range(1, 151), 1, 0, 0)
    assert repository.query_count - before == 1
    assert len(chunks) == 150


class _FakeRepository:
    query_count = 0

    def resolve_contig(self, name):
        return 1, name, 100

    def available_samples(self, contig_id):
        return [(1, "sample")]

    def order_samples(self, contig_id, samples, ordering):
        return samples

    def feature_styles(self, subplot, table=None):
        if table == "Feature_blob":
            return [FeatureStyle("curve", "#000", 1, 0.2, 1, subplot, "coverage", table, "dense")]
        return []

    def feature_id(self, name):
        return 1

    def blob_settings(self, name):
        return 1.0, 65536, [100, 1000], 1

    def sample_chunks(self, *args):
        return {1: [(0, b"raw")]}


def test_service_returns_plain_data_without_bokeh(monkeypatch):
    import thebigbam.plotting.services.all_samples as service_module

    monkeypatch.setattr(
        service_module,
        "decode_raw_chunks",
        lambda rows, scale, chunk_size: {
            "x": np.array([0, 1, 2]),
            "y": np.array([1.0, 2.0, 3.0]),
            "sparse": False,
        },
    )
    request = AllSamplesPlotRequest("c1", "Coverage", GenomicRegion(1, 10))
    result = AllSamplesDataService(_FakeRepository()).load(request)
    assert result.sample_tracks[0].series[0].data["y"] == (1.0, 2.0, 3.0)


def test_all_samples_minimum_frequency_is_absolute(monkeypatch):
    import thebigbam.plotting.services.all_samples as service_module

    monkeypatch.setattr(
        service_module,
        "decode_raw_chunks",
        lambda rows, scale, chunk_size: {
            "x": np.array([0, 1, 2]),
            "y": np.array([0.2, 0.5, 0.8]),
            "sparse": False,
        },
    )
    request = AllSamplesPlotRequest(
        "c1",
        "Coverage",
        GenomicRegion(1, 10),
        display=DisplayOptions(min_relative_value=0.5),
    )

    result = AllSamplesDataService(_FakeRepository()).load(request)

    assert result.sample_tracks[0].series[0].data["y"] == (0, 0.5, 0.8)


def test_renderer_accepts_plain_data_without_database_calls():
    style = FeatureStyle("curve", "#000", 1, 0.2, 1, "Coverage", "coverage", "Feature_blob", "dense")
    series = FeatureSeries(style, {"x": (1, 2, 3), "y": (1.0, 2.0, 3.0)})
    data = AllSamplesPlotData(
        1,
        "c1",
        100,
        GenomicRegion(1, 10),
        (SampleTrack(1, "sample", (series,)),),
        y_max=3.0,
    )
    _, _, figures = AllSamplesRenderer().render_figures(data, 100, True)
    assert len(figures) == 1
    assert figures[0].title.text == ""
