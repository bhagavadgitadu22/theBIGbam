from concurrent.futures import ThreadPoolExecutor

import duckdb

from thebigbam.database.database_getters import HistogramResult, resolve_histogram_bins
from thebigbam.plotting.shared.histogram_cache import SharedHistogramCache


def _metadata(source="Metric"):
    return {"Coverage": {"source": source, "columns": {"value": {"type": "numeric"}}}}


def test_exact_histogram_uses_all_values(tmp_path):
    path = tmp_path / "small.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute('CREATE TABLE Metric AS SELECT range::DOUBLE AS "value" FROM range(100)')
    conn.close()

    result = resolve_histogram_bins(str(path), _metadata(), "Coverage", "value", n_bins=10)

    assert isinstance(result, HistogramResult)
    assert result.approximate is False
    assert result.sampled_rows == 100
    assert result.bin_count == 10
    assert sum(result.counts) == 100


def test_large_histogram_is_deterministic_and_bounded(tmp_path):
    path = tmp_path / "large.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute('CREATE TABLE Metric AS SELECT range::DOUBLE AS "value" FROM range(1000000)')
    conn.close()

    first = resolve_histogram_bins(str(path), _metadata(), "Coverage", "value")
    second = resolve_histogram_bins(str(path), _metadata(), "Coverage", "value")

    assert first.approximate is True
    assert 0 < first.sampled_rows <= 100_000
    assert first.edges.tolist() == second.edges.tolist()
    assert first.counts.tolist() == second.counts.tolist()


def test_explicit_metric_prefers_compact_base_table(tmp_path):
    path = tmp_path / "base.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute("CREATE TABLE coverage(value DOUBLE)")
    conn.execute("INSERT INTO coverage VALUES (1), (2), (3)")
    conn.execute('CREATE VIEW Explicit_coverage AS SELECT "value" * 100 AS "value" FROM coverage')
    conn.close()

    result = resolve_histogram_bins(
        str(path), _metadata("Explicit_coverage"), "Coverage", "value"
    )

    assert result.edges[0] == 1
    assert result.edges[-1] == 3


def test_shared_cache_deduplicates_inflight_and_retries_failures():
    cache = SharedHistogramCache()
    calls = []

    def load():
        calls.append(1)
        return "value"

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: cache.get_or_load(("key",), load)[0], range(2)))
    assert results == ["value", "value"]
    assert len(calls) == 1

    attempts = []

    def flaky():
        attempts.append(1)
        if len(attempts) == 1:
            raise RuntimeError("temporary")
        return "recovered"

    try:
        cache.get_or_load(("failure",), flaky)
    except RuntimeError:
        pass
    assert cache.get_or_load(("failure",), flaky)[0] == "recovered"
    assert len(attempts) == 2

    none_calls = []
    assert cache.get_or_load(("none",), lambda: none_calls.append(1))[0] is None
    assert cache.get_or_load(("none",), lambda: none_calls.append(1))[0] is None
    assert len(none_calls) == 2
