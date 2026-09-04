import sys
import types

if "thebigbam_rs" not in sys.modules:
    rust_stub = types.ModuleType("thebigbam_rs")
    rust_stub.decode_dense_chunk = lambda data: []
    rust_stub.decode_sparse_chunk = lambda data: ([], [])
    sys.modules["thebigbam_rs"] = rust_stub

from thebigbam.plotting.models.plots import AllSamplesPlotRequest, GenomicRegion
from thebigbam.plotting.services.all_samples import AllSamplesDataService
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.timing import ApplyProfiler


class EmptyRepository:
    def __init__(self):
        self.resolve_calls = 0

    def resolve_contig(self, name):
        self.resolve_calls += 1
        return 1, name, 100

    def available_samples(self, contig_id):
        return [(1, "sample")]


def test_session_cache_is_bounded_and_records_invalidation_reason():
    cache = SessionDataCache(max_entries=2)
    cache.put("a", (1,))
    cache.put("b", (2,))
    assert cache.get("a") == (True, (1,))
    cache.put("c", (3,))
    assert cache.get("b") == (False, None)

    cache.invalidate("filter_change")

    stats = cache.stats()
    assert stats.entries == 0
    assert stats.hits == 1
    assert stats.misses == 1
    assert stats.invalidations == 1
    assert stats.last_invalidation == "filter_change"


def test_session_cache_preserves_selected_namespaces_on_invalidation():
    cache = SessionDataCache()
    cache.put(("mag_context", "m1"), "context")
    cache.put(("mag_overview", "positions"), "annotations")
    cache.put(("feature_data", "m1"), "features")

    cache.invalidate("filter_change", preserve_prefixes=("mag_context", "mag_overview"))

    assert cache.get(("mag_context", "m1")) == (True, "context")
    assert cache.get(("mag_overview", "positions")) == (True, "annotations")
    assert cache.get(("feature_data", "m1")) == (False, None)


def test_all_samples_plain_data_cache_avoids_repeated_repository_work():
    repository = EmptyRepository()
    cache = SessionDataCache()
    profiler = ApplyProfiler("contig_all", "c1")
    request = AllSamplesPlotRequest(
        "c1",
        None,
        GenomicRegion(1, 10),
        data_cache=cache,
        profiler=profiler,
    )

    first = AllSamplesDataService(repository).load(request)
    second = AllSamplesDataService(repository).load(request)

    assert first is second
    assert repository.resolve_calls == 1
    assert profiler.cache_misses == {"all_samples_data": 1}
    assert profiler.cache_hits == {"all_samples_data": 1}
    assert not hasattr(first, "references")


def test_region_changes_miss_and_session_caches_do_not_share_values():
    first_repository = EmptyRepository()
    first_cache = SessionDataCache()
    service = AllSamplesDataService(first_repository)
    service.load(AllSamplesPlotRequest("c1", None, GenomicRegion(1, 10), data_cache=first_cache))
    service.load(AllSamplesPlotRequest("c1", None, GenomicRegion(10, 20), data_cache=first_cache))
    assert first_repository.resolve_calls == 2

    second_repository = EmptyRepository()
    AllSamplesDataService(second_repository).load(
        AllSamplesPlotRequest("c1", None, GenomicRegion(1, 10), data_cache=SessionDataCache())
    )
    assert second_repository.resolve_calls == 1


def test_profiler_reports_phases_and_cache_metrics(capsys):
    profiler = ApplyProfiler("mag_one", "m1", "s1", enabled=True)
    with profiler.phase("blob_retrieval"):
        pass
    profiler.cache("feature_data", False)
    profiler.cache("feature_data", True)
    profiler.report(SessionDataCache().stats())

    output = capsys.readouterr().out
    assert "mode=mag_one subject=m1 sample=s1" in output
    assert "blob_retrieval=" in output
    assert "feature_data=hit:1/miss:1" in output
