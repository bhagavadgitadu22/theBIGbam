"""Candidate reuse preserves spatial sampling while bounding retained memory."""

import math
import random
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import pysam
import pytest

from thebigbam.plotting.models.read_candidate import ReadCandidate
from thebigbam.plotting.models.read_plot import ReadPlotRequest
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.services.read_sampling import build_read_pool, sample_reads
from thebigbam.plotting.services.read_summaries import summarize_candidate
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler
from thebigbam.plotting.shared.read_identity import read_identity, read_score

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"


def _candidates():
    reads = [ReadCandidate(f"r{i}", (i * i * 173) % 1_000_000, (i * i * 173) % 1_000_000 + 34,
                           "10M2I5M4D15M", "primary", 0, 60, None) for i in range(500)]
    random.Random(17).shuffle(reads)
    return reads


def _reference_sample(reads, end, limit):
    """Full-population reference for the previous spatial round-robin rule."""
    bins = defaultdict(list)
    count = max(1, math.ceil(end / 1000))
    quota = math.ceil(limit / count)
    for read in reads:
        bins[min(count - 1, read.start // 1000)].append(read)
    def key(read):
        return read_score(read), read_identity(read)
    ranked = [sorted(bins[index], key=key)[:quota] for index in sorted(bins)]
    selected = {}
    for rank in range(quota):
        for group in ranked:
            if len(group) < rank + 1:
                continue
            if len(selected) < limit:
                selected[read_identity(group[rank])] = group[rank]
    for read in sorted(reads, key=key)[:limit]:
        if len(selected) < limit:
            selected.setdefault(read_identity(read), read)
    return tuple(sorted(selected.values(), key=lambda read: (read.start, read.end, read.name, read.is_read1)))


@pytest.mark.parametrize("end", [1_000_000, 100_000_000])
@pytest.mark.parametrize("limit", [1, 10, 64, 128])
def test_reservoir_preserves_full_population_selection(end, limit):
    reads = _candidates()
    pool = build_read_pool(reads, 0, end, 128)
    direct, total = sample_reads(reads, 0, end, limit)
    reused, _ = sample_reads(pool.candidates, 0, end, limit)
    assert reused == direct == _reference_sample(reads, end, limit)
    assert total == pool.fetched == len(reads)
    assert len(pool.candidates) <= 3 * pool.capacity


def test_count_changes_reuse_candidates_and_build_only_selected_summaries():
    cache = SessionDataCache()
    request = ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 39732, 10, False,
                              frozenset({"primary"}))
    first = replace(request, profiler=ReadProfiler())
    prepare_read_plot(first, cache)
    assert first.profiler.counts["bam_examined"] > 0
    assert first.profiler.counts["summaries_built"] == 10
    for limit in (100, 1000):
        changed = replace(request, max_reads=limit, profiler=ReadProfiler())
        reused = prepare_read_plot(changed, cache)
        assert reused == prepare_read_plot(replace(changed, profiler=None))
        assert changed.profiler.counts["candidate_cache_hits"] == 1
        assert changed.profiler.counts.get("bam_examined", 0) == 0
        assert changed.profiler.counts["summaries_built"] == limit
    assert cache.retained_bytes <= cache.max_retained_bytes


def test_oversize_candidates_are_not_cached():
    cache = SessionDataCache(max_retained_bytes=0)
    request = ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 500, 10, False,
                              frozenset({"primary"}), profiler=ReadProfiler())
    prepare_read_plot(request, cache)
    assert request.profiler.counts["candidate_cache_oversize"] == 1
    changed = replace(request, max_reads=20, profiler=ReadProfiler())
    prepare_read_plot(changed, cache)
    assert changed.profiler.counts["bam_examined"] > 0
    assert cache.retained_bytes == 0


def test_fingerprint_change_invalidates_candidate_reuse(monkeypatch):
    fingerprint = [1, 1]
    monkeypatch.setattr("thebigbam.plotting.services.read_plots._fingerprint", lambda _: tuple(fingerprint))
    request = ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 500, 10, False, frozenset({"primary"}))
    cache = SessionDataCache()
    prepare_read_plot(request, cache)
    fingerprint[1] = 2
    changed = replace(request, max_reads=20, profiler=ReadProfiler())
    prepare_read_plot(changed, cache)
    assert changed.profiler.counts["candidate_cache_misses"] == 1
    assert changed.profiler.counts["bam_examined"] > 0


def test_request_above_reservoir_capacity_rescans(monkeypatch):
    monkeypatch.setattr("thebigbam.plotting.services.read_plots.REUSABLE_READ_CAPACITY", 20)
    request = ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 500, 10, False, frozenset({"primary"}))
    cache = SessionDataCache()
    prepare_read_plot(request, cache)
    changed = replace(request, max_reads=100, profiler=ReadProfiler())
    data = prepare_read_plot(changed, cache)
    assert changed.profiler.counts["bam_examined"] > 0
    assert changed.profiler.counts["candidate_cache_misses"] == 1
    assert len(data.alignments) == 100


def test_candidate_hashes_are_reused(monkeypatch):
    calls = []

    def score(read):
        calls.append(read.name)
        return read_score(read)

    monkeypatch.setattr("thebigbam.plotting.models.read_candidate.read_score", score)
    reads = _candidates()
    pool = build_read_pool(reads, 0, 1_000_000, 128)
    sample_reads(pool.candidates, 0, 1_000_000, 10)
    sample_reads(pool.candidates, 0, 1_000_000, 100)
    assert len(calls) == len(reads)


def test_deferred_blocks_match_pysam():
    record = pysam.AlignedSegment()
    record.reference_start = 7
    record.cigarstring = "5S10M2I5M4D3M6N2=1X4H"
    read = ReadCandidate("test", 7, record.reference_end, record.cigarstring, "primary", 81, 42, 100)
    summary = summarize_candidate(read)
    assert summary.blocks == tuple(record.get_blocks())
    assert summary.is_paired and summary.is_read1 and not summary.forward


def test_weighted_lru_eviction_and_invalidation():
    cache = SessionDataCache(max_retained_bytes=100)
    assert cache.put(("read_candidates", "a"), "a", retained_bytes=40)
    assert cache.put(("read_candidates", "b"), "b", retained_bytes=40)
    cache.get(("read_candidates", "a"))
    assert cache.put(("read_candidates", "c"), "c", retained_bytes=40)
    assert not cache.get(("read_candidates", "b"))[0]
    assert cache.retained_bytes == 80
    assert not cache.put("oversize", "x", retained_bytes=101)
    cache.invalidate("test", preserve_prefixes=("read_candidates",))
    assert cache.retained_bytes == 80
    cache.invalidate("clear")
    assert cache.retained_bytes == 0
