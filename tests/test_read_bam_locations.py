"""Indexed queries may skip physical records; selected offsets must remain exact."""

from dataclasses import replace

import pysam
import pytest

from thebigbam.plotting.models.read_plot import ReadPlotRequest
from thebigbam.plotting.repositories.bam_reads import fetch_read_candidates
from thebigbam.plotting.repositories.read_details import selected_read_details
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler


@pytest.fixture
def gapped_query_bam(tmp_path):
    path = tmp_path / "gapped.bam"
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr", "LN": 1000}, {"SN": "other", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, start, length, md, reference in (
            ("long-a", 0, 200, "100C99", 0), ("skip-a", 10, 5, "5", 0),
            ("long-b", 20, 200, "80C119", 0), ("skip-b", 30, 5, "5", 0),
            ("long-c", 40, 200, "60C139", 0), ("other", 0, 200, "200", 1),
        ):
            record = pysam.AlignedSegment()
            record.query_name, record.query_sequence = name, "A" * length
            record.reference_id, record.reference_start, record.mapping_quality = reference, start, 60
            record.cigarstring = f"{length}M"
            record.set_tag("MD", md)
            bam.write(record)
    pysam.index(str(path))
    return path


@pytest.mark.parametrize("diagnose", [False, True])
def test_offsets_match_when_regional_query_skips_intervening_records(gapped_query_bam, diagnose):
    profiler = ReadProfiler() if diagnose else None
    candidates = tuple(fetch_read_candidates(str(gapped_query_bam), "chr", 100, 110,
                                            frozenset({"primary"}), profiler=profiler))
    assert [read.name for read in candidates] == ["long-a", "long-b", "long-c"]
    details = tuple(selected_read_details(str(gapped_query_bam), "chr", candidates))
    assert [read.name for read, _, _ in details] == [read.name for read in candidates]
    assert all(sequence == "A" * 200 for _, sequence, _ in details)


def test_failed_discrepancy_result_is_not_cached_and_next_apply_recovers(gapped_query_bam, monkeypatch):
    original = selected_read_details
    calls = []

    def load(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ValueError("temporary detail lookup failure")
        yield from original(*args, **kwargs)

    monkeypatch.setattr("thebigbam.plotting.services.read_discrepancies.selected_read_details", load)
    request = ReadPlotRequest(str(gapped_query_bam), "chr", 1000, 101, 110, 100, False,
                              frozenset({"primary"}), event_variables=frozenset({"mismatches"}))
    cache = SessionDataCache()
    failed = prepare_read_plot(request, cache)
    assert not failed.details_loaded and failed.warnings
    assert cache.retained_bytes == 0
    retry = replace(request, profiler=ReadProfiler())
    recovered = prepare_read_plot(retry, cache)
    assert recovered.details_loaded and recovered.warnings == ()
    assert retry.profiler.counts["bam_examined"] == 3
    assert len(recovered.alignments) == 3
    assert all(len(read.events) == 1 for read in recovered.alignments)
    assert prepare_read_plot(request, cache) is recovered
