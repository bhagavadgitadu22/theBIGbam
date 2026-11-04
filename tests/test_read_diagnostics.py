"""Read diagnostics distinguish APPLY scans, cache reuse and static zooms."""

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from bokeh.models import Range1d

from thebigbam.plotting.application.read_track_pipeline import build_read_tracks
from thebigbam.plotting.models.read_plot import ReadPlotOptions, ReadPlotRequest, ReadTrackConfig
from thebigbam.plotting.renderers.read_alignments import render_read_track
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"


def _request(**kwargs):
    return ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 500, 30, False,
                           frozenset({"primary"}), sample_name="test", **kwargs)


def _reports(capsys):
    return [json.loads(line.removeprefix("[read-profile] "))
            for line in capsys.readouterr().out.splitlines() if line.startswith("[read-profile] ")]


def test_read_profile_scan_sampling_layout_and_cache(capsys):
    cache = SessionDataCache()
    request = _request(profiler=ReadProfiler())
    data = prepare_read_plot(request, cache)
    render_read_track(data, request, Range1d(1, 500))
    report, = _reports(capsys)
    assert report["operation"] == "apply"
    assert report["sample"] == "test"
    assert report["counts"]["bam_examined"] >= report["counts"]["bam_eligible"] == data.fetched
    assert report["plotted_alignments"] <= 30
    assert report["counts"]["cache_misses"] == 1
    assert report["seconds"]["bam_summary"] > 0
    scan_parts = ("bam_open", "bam_region_setup", "bam_iterator", "bam_eligibility", "bam_candidate", "bam_close")
    assert all(report["seconds"][key] > 0 for key in scan_parts)
    assert sum(report["seconds"][key] for key in scan_parts) <= report["seconds"]["bam_scan_including_candidates"]
    assert report["seconds"]["bam_scan_other"] >= 0
    assert report["counts"]["summaries_built"] == len(data.alignments)
    assert data == prepare_read_plot(replace(request, profiler=None))
    assert all(report["seconds"][key] > 0 for key in ("sampling", "row_packing", "glyph_construction", "figure_setup"))

    cached_request = replace(request, profiler=ReadProfiler())
    cached = prepare_read_plot(cached_request, cache)
    render_read_track(cached, cached_request, Range1d(1, 500))
    hit, = _reports(capsys)
    assert cached is data
    assert hit["counts"]["cache_hits"] == 1
    assert hit["counts"].get("bam_examined", 0) == 0
    assert hit["seconds"].get("sampling", 0) == 0
    assert all(key not in hit["seconds"] for key in scan_parts)


def test_zoom_does_not_trigger_read_work_or_profiles(capsys):
    request = _request(profiler=ReadProfiler())
    shared = Range1d(1, 500)
    cache = SessionDataCache()
    render_read_track(prepare_read_plot(request, cache), request, shared)
    _reports(capsys)
    before = dict(request.profiler.counts)
    shared.end = 450
    shared.start, shared.end = 1000, 1500
    assert _reports(capsys) == []
    assert request.profiler.counts == before


def test_read_diagnostics_are_silent_when_disabled(capsys):
    request = _request()
    render_read_track(prepare_read_plot(request), request, Range1d(1, 500))
    assert _reports(capsys) == []


def test_member_tracks_report_scans_and_second_sampling_stage(capsys):
    catalog = SimpleNamespace(by_sample={"sample": SimpleNamespace(path=FIXTURE, available=True)})
    config = ReadTrackConfig(ReadPlotOptions(enabled=True, max_reads=30), catalog,
                             frozenset({"primary"}), {}, diagnose=True)
    build_read_tracks(config, ("sample",), (("NC_002167.1", 39732, 1000),),
                      Range1d(1001, 1500), SessionDataCache())
    report, = _reports(capsys)
    assert report["sample"] == "sample"
    assert report["counts"]["bam_examined"] > 0
    assert report["seconds"]["member_coordinates"] > 0
    assert report["seconds"]["member_resampling"] > 0
    assert report["plotted_alignments"] <= 30
