"""Deletion letters remain complete and discrepancy windows bound work and glyphs."""

from dataclasses import replace
from pathlib import Path

from bokeh.models import Range1d

from thebigbam.plotting.models.read_candidate import ReadCandidate
from thebigbam.plotting.models.read_events import ReadEvent
from thebigbam.plotting.models.read_plot import ReadPlotData, ReadPlotRequest
from thebigbam.plotting.renderers.read_alignments import render_read_track
from thebigbam.plotting.repositories.read_details import selected_read_details
from thebigbam.plotting.services.read_event_extraction import extract_read_events
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.services.read_reference import ReadReference
from thebigbam.plotting.services.read_summaries import summarize_candidate
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"


def _request(**kwargs):
    return ReadPlotRequest(str(FIXTURE), "NC_002167.1", 39732, 1, 1000, 20, False,
                           frozenset({"primary"}), event_variables=frozenset({"deletions"}), **kwargs)


def test_adjacent_deletions_share_the_md_token_without_losing_letters():
    request = _request()
    for cigar, sequence in (("1M2D1P3D1M", "AA"), ("1M2D1I3D1M", "AGA")):
        read = ReadCandidate("read", 10, 17, cigar, "primary", 0, 60, None)
        events, missing, _ = extract_read_events(read, sequence, "1^ACGTA1", ReadReference(request), request, 100)
        assert [event.reference for event in events] == ["AC", "GTA"]
        assert not missing


def test_incomplete_md_uses_explicit_unknown_bases_instead_of_blank_letters():
    read = ReadCandidate("read", 10, 17, "1M5D1M", "primary", 0, 60, None)
    request = _request()
    events, missing, _ = extract_read_events(read, "AA", "1^AC1", ReadReference(request), request, 100)
    assert events[0].reference == "AC???" and missing


def test_window_limit_skips_details_and_overlays_but_keeps_read_arrows(monkeypatch):
    calls = []

    def load(*args, **kwargs):
        calls.append(1)
        yield from selected_read_details(*args, **kwargs)

    monkeypatch.setattr("thebigbam.plotting.services.read_discrepancies.selected_read_details", load)
    cache = SessionDataCache()
    request = _request(max_discrepancy_window=100)
    data = prepare_read_plot(request, cache)
    assert data.alignments and not any(read.events for read in data.alignments)
    assert calls == []
    plot, _ = render_read_track(data, request, Range1d(1, 1000))
    assert plot.select_one({"name": "read-event-bases"}) is None
    assert plot.renderers[0].data_source.data["name"]
    enabled = replace(request, max_discrepancy_window=2000, profiler=ReadProfiler())
    enriched = prepare_read_plot(enabled, cache)
    assert calls == [1]
    assert enabled.profiler.counts["candidate_cache_hits"] == 1
    assert enabled.profiler.counts.get("bam_examined", 0) == 0
    assert [read.name for read in enriched.alignments] == [read.name for read in data.alignments]


def test_mag_window_limit_is_checked_before_member_local_windows(monkeypatch):
    read = ReadCandidate("read", 10, 30, "20M", "primary", 0, 60, None)
    monkeypatch.setattr("thebigbam.plotting.services.read_plots._fingerprint", lambda _: (1, 1))
    monkeypatch.setattr("thebigbam.plotting.services.read_plots.fetch_read_candidates", lambda *args, **kwargs: iter((read,)))

    def unexpected(*args, **kwargs):
        raise AssertionError("Wide MAG windows must not fetch discrepancy details")

    monkeypatch.setattr("thebigbam.plotting.services.read_discrepancies.selected_read_details", unexpected)
    request = replace(_request(max_discrepancy_window=100), contig_length=550, end=550,
                      members=(("a", 50, 0), ("b", 50, 500)))
    data = prepare_read_plot(request)
    assert len(data.alignments) == 2
    assert not any(read.events for read in data.alignments)


def test_letter_budget_hides_letters_consistently_instead_of_leaving_later_boxes_blank(monkeypatch):
    monkeypatch.setattr("thebigbam.plotting.renderers.read_events.MAX_LETTERS", 3)
    first = summarize_candidate(ReadCandidate("first", 0, 10, "10M", "primary", 0, 60, None))
    first = replace(first, events=(ReadEvent("deletions", 2, 4, 2, reference="AC"),))
    second = replace(first, name="second", alignment_id=("second",))
    request = _request(max_discrepancy_window=100)
    plot, _ = render_read_track(ReadPlotData((first, second), 2), request, Range1d(1, 10))
    letters = plot.select_one({"name": "read-event-letters"})
    assert letters.data_source.data["text"] == [] and not letters.visible
    assert len(plot.select_one({"name": "read-event-bases"}).data_source.data["read"]) == 2
    callback = plot.x_range.js_property_callbacks["change:start"][0]
    assert callback.args["discrepancy_threshold"] == 100
    assert not callback.args["letters_available"]
