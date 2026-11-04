"""CIGAR/MD discrepancies, circular projection and batched read overlays."""

from dataclasses import replace
from pathlib import Path

import duckdb
from bokeh.models import Range1d

from thebigbam.plotting.models.read_candidate import ReadCandidate
from thebigbam.plotting.models.read_events import EVENT_VARIABLES, ReadEvent
from thebigbam.plotting.models.read_plot import ReadPlotData, ReadPlotRequest
from thebigbam.plotting.renderers.read_alignments import render_read_track
from thebigbam.plotting.renderers.read_events import sequence_preview
from thebigbam.plotting.repositories.bam_reads import fetch_read_candidates
from thebigbam.plotting.repositories.genome_tracks import GenomeTrackRepository
from thebigbam.plotting.repositories.read_details import selected_read_details
from thebigbam.plotting.services.read_event_coordinates import circular_events, project_events
from thebigbam.plotting.services.read_event_extraction import extract_read_events
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.services.read_reference import ReadReference
from thebigbam.plotting.services.read_summaries import summarize_candidate
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"


def _request(**kwargs):
    return ReadPlotRequest(str(FIXTURE), "NC_002167.1", 100, 1, 100, 20, False,
                           frozenset({"primary"}), event_variables=EVENT_VARIABLES, **kwargs)


def test_cigar_and_md_extract_all_event_types_in_reference_orientation():
    read = ReadCandidate("read", 10, 19, "2S3M2I2M2D2M3S", "primary", 16, 60, None)
    request = _request()
    events, missing, limited = extract_read_events(read, "TTACGGGTACCAAA", "3A1^GT2", ReadReference(request), request, 100)
    assert not missing and not limited
    assert [(event.kind, event.start, event.end) for event in events] == [
        ("left_clippings", 10, 10), ("insertions", 13, 13), ("mismatches", 13, 14),
        ("deletions", 15, 17), ("right_clippings", 19, 19),
    ]
    assert events[1].sequence == "GG"
    assert events[2].sequence == "T" and events[2].reference == "A"
    assert events[3].reference == "GT"
    assert events[4].sequence == "AAA"


def test_reference_skip_does_not_shift_md_mismatches():
    read = ReadCandidate("spliced", 0, 11, "3M5N3M", "primary", 0, 60, None)
    request = _request()
    events, _, _ = extract_read_events(read, "AAAGAA", "3A2", ReadReference(request), request, 100)
    assert [(event.start, event.sequence, event.reference) for event in events] == [(8, "G", "A")]


def test_explicit_mismatch_cigar_without_md_or_reference():
    read = ReadCandidate("read", 0, 3, "1=1X1=", "primary", 0, 60, None)
    request = _request()
    events, _, _ = extract_read_events(read, "AGA", None, ReadReference(request), request, 100)
    assert len(events) == 1 and events[0].sequence == "G"


def test_reference_fallback_and_event_budget():
    class Reference:
        def base(self, position):
            return "A"

    read = ReadCandidate("read", 0, 5, "5M", "primary", 0, 60, None)
    events, missing, limited = extract_read_events(read, "GGGGG", None, Reference(), _request(), 2)
    assert len(events) == 2 and limited and not missing


def test_database_reference_fallback_is_circular_and_cached():
    conn = duckdb.connect(":memory:")
    try:
        conn.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR)")
        conn.execute("CREATE TABLE Contig_sequence (Contig_id INTEGER, Sequence VARCHAR)")
        conn.execute("INSERT INTO Contig VALUES (1, 'NC_002167.1')")
        conn.execute("INSERT INTO Contig_sequence VALUES (1, ?)", ["A" * 100])
        source = GenomeTrackRepository(conn)
        request = _request(reference_source=source)
        reference = ReadReference(request, SessionDataCache())
        read = ReadCandidate("read", 0, 5, "1M2D2M", "primary", 0, 60, None)
        events, missing, _ = extract_read_events(read, "GAA", None, reference, request, 100)
        assert not missing
        assert [(event.kind, event.reference) for event in events] == [("mismatches", "A"), ("deletions", "AA")]
        assert reference.base(100) == "A"
        assert source.query_count == 1
    finally:
        conn.close()


def test_missing_md_and_reference_do_not_fabricate_mismatches():
    request = _request()
    read = ReadCandidate("read", 0, 3, "3M", "primary", 0, 60, None)
    events, missing, _ = extract_read_events(read, "AAA", None, ReadReference(request), request, 100)
    assert events == () and missing


def test_hard_clips_have_sizes_without_fabricated_sequence():
    read = ReadCandidate("read", 0, 3, "25H3M30H", "primary", 0, 60, None)
    request = _request()
    events, _, _ = extract_read_events(read, "AAA", "3", ReadReference(request), request, 100)
    assert [event.length for event in events] == [25, 30]
    assert all(event.clipping == "hard" and event.sequence == "" for event in events)


def test_circular_and_mag_event_coordinates_keep_deleted_bases():
    read = summarize_candidate(ReadCandidate("read", 95, 110, "15M", "primary", 0, 60, None))
    read = replace(read, events=(ReadEvent("deletions", 98, 103, 5, reference="ACGTA"),
                                 ReadEvent("mismatches", 104, 105, 1, "G", "A")))
    assert [(event.start, event.end, event.reference) for event in circular_events(read, 100)] == [
        (98, 100, "AC"), (0, 3, "GTA"), (4, 5, "A"),
    ]
    shifted = project_events(read.events, 100, 110, 100, 500, read_start=95)
    assert [(event.start, event.end) for event in shifted] == [(500, 503), (504, 505)]


def test_tooltip_previews_are_limited_to_twenty_bases():
    assert sequence_preview("A" * 20) == "A" * 20
    assert sequence_preview("A" * 21) == "A" * 20 + "..."
    assert sequence_preview("A" * 20, 100) == "A" * 20 + "..."


def test_overlays_follow_variables_colors_and_discrepancy_window():
    read = summarize_candidate(ReadCandidate("read", 10, 30, "20M", "primary", 0, 60, None))
    read = replace(read, events=(ReadEvent("insertions", 15, 15, 40, "A" * 40),
                                 ReadEvent("deletions", 20, 22, 2, reference="GT")))
    request = replace(_request(max_discrepancy_window=100), colors={"insertions": "#123456"})
    plot, _ = render_read_track(ReadPlotData((read,), 1), request, Range1d(1, 100))
    bars = plot.select_one({"name": "read-event-boundaries"}).data_source.data
    assert bars["color"] == ["#123456"]
    assert bars["sequence"] == ["A" * 20 + "..."]
    marker = plot.select_one({"name": "read-event-boundaries"}).glyph
    assert marker.width_units == "screen" and marker.width == 10
    thin, _ = render_read_track(ReadPlotData((read,), 1), replace(request, read_height=1), Range1d(1, 100))
    assert thin.select_one({"name": "read-event-boundaries"}).glyph.width == marker.width
    letters = plot.select_one({"name": "read-event-letters"})
    assert letters.data_source.data["text"] == ["G", "T"]
    assert letters.visible
    assert plot.x_range.js_property_callbacks["change:start"]
    assert not plot.x_range._callbacks.get("start")
    zoomed, _ = render_read_track(ReadPlotData((read,), 1), request, Range1d(10, 15))
    assert zoomed.select_one({"name": "read-event-letters"}).visible


def test_selected_bam_details_use_offsets_and_candidate_cache():
    candidates = tuple(fetch_read_candidates(str(FIXTURE), "NC_002167.1", 0, 500, frozenset({"primary"})))
    chosen = (candidates[0], candidates[5], candidates[-1])
    details = tuple(selected_read_details(str(FIXTURE), "NC_002167.1", chosen))
    assert {read.name for read, _, _ in details} == {read.name for read in chosen}
    assert all(sequence and md is not None for _, sequence, md in details)
    request = replace(_request(), contig_length=39732, end=500, event_variables=frozenset())
    cache = SessionDataCache()
    base = prepare_read_plot(request, cache)
    enriched_request = replace(request, event_variables=EVENT_VARIABLES, profiler=ReadProfiler())
    enriched = prepare_read_plot(enriched_request, cache)
    assert [read.name for read in enriched.alignments] == [read.name for read in base.alignments]
    assert enriched_request.profiler.counts["candidate_cache_hits"] == 1
    assert enriched_request.profiler.counts.get("bam_examined", 0) == 0
    assert not any("Could not load" in warning for warning in enriched.warnings)
