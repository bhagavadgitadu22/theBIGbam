"""Indexed BAM event retrieval, including origin-only circular and MAG views."""

from dataclasses import replace

import pysam
import pytest

from thebigbam.plotting.models.read_events import EVENT_VARIABLES
from thebigbam.plotting.models.read_plot import ReadPlotRequest
from thebigbam.plotting.services.read_event_coordinates import circular_events
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.shared.data_cache import SessionDataCache


@pytest.fixture
def event_bam(tmp_path):
    path = tmp_path / "events.bam"
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr", "LN": 100}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, start, cigar, sequence, md in (
            ("clipped", 10, "2S3M2I2M2D2M3S", "TTACGGGTACCAAA", "3A1^GT2"),
            ("deletion", 90, "8M5D5M", "A" * 13, "8^ACGTA5"),
            ("crossing", 95, "10M", "AAAAAGAAAA", "5C4"),
        ):
            read = pysam.AlignedSegment()
            read.query_name, read.query_sequence = name, sequence
            read.reference_id, read.reference_start, read.mapping_quality = 0, start, 60
            read.cigarstring = cigar
            read.set_tag("MD", md)
            bam.write(read)
    pysam.index(str(path))
    return path


def _request(path, **kwargs):
    return ReadPlotRequest(str(path), "chr", 100, 1, 100, 20, False,
                           frozenset({"primary"}), event_variables=EVENT_VARIABLES, circular=True, **kwargs)


def test_selected_read_sequences_and_clips_are_loaded_from_indexed_bam(event_bam):
    data = prepare_read_plot(_request(event_bam))
    assert not data.warnings
    clipped = next(read for read in data.alignments if read.name == "clipped")
    assert {event.kind for event in clipped.events} == EVENT_VARIABLES
    insertion = next(event for event in clipped.events if event.kind == "insertions")
    assert insertion.sequence == "GG" and insertion.length == 2


def test_origin_window_includes_crossing_alignments_and_discrepancies(event_bam):
    request = replace(_request(event_bam), end=5)
    data = prepare_read_plot(request)
    assert {read.name for read in data.alignments} == {"deletion", "crossing"}
    assert data.fetched == 2 and not data.warnings
    events = {read.name: circular_events(read, 100) for read in data.alignments}
    assert [(event.start, event.end, event.reference) for event in events["deletion"]] == [(0, 3, "GTA")]
    assert [(event.start, event.sequence, event.reference) for event in events["crossing"]] == [(0, "G", "C")]


def test_circular_discrepancies_follow_mag_offsets(event_bam):
    request = replace(_request(event_bam), contig_length=600, start=501, end=505,
                      members=(("chr", 100, 500),))
    data = prepare_read_plot(request, SessionDataCache())
    events = [event for read in data.alignments for event in read.events]
    assert {(event.kind, event.start, event.end) for event in events} == {
        ("deletions", 500, 503), ("mismatches", 500, 501),
    }
