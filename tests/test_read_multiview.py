"""Read tracks across multiple samples and concatenated MAG contigs."""

from dataclasses import replace
from types import SimpleNamespace

from bokeh.models import Range1d

from thebigbam.plotting.application.read_track_pipeline import build_read_tracks
from thebigbam.plotting.models.read_plot import (
    ReadAlignment,
    ReadPlotData,
    ReadPlotOptions,
    ReadPlotRequest,
    ReadTrackConfig,
)
from thebigbam.plotting.services.read_plots import prepare_read_plot


def test_mag_reads_query_local_regions_and_wrap_before_offset(monkeypatch):
    calls = []

    def fetch(path, contig, start, end, classes):
        calls.append((contig, start, end))
        return iter((ReadAlignment("read", 90, 110, ((90, 110),), True, "primary", 60,
                                   "20M", None, False, False),))

    monkeypatch.setattr("thebigbam.plotting.services.read_plots.fetch_read_candidates", fetch)
    monkeypatch.setattr("thebigbam.plotting.services.read_plots._fingerprint", lambda _: (1, 1))
    request = ReadPlotRequest("fake.bam", "mag", 1100, 1001, 1100, 100, False, frozenset({"primary"}),
                              members=(("first", 1000, 0), ("second", 100, 1000)))
    data = prepare_read_plot(request)
    assert calls == [("second", 0, 100)]
    assert {(read.start, read.end) for read in data.alignments} == {(1000, 1010), (1090, 1100)}
    assert data.fetched == 1
    limited = prepare_read_plot(replace(request, max_reads=1))
    assert len(limited.alignments) == 1


def test_multiple_samples_share_range_without_refresh_callbacks(monkeypatch):
    monkeypatch.setattr("thebigbam.plotting.application.read_track_pipeline.prepare_read_plot",
                        lambda *args: ReadPlotData((), 0))
    catalog = SimpleNamespace(by_sample={
        name: SimpleNamespace(path=f"{name}.bam", available=True) for name in ("one", "two")
    })
    config = ReadTrackConfig(ReadPlotOptions(enabled=True), catalog, frozenset({"primary"}), {})
    xrange = Range1d(1, 100)
    figures = build_read_tracks(config, ("two", "one"), (("contig", 100, 0),), xrange, None)
    assert len(figures) == 2
    assert "two" in figures[0][0].title.text
    assert "one" in figures[1][0].title.text
    assert all(track[0].x_range is xrange for track in figures)
    assert not xrange._callbacks.get("start")
    assert not xrange._callbacks.get("end")
