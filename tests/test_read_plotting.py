"""Read plot selection, sampling, packing and batched rendering."""

from __future__ import annotations

from pathlib import Path
from bokeh.models import CheckboxGroup, Div, HoverTool, RadioButtonGroup, Range1d, Select, Spinner

from thebigbam.plotting.composers.layout import assemble_grid
from thebigbam.plotting.controls.read_plotting import build_read_plot_controls
from thebigbam.plotting.models.read_plot import ReadAlignment, ReadPlotRequest
from thebigbam.plotting.renderers.read_alignments import SELECTED_COLOR, render_read_alignments, render_read_track
from thebigbam.plotting.services.read_layout import pack_read_rows
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.lifecycle import PlotLifecycle

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"
CONTIG = "NC_002167.1"


def _read(name: str, start: int, end: int, *, first: bool = False, paired: bool = False) -> ReadAlignment:
    return ReadAlignment(name, start, end, ((start, end),), True, "primary", 60, f"{end-start}M", None, first, paired)


def test_read_controls_default_off_and_follow_sample_availability():
    sample = Select(value=FIXTURE.stem, options=[FIXTURE.stem, "missing"])
    scope = RadioButtonGroup(labels=["ONE", "ALL"], active=0)
    view = RadioButtonGroup(labels=["MAG", "CONTIG"], active=1)
    controls = build_read_plot_controls(
        max_reads=Spinner(value=1000), link_pairs=CheckboxGroup(labels=["Link paired reads"]),
        plot_height=Spinner(value=400), read_height=Spinner(value=14),
    )

    assert controls.add_reads.active == []
    assert controls.max_reads.value == 1000
    assert controls.link_pairs.active == []
    assert not controls.add_reads.disabled
    sample.value = "missing"
    assert not controls.add_reads.disabled
    sample.value = FIXTURE.stem
    scope.active = 1
    assert controls.section.visible
    scope.active = 0
    assert controls.section.visible
    view.active = 0
    assert controls.section.visible
    view.active = 1
    assert controls.section.visible


def test_indexed_fetch_has_stable_global_cap_and_shared_bokeh_range():
    request = ReadPlotRequest(
        str(FIXTURE), CONTIG, 39732, 1, 500, 40, False, frozenset({"primary"}), {}
    )
    cache = SessionDataCache()
    first = prepare_read_plot(request, cache)
    second = prepare_read_plot(request, cache)
    shared = Range1d(1, 500)
    plot = render_read_alignments(first, request, shared)

    assert first == second
    assert first.fetched > 40
    assert len(first.alignments) == 40
    assert cache.stats().hits == 1
    assert plot.x_range is shared
    assert len(plot.renderers) == 4  # bodies, arrowheads, connectors and messages stay batched
    assert all(read.start < 500 and read.end > 0 for read in first.alignments)


def test_pair_linking_places_mates_on_one_row_with_black_connector():
    reads = (_read("pair", 10, 20, first=True, paired=True), _read("other", 12, 16),
             _read("pair", 30, 40, paired=True))

    layout = pack_read_rows(reads, link_pairs=True)

    pair_rows = {item.row for item in layout.reads if item.read.name == "pair"}
    assert len(pair_rows) == 1
    assert len(layout.connectors) == 1
    assert (layout.connectors[0].left, layout.connectors[0].right) == (20, 30)
    assert layout.row_count == 2


def test_no_selected_read_class_shows_actionable_empty_state():
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 39732, 1, 500, 1000, False, frozenset(), {})
    data = prepare_read_plot(request)

    assert data.alignments == ()
    assert data.message == "No read alignment classes were requested."
    plot, scrollbar = render_read_track(data, request, Range1d(1, 500))
    message = plot.renderers[3]
    assert message.data_source.data["y"] == [(plot.y_range.start + plot.y_range.end) / 2]
    assert message.glyph.text_baseline == "middle"
    assert not scrollbar.visible


def test_read_height_arrow_geometry_and_unlimited_window(monkeypatch):
    forward = _read("forward", 10, 30)
    reverse = ReadAlignment("reverse", 40, 60, ((40, 60),), False, "primary", 60, "20M", None, False, False)
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 39732, 1, 100, 100, False,
                              frozenset({"primary"}), {}, plot_height=450, read_height=20)
    from thebigbam.plotting.models.read_plot import ReadPlotData

    plot = render_read_alignments(ReadPlotData((forward, reverse), 2), request, Range1d(1, 100))
    heads = plot.renderers[1].data_source.data
    assert plot.height == 450
    assert heads["xs"][0][1] > heads["xs"][0][0]  # forward tip is right
    assert heads["xs"][1][1] < heads["xs"][1][0]  # reverse tip is left
    assert heads["ys"][0][0] < heads["ys"][0][1] < heads["ys"][0][2]
    body = plot.renderers[0].data_source.data
    assert body["left"][0] < heads["xs"][0][3] < body["right"][0]
    assert body["left"][1] < heads["xs"][1][3] < body["right"][1]
    assert any("name" in callback.code for callbacks in plot.js_event_callbacks.values() for callback in callbacks)

    monkeypatch.setattr("thebigbam.plotting.services.read_plots._fingerprint", lambda _: (1, 1))
    monkeypatch.setattr("thebigbam.plotting.services.read_plots.fetch_read_candidates",
                        lambda *args: iter((forward, reverse)))
    wide = ReadPlotRequest(str(FIXTURE), CONTIG, 1_000_000, 1, 1_000_000, 1, False,
                           frozenset({"primary"}), {})
    data = prepare_read_plot(wide)
    assert data.message is None
    assert len(data.alignments) == 1
    assert data.fetched == 2


def test_panning_read_track_preserves_sources_and_row_layout():
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 39732, 1, 500, 30, False,
                              frozenset({"primary"}), {})
    cache = SessionDataCache()
    shared = Range1d(1, 500)
    plot = render_read_alignments(prepare_read_plot(request, cache), request, shared)
    before = [renderer.data_source.data for renderer in plot.renderers]

    shared.start = 1000
    shared.end = 1500
    assert all(renderer.data_source.data is source for renderer, source in zip(plot.renderers, before))
    assert not shared._callbacks.get("start")
    assert not shared._callbacks.get("end")


def test_read_scrollbar_and_one_pixel_reads_are_lines():
    from thebigbam.plotting.models.read_plot import ReadPlotData

    reads = tuple(_read(f"read-{index}", 10, 30) for index in range(500))
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 39732, 1, 100, 500, False,
                              frozenset({"primary"}), {}, plot_height=400, read_height=1)
    plot, scrollbar = render_read_track(ReadPlotData(reads, 500), request, Range1d(1, 100))

    assert isinstance(scrollbar, Div)
    assert scrollbar.visible
    assert "overflow-y:scroll" in scrollbar.text
    assert "scrollTop" in scrollbar.text
    assert plot.renderers[1].data_source.data["xs"] == []
    body = plot.renderers[0].data_source.data
    assert body["top"][0] - body["bottom"][0] == 1
    from bokeh.plotting import figure

    coverage = figure(x_range=plot.x_range)
    grid = assemble_grid((coverage, (plot, scrollbar)), show_subplot_titles=False)
    cells = {id(child): (row_index, column) for child, row_index, column in grid.children}
    assert cells[id(coverage)] == (0, 1)
    assert cells[id(plot)] == (1, 1)
    assert cells[id(scrollbar)] == (1, 0)
    assert grid.cols == ["24px", "minmax(0, 1fr)"]
    assert scrollbar.margin == 0
    assert scrollbar in grid.select({"type": Div})
    assert PlotLifecycle.shared_xrange(grid) is plot.x_range


def test_read_plot_has_no_vertical_grid_and_compact_cigar_tooltips():
    from thebigbam.plotting.models.read_plot import ReadPlotData

    reads = (
        ReadAlignment("long", 10, 30, ((10, 30),), True, "primary", 60, "M" * 21, None, False, False),
        ReadAlignment("exact", 40, 60, ((40, 60),), True, "primary", 60, "M" * 20, None, False, False),
    )
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 39732, 1, 100, 2, False, frozenset({"primary"}))
    plot = render_read_alignments(ReadPlotData(reads, 2), request, Range1d(1, 100))

    assert all(not grid.visible for grid in plot.xgrid)
    assert plot.select_one(HoverTool).tooltips == [
        ("Read", "@name"), ("MAPQ", "@mapq"),
        ("Alignment start", "@start{0,0}"), ("Alignment end", "@end{0,0}"),
        ("CIGAR", "@cigar_preview"),
    ]
    cigar_values = plot.renderers[0].data_source.data["cigar_preview"]
    assert cigar_values == ["M" * 20 + "...", "M" * 20]


def test_click_selection_colors_all_plotted_alignments_and_circular_sides():
    from thebigbam.plotting.models.read_plot import ReadPlotData

    def member(start, end, alignment_class="primary"):
        return ReadAlignment("same-name", start, end, ((start, end),), True, alignment_class,
                             60, "15M", None, False, False)

    reads = (
        member(10, 25), member(30, 45, "secondary"), member(50, 65, "supplementary"),
        member(70, 85), member(95, 110), _read("different", 40, 55),
    )
    request = ReadPlotRequest(str(FIXTURE), CONTIG, 100, 1, 100, 100, False,
                              frozenset({"primary", "secondary", "supplementary"}))
    shared = Range1d(1, 100)
    plot, _slider = render_read_track(ReadPlotData(reads, len(reads)), request, shared)
    body = plot.renderers[0].data_source
    selection = plot.js_event_callbacks["tap"][0].args["selection"]
    assert {"tap", "doubletap"} <= set(plot.js_event_callbacks)
    assert any(left < 10 for left, name in zip(body.data["left"], body.data["name"]) if name == "same-name")
    assert any(left > 90 for left, name in zip(body.data["left"], body.data["name"]) if name == "same-name")

    selection.data = {"name": ["same-name"]}
    colors = list(body.data["color"])
    shared.start = 2
    assert selection.data == {"name": ["same-name"]}
    assert body.data["color"] == colors
    callback = selection.js_property_callbacks["change:data"][0]
    assert SELECTED_COLOR in callback.code
    assert "name === selected" in callback.code
