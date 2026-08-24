import time

from bokeh.models import ColumnDataSource, Range1d

from thebigbam.plotting.renderers.mag_tracks import MagTrackRenderer, visible_mag_track_data


def _members(count=1200, length=1000):
    return tuple((f"contig_{index}", length, index * length) for index in range(count))


def test_large_mag_uses_constant_number_of_renderers_and_sources():
    started = time.perf_counter()
    plot = MagTrackRenderer().render(_members(), 40, Range1d(0, 1_200_000))
    elapsed = time.perf_counter() - started
    sources = plot.select({"type": ColumnDataSource})
    overview = plot.tags[-1]
    assert len(sources) <= 4
    assert len(plot.renderers) <= 4
    assert len(plot.references()) < 100
    assert overview["boundaries"] == 1201
    assert overview["visible_segments"] == 1200
    assert elapsed < 1.0


def test_visible_window_clips_segments_boundaries_and_dots_before_rendering():
    data = visible_mag_track_data(
        _members(10, 100),
        250,
        450,
        {"xs": [100, 300, 500], "colors": ["red", "green", "blue"]},
    )
    assert [segment[0] for segment in data.segments] == ["contig_2", "contig_3", "contig_4"]
    assert data.boundaries == (300, 400)
    assert data.dot_xs == (300.0,)
    assert data.dot_colors == ("green",)


def test_empty_members_return_no_plot():
    assert MagTrackRenderer().render(()) is None


def test_mag_track_double_click_copies_contig_name():
    plot = MagTrackRenderer().render(_members(2), 40, Range1d(0, 2_000))

    callbacks = plot.js_event_callbacks["doubletap"]

    assert len(callbacks) == 1
    assert "navigator.clipboard.writeText(name)" in callbacks[0].code
    assert "tip.textContent = 'Copied'" in callbacks[0].code
    hover = next(tool for tool in plot.tools if tool.__class__.__name__ == "HoverTool")
    assert ("Copy", "Double-click this segment") in hover.tooltips
