"""Constant-renderer Bokeh overview for MAG contig boundaries."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass

from bokeh.models import (
    BoxZoomTool,
    ColumnDataSource,
    HoverTool,
    Range1d,
    SaveTool,
    WheelZoomTool,
)
from bokeh.plotting import figure


@dataclass(frozen=True)
class VisibleMagTrackData:
    boundaries: tuple[float, ...]
    segments: tuple[tuple[str, int, int, int], ...]
    dot_xs: tuple[float, ...] = ()
    dot_colors: tuple[str, ...] = ()


def visible_mag_track_data(members, visible_start, visible_end, track_dots=None):
    """Select visible members using their sorted offsets before creating Bokeh models."""
    members = tuple(members)
    if not members:
        return VisibleMagTrackData((), ())
    offsets = [int(member[2]) for member in members]
    ends = [int(offset + length) for _name, length, offset in members]
    first = bisect_right(ends, visible_start)
    last = bisect_left(offsets, visible_end, lo=first)
    visible = tuple(
        (str(name), int(length), int(offset), int(offset + length)) for name, length, offset in members[first:last]
    )
    boundary_values = []
    for _name, _length, start, end in visible:
        if visible_start <= start <= visible_end:
            boundary_values.append(start)
        if visible_start <= end <= visible_end:
            boundary_values.append(end)
    boundaries = tuple(dict.fromkeys(boundary_values))
    dot_xs, dot_colors = (), ()
    if track_dots:
        pairs = [
            (float(x), str(color))
            for x, color in zip(track_dots.get("xs", ()), track_dots.get("colors", ()))
            if visible_start <= x <= visible_end
        ]
        dot_xs = tuple(x for x, _color in pairs)
        dot_colors = tuple(color for _x, color in pairs)
    return VisibleMagTrackData(boundaries, visible, dot_xs, dot_colors)


class MagTrackRenderer:
    """Render a MAG overview with a constant number of glyph renderers."""

    def render(self, members, height=30, x_range=None, track_dots=None):
        members = tuple(members)
        if not members:
            return None
        total_length = max(int(offset + length) for _name, length, offset in members)
        shared_range = x_range if x_range is not None else Range1d(0, total_length)
        visible_start, visible_end = float(shared_range.start), float(shared_range.end)
        data = visible_mag_track_data(members, visible_start, visible_end, track_dots)
        tools = "xpan,reset" if x_range is not None else ""
        plot = figure(
            height=int(height),
            sizing_mode="stretch_width",
            x_range=shared_range,
            y_range=Range1d(-1, 1),
            tools=tools,
            toolbar_location=None,
            outline_line_color=None,
        )
        if x_range is not None:
            wheel = WheelZoomTool(dimensions="width")
            plot.add_tools(wheel, BoxZoomTool(dimensions="width"))
            plot.toolbar.active_scroll = wheel
        plot.add_tools(SaveTool())
        plot.yaxis.visible = False
        plot.grid.grid_line_color = None

        backbone = ColumnDataSource({"left": [visible_start], "right": [visible_end]})
        plot.quad(left="left", right="right", bottom=-0.8, top=0.8, source=backbone, color="#DCDCDC", line_color=None)

        boundary_source = ColumnDataSource(
            {"xs": [[x, x] for x in data.boundaries], "ys": [[-0.8, 0.8] for _x in data.boundaries]}
        )
        plot.multi_line(xs="xs", ys="ys", source=boundary_source, color="black", line_width=3)

        if data.dot_xs:
            dot_source = ColumnDataSource(
                {"x": list(data.dot_xs), "y": [0] * len(data.dot_xs), "color": list(data.dot_colors)}
            )
            plot.scatter(
                x="x", y="y", color="color", source=dot_source, size=8, alpha=0.85, line_color="white", line_width=0.5
            )

        if data.segments:
            segment_source = ColumnDataSource(
                {
                    "name": [segment[0] for segment in data.segments],
                    "length": [segment[1] for segment in data.segments],
                    "x0": [segment[2] for segment in data.segments],
                    "x1": [segment[3] for segment in data.segments],
                }
            )
            hover_regions = plot.quad(
                left="x0", right="x1", bottom=-0.8, top=0.8, source=segment_source, fill_alpha=0.0, line_color=None
            )
            plot.add_tools(
                HoverTool(
                    renderers=[hover_regions],
                    tooltips=[
                        ("Contig", "@name"),
                        ("Length", "@length{0,0}"),
                        ("Start", "@x0{0,0}"),
                        ("End", "@x1{0,0}"),
                    ],
                )
            )
        plot.tags.append(
            {
                "component": "mag_overview",
                "boundaries": len(data.boundaries),
                "visible_segments": len(data.segments),
                "data_sources": len(plot.select({"type": ColumnDataSource})),
            }
        )
        return plot
