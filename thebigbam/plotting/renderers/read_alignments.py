"""Batched Bokeh glyphs for BAM alignments and their interactions."""

from __future__ import annotations

from time import perf_counter

from bokeh.events import DoubleTap, Tap
from bokeh.models import ColumnDataSource, CustomJS, Div, HoverTool, NumeralTickFormatter, Range1d, WheelZoomTool
from bokeh.plotting import figure

from ..models.read_plot import ReadPlotData, ReadPlotRequest
from ..services.read_coordinates import wrap_interval
from ..services.read_layout import pack_read_rows
from .clipboard import COPY_FEEDBACK_JS
from .read_events import render_read_events
from .read_scrollbar import scrollbar_html

COLOR_VARIABLES = {"secondary": "secondary_reads", "supplementary": "supplementary_reads"}
READ_CLASS_COLORS = {"secondary_reads": "#1f77b4", "supplementary_reads": "#B45C1F"}
SELECTED_COLOR = "#e03131"
READ_HIT_TEST_JS = """
    const x = cb_obj.x, y = cb_obj.y;
    let name = null;
    for (let i = 0; i < body.data.name.length; i++) {
        if (body.data.left[i] <= x && x <= body.data.right[i] &&
            body.data.bottom[i] <= y && y <= body.data.top[i]) {
            name = body.data.name[i]; break;
        }
    }
    if (name === null) {
        for (let i = 0; i < heads.data.name.length; i++) {
            if (heads.data.x0[i] <= x && x <= heads.data.x1[i] &&
                heads.data.bottom[i] <= y && y <= heads.data.top[i]) {
                const base = heads.data.xs[i][0], tip = heads.data.xs[i][1];
                const fraction = (x - base) / (tip - base);
                const halfHeight = (heads.data.top[i] - heads.data.bottom[i]) / 2;
                const center = heads.data.ys[i][1];
                if (Math.abs(y - center) <= halfHeight * fraction) {
                    name = heads.data.name[i]; break;
                }
            }
        }
    }
"""


def _color(alignment_class: str, colors) -> str:
    variable = COLOR_VARIABLES.get(alignment_class)
    return colors.get(variable, READ_CLASS_COLORS[variable]) if variable else "#999999"


def _cigar_preview(cigar: str) -> str:
    return cigar[:20] + ("..." if len(cigar) > 20 else "")


def _display_blocks(blocks: tuple[tuple[int, int], ...], contig_length: int) -> tuple[tuple[int, int], ...]:
    """Show both ends of a SAM circular alignment crossing the reference origin."""
    return tuple(span for start, end in blocks for span in wrap_interval(start, end, contig_length))


def render_read_track(data: ReadPlotData, request: ReadPlotRequest, x_range: Range1d):
    setup_started = perf_counter() if request.profiler is not None else 0
    plot = figure(
        title=request.sample_name or "Reads", height=request.plot_height, x_range=x_range, y_range=Range1d(0, 1),
        tools="xpan,ypan,xwheel_zoom,xbox_zoom,reset,save", sizing_mode="stretch_width",
    )
    plot.toolbar.logo = None
    plot.toolbar.active_scroll = plot.select_one(WheelZoomTool)
    plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
    plot.xgrid.visible = False
    plot.yaxis.visible = False
    plot.ygrid.visible = False
    plot.min_border_left = 40
    plot.min_border_right = 10
    row_scrollbar = Div(
        text="",
        width=24,
        min_width=24,
        max_width=24,
        height=request.plot_height,
        sizing_mode="fixed",
        margin=0,
        styles={"padding": "0", "margin": "0", "position": "relative", "z-index": "1"},
    )
    body = ColumnDataSource({key: [] for key in (
        "left", "right", "bottom", "top", "color", "base_color", "name", "mapq",
        "cigar_preview", "start", "end",
    )})
    heads = ColumnDataSource({"xs": [], "ys": [], "color": [], "base_color": [], "name": [], "x0": [], "x1": [],
                              "bottom": [], "top": []})
    selected_read = ColumnDataSource({"name": [""]})
    connectors = ColumnDataSource({"x0": [], "x1": [], "y": []})
    messages = ColumnDataSource({"x": [], "y": [], "text": []})
    bodies = plot.quad(left="left", right="right", bottom="bottom", top="top",
                       fill_color="color", line_color=None, source=body)
    plot.patches(xs="xs", ys="ys", fill_color="color", line_color=None, source=heads)
    plot.segment(x0="x0", y0="y", x1="x1", y1="y", color="black", line_width=1, source=connectors)
    plot.text(x="x", y="y", text="text", text_align="center", text_baseline="middle", source=messages)
    plot.add_tools(HoverTool(renderers=[bodies], tooltips=[
        ("Read", "@name"), ("MAPQ", "@mapq"),
        ("Alignment start", "@start{0,0}"), ("Alignment end", "@end{0,0}"),
        ("CIGAR", "@cigar_preview"),
    ]))
    selected_read.js_on_change("data", CustomJS(args={"body": body, "heads": heads, "selection": selected_read}, code=f"""
        const selected = selection.data.name[0];
        body.data = {{...body.data, color: body.data.name.map((name, i) =>
            selected && name === selected ? '{SELECTED_COLOR}' : body.data.base_color[i])}};
        heads.data = {{...heads.data, color: heads.data.name.map((name, i) =>
            selected && name === selected ? '{SELECTED_COLOR}' : heads.data.base_color[i])}};
    """))

    def redraw(current: ReadPlotData, visible: ReadPlotRequest) -> None:
        profiler = visible.profiler
        started = perf_counter() if profiler is not None else 0
        layout = pack_read_rows(current.alignments, visible.link_pairs, visible.contig_length)
        if profiler is not None:
            profiler.add_time("row_packing", perf_counter() - started)
            profiler.count("rows", layout.row_count)
        started = perf_counter() if profiler is not None else 0
        row_gap = 0 if visible.read_height <= 2 else 1 if visible.read_height < 6 else 3
        row_pitch = visible.read_height + row_gap
        rows_visible = max(1, (visible.plot_height - 55) / row_pitch)
        plot.y_range.end = layout.row_count + 1
        plot.y_range.start = plot.y_range.end - rows_visible
        row_scrollbar.height = visible.plot_height
        row_scrollbar.visible = layout.row_count + 1 > rows_visible
        row_scrollbar.text = scrollbar_html(
            plot.y_range.id, visible.plot_height, layout.row_count, rows_visible,
        )
        sample_label = visible.sample_name or "Reads"
        plot.title.text = f"{sample_label} ({len(current.alignments):,} of {current.fetched:,} eligible reads)"
        messages.data = {
            "x": [(visible.start + visible.end) / 2] if current.message else [],
            "y": [(plot.y_range.start + plot.y_range.end) / 2] if current.message else [],
            "text": [current.message] if current.message else [],
        }
        columns = {key: [] for key in body.data}
        head_columns = {key: [] for key in heads.data}
        span = max(1, visible.end - visible.start)
        for placed in layout.reads:
            read = placed.read
            y = layout.row_count - placed.row - 0.5
            half_body = visible.read_height / (2 * row_pitch)
            bottom, top = y - half_body, y + half_body
            base_color = _color(read.alignment_class, visible.colors)
            color = base_color
            # Keep the alignment body continuous; discrepancy overlays can
            # later render CIGAR-specific gaps and events on top of it.
            blocks = _display_blocks(((read.start, read.end),), visible.contig_length)
            if not blocks:
                continue
            terminal = blocks[-1] if read.forward else blocks[0]
            tip = (terminal[1] if read.forward else terminal[0]) + 0.5
            show_head = visible.read_height >= 6
            if show_head:
                width = min(max(0.5, span / 100), max(0.5, (terminal[1] - terminal[0]) / 2))
                base = tip - width if read.forward else tip + width
                # Extend the head beneath the body: abutting canvas shapes
                # can leave a light antialiasing seam at their shared edge.
                join = base - width / 2 if read.forward else base + width / 2
                head_columns["xs"].append([base, tip, base, join, join])
                head_columns["ys"].append([bottom, y, top, top, bottom])
                head_columns["color"].append(color)
                head_columns["base_color"].append(base_color)
                head_columns["name"].append(read.name)
                head_columns["x0"].append(min(base, tip))
                head_columns["x1"].append(max(base, tip))
                head_columns["bottom"].append(bottom)
                head_columns["top"].append(top)
            for block in blocks:
                left, right = block[0] + 0.5, block[1] + 0.5
                if show_head and block == terminal:
                    right = base if read.forward else right
                    left = base if not read.forward else left
                if right <= left:
                    continue
                values = (left, right, bottom, top, color, base_color, read.name, read.mapq,
                          _cigar_preview(read.cigar), read.start + 1, read.end)
                for key, value in zip(columns, values):
                    columns[key].append(value)
        body.data = columns
        heads.data = head_columns
        connectors.data = {
            "x0": [item.left + 0.5 for item in layout.connectors],
            "x1": [item.right + 0.5 for item in layout.connectors],
            "y": [layout.row_count - item.row - 0.5 for item in layout.connectors],
        }
        render_read_events(plot, layout, visible, row_pitch, current.warnings)
        if profiler is not None:
            profiler.add_time("glyph_construction", perf_counter() - started)
            profiler.count("body_glyphs", len(columns["name"]))
            profiler.count("arrow_glyphs", len(head_columns["name"]))
            profiler.count("pair_links", len(layout.connectors))
            profiler.report(visible, current)

    if request.profiler is not None:
        request.profiler.add_time("figure_setup", perf_counter() - setup_started)
    redraw(data, request)
    plot.js_on_event(Tap, CustomJS(args={"body": body, "heads": heads, "selection": selected_read},
                                   code=READ_HIT_TEST_JS + "selection.data = {name: [name || '']};"))
    plot.js_on_event(DoubleTap, CustomJS(args={"body": body, "heads": heads},
                                         code=READ_HIT_TEST_JS + "if (name === null) return;" + COPY_FEEDBACK_JS))
    return plot, row_scrollbar


def render_read_alignments(data: ReadPlotData, request: ReadPlotRequest, x_range: Range1d):
    """Return the read figure for callers that do not need its scroll control."""
    return render_read_track(data, request, x_range)[0]
