"""Batched discrepancy rectangles, boundary bars and optional base letters."""

from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Title

from ..services.read_discrepancies import MAX_READ_EVENTS
from ..services.read_event_coordinates import circular_events
from .genome_tracks import NucleotideRenderer

EVENT_COLORS = {"insertions": "#e50001", "deletions": "#00E5E4",
                "left_clippings": "#8e43e7", "right_clippings": "#9CE743"}
EVENT_LABELS = {"insertions": "Insertion", "deletions": "Deletion", "mismatches": "Mismatch",
                "left_clippings": "Left clipping", "right_clippings": "Right clipping"}
MAX_LETTERS = 50_000
BOUNDARY_WIDTH_PX = 10


def sequence_preview(sequence, length=None):
    if not sequence:
        return "Unavailable"
    return sequence[:20] + ("..." if max(len(sequence), length or 0) > 20 else "")


def render_read_events(plot, layout, request, row_pitch, warnings=()):
    if not request.event_variables or abs(plot.x_range.end - plot.x_range.start) > request.max_discrepancy_window:
        return
    tooltip_keys = ("read", "kind", "size", "sequence", "reference", "clipping", "position")
    rectangles = {key: [] for key in (*tooltip_keys, "left", "right", "bottom", "top", "fill", "color")}
    bars = {key: [] for key in (*tooltip_keys, "x", "y", "height", "color")}
    letters = {key: [] for key in ("x", "y", "text", "color")}
    can_show_letters = request.read_height >= 6
    limited, events_limited = False, False
    for placed in layout.reads:
        if len(rectangles["read"]) + len(bars["read"]) >= MAX_READ_EVENTS:
            events_limited = True
            break
        read = placed.read
        y = layout.row_count - placed.row - 0.5
        half = request.read_height / (2 * row_pitch)
        for event in circular_events(read, request.contig_length):
            if event.kind not in request.event_variables:
                continue
            if len(rectangles["read"]) + len(bars["read"]) >= MAX_READ_EVENTS:
                events_limited = True
                break
            color = (NucleotideRenderer.COLORS.get(event.sequence.upper(), "#999999")
                     if event.kind == "mismatches" else request.colors.get(event.kind, EVENT_COLORS[event.kind]))
            details = dict(read=read.name, kind=EVENT_LABELS[event.kind], size=event.length,
                           sequence=sequence_preview(event.sequence, event.length), reference=sequence_preview(event.reference),
                           clipping=event.clipping or "—", position=event.start + 1)
            if event.start == event.end:
                values = dict(details, x=event.start + 0.5, y=y, height=2 * half, color=color)
                for key, value in values.items():
                    bars[key].append(value)
                continue
            values = dict(details, left=event.start + 0.5, right=event.end + 0.5, bottom=y - half, top=y + half,
                          fill="white" if event.kind == "deletions" else color, color=color)
            for key, value in values.items():
                rectangles[key].append(value)
            if not can_show_letters:
                continue
            text = event.reference if event.kind == "deletions" else event.sequence or "?"
            limited |= len(text) < event.end - event.start
            if limited:
                continue
            for index in range(min(event.end - event.start, len(text))):
                if len(letters["text"]) >= MAX_LETTERS:
                    limited = True
                    break
                letters["x"].append(event.start + index + 1)
                letters["y"].append(y)
                letters["text"].append(text[index])
                letters["color"].append(color if event.kind == "deletions" else "white")
    boxes = plot.quad(left="left", right="right", bottom="bottom", top="top", fill_color="fill",
                      line_color="color", line_width=1, source=ColumnDataSource(rectangles), name="read-event-bases")
    boundaries = plot.rect(x="x", y="y", width=BOUNDARY_WIDTH_PX, width_units="screen", height="height",
                           height_units="data", fill_color="color", line_color=None,
                           source=ColumnDataSource(bars), name="read-event-boundaries")
    if limited:
        letters = {key: [] for key in letters}
    text = plot.text(x="x", y="y", text="text", text_color="color", text_align="center", text_baseline="middle",
                     text_font_size=f"{min(10, max(5, request.read_height * .65))}px",
                     source=ColumnDataSource(letters), name="read-event-letters")
    text.visible = not limited and can_show_letters
    visibility = CustomJS(args={"letters": text, "range": plot.x_range,
                                 "tall_enough": can_show_letters, "letters_available": not limited,
                                 "boxes": boxes, "boundaries": boundaries, "discrepancy_threshold": request.max_discrepancy_window}, code="""
        const span = Math.abs(range.end - range.start);
        const showEvents = span <= discrepancy_threshold;
        if (boxes.visible !== showEvents) boxes.visible = showEvents;
        if (boundaries.visible !== showEvents) boundaries.visible = showEvents;
        const visible = showEvents && letters_available && tall_enough;
        if (letters.visible !== visible) letters.visible = visible;
    """)
    plot.x_range.js_on_change("start", visibility)
    plot.x_range.js_on_change("end", visibility)
    plot.add_tools(HoverTool(renderers=[boxes, boundaries], tooltips=[
        ("Read", "@read"), ("Event", "@kind"), ("Position", "@position{0,0}"), ("Length (bp)", "@size{0,0}"),
        ("Sequence", "@sequence"), ("Reference", "@reference"), ("Clipping", "@clipping"),
    ]))
    messages = list(warnings)
    if limited:
        messages.append("Base letters hidden at this density; reduce reads or window and click APPLY.")
    if events_limited:
        messages.append(f"Displayed discrepancies limited to {MAX_READ_EVENTS:,} events.")
    if messages:
        plot.add_layout(Title(text=" ".join(messages), text_font_size="9pt"), "below")
    if request.profiler is not None:
        request.profiler.count("event_rectangles", len(rectangles["read"]))
        request.profiler.count("event_bars", len(bars["read"]))
        request.profiler.count("event_letters", len(letters["text"]))
