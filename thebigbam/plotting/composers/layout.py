"""Final assembly of already-rendered Bokeh figures."""

from __future__ import annotations

import time

from bokeh.layouts import gridplot


def assemble_grid(
    plots,
    *,
    empty_message="No plots to display",
    enable_timing=False,
    show_subplot_titles=True,
    gene_map=None,
    freeze_gene_map=False,
):
    """Combine individual figures into the application's shared plot grid."""
    figures = [plot for plot in plots if plot is not None]
    if not figures:
        raise ValueError(empty_message)
    frozen_gene_map = freeze_gene_map and gene_map is not None and any(plot is gene_map for plot in figures)
    if frozen_gene_map:
        # Keep the map inside the same Bokeh grid as the tracks so Bokeh
        # computes their shared x-axis and frame geometry together.
        gene_map.styles = {
            **gene_map.styles,
            "position": "sticky",
            "top": "0px",
            "z-index": "10",
            "background-color": "white",
        }
    # Keep the scrollbar at the visible left edge, with all plotting frames
    # in the next column so their genomic axes remain aligned.
    has_scrollbars = any(isinstance(plot, tuple) for plot in figures)
    rows = [
        [plot[1], plot[0]] if isinstance(plot, tuple)
        else [None, plot] if has_scrollbars else [plot]
        for plot in figures
    ]
    if not show_subplot_titles:
        for cells in rows:
            for figure in cells:
                for legend in getattr(figure, "legend", ()):
                    legend.visible = False
    started = time.perf_counter()
    grid = gridplot(rows, merge_tools=True, sizing_mode="stretch_width")
    if has_scrollbars:
        scrollbar_width = max(plot[1].width for plot in figures if isinstance(plot, tuple))
        grid.cols = [f"{scrollbar_width}px", "minmax(0, 1fr)"]
    if enable_timing:
        print(f"[timing]   gridplot ({len(figures)} figures): {time.perf_counter() - started:.3f}s", flush=True)
    return grid
