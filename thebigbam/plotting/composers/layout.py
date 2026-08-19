"""Final assembly of already-rendered Bokeh figures."""

from __future__ import annotations

import time

from bokeh.layouts import gridplot


def assemble_grid(plots, *, empty_message="No plots to display", enable_timing=False):
    """Combine individual figures into the application's shared plot grid."""
    figures = [plot for plot in plots if plot is not None]
    if not figures:
        raise ValueError(empty_message)
    started = time.perf_counter()
    grid = gridplot([[figure] for figure in figures], merge_tools=True, sizing_mode="stretch_width")
    if enable_timing:
        print(f"[timing]   gridplot ({len(figures)} figures): {time.perf_counter() - started:.3f}s", flush=True)
    return grid
