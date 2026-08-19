"""SQL-free rendering helpers shared by single-sample and MAG composers."""

from __future__ import annotations

import time
from collections import defaultdict

from bokeh.layouts import gridplot
from bokeh.models import Range1d

PRIMARY_RELATIVE_SUBPLOTS = frozenset(
    {
        "Primary alignments",
        "Alignments by strand",
        "Other alignments",
        "Clippings",
        "Indels",
        "Mismatches",
        "Bad orientations",
        "Coverage reduced",
        "Reads termini",
        "Non-inward pairs",
        "Missing mates",
    }
)


def assemble_grid(plots, *, empty_message="No plots to display", enable_timing=False):
    plots = [plot for plot in plots if plot is not None]
    if not plots:
        raise ValueError(empty_message)
    started = time.perf_counter()
    grid = gridplot([[plot] for plot in plots], merge_tools=True, sizing_mode="stretch_width")
    if enable_timing:
        print(f"[timing]   gridplot ({len(plots)} figures): {time.perf_counter() - started:.3f}s", flush=True)
    return grid


def apply_per_feature_y_ranges(feature_subplots) -> None:
    groups = defaultdict(list)
    for feature, subplot, maximum in feature_subplots:
        groups[feature].append((subplot, maximum))
    for items in groups.values():
        maximum = max(value for _subplot, value in items)
        if maximum > 0:
            for subplot, _value in items:
                subplot.y_range = Range1d(0, maximum)


def apply_primary_relative_y_range(feature_subplots, primary_maximum) -> None:
    if primary_maximum <= 0:
        return
    for feature, subplot, _maximum in feature_subplots:
        if feature in PRIMARY_RELATIVE_SUBPLOTS:
            subplot.y_range = Range1d(0, primary_maximum)
