"""Bokeh range projection shared by single-sample and MAG renderers."""

from __future__ import annotations

from collections import defaultdict

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
