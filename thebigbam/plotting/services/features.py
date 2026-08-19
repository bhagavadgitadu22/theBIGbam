"""Pure transformations shared by feature-data loading paths."""

from __future__ import annotations


def apply_relative_threshold(values, minimum_relative_value):
    """Zero values below a fraction of the series maximum without mutating input."""
    values = list(values)
    if minimum_relative_value <= 0 or not values:
        return values
    maximum = max(values)
    if maximum <= 0:
        return values
    threshold = minimum_relative_value * maximum
    return [value if value >= threshold else 0 for value in values]
