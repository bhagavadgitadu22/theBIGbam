"""Pure transformations shared by feature-data loading paths."""

from __future__ import annotations

from bisect import bisect_left


def add_valid_termini(values, centers):
    """Mark bars whose represented position or bin contains a passing terminus."""
    ordered = sorted(centers)
    first = values.get("first_pos", values["x"])
    last = values.get("last_pos", values["x"])
    values["terminus_valid"] = [
        "Yes" if (index := bisect_left(ordered, lo)) < len(ordered) and ordered[index] <= hi else "No"
        for lo, hi in zip(first, last)
    ]
    return values


def apply_minimum_frequency(values, minimum_frequency):
    """Zero values below an absolute frequency cutoff without mutating input."""
    values = list(values)
    if minimum_frequency <= 0 or not values:
        return values
    return [value if value >= minimum_frequency else 0 for value in values]
