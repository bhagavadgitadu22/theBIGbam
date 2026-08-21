"""Pure transformations shared by feature-data loading paths."""

from __future__ import annotations


def apply_minimum_frequency(values, minimum_frequency):
    """Zero values below an absolute frequency cutoff without mutating input."""
    values = list(values)
    if minimum_frequency <= 0 or not values:
        return values
    return [value if value >= minimum_frequency else 0 for value in values]


# Kept as an import-compatible alias for callers outside the plotting package.
apply_relative_threshold = apply_minimum_frequency
