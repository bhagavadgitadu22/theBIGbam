"""Map selected variable subplots to read classes and DB-backed plots."""

from __future__ import annotations

from typing import Iterable, Mapping

READ_CLASS_BY_VARIABLE = {
    "primary_reads": "primary",
    "primary_reads_plus_only": "primary",
    "primary_reads_minus_only": "primary",
    "secondary_reads": "secondary",
    "supplementary_reads": "supplementary",
}


def read_classes_for_subplots(
    subplots: Iterable[str], subplot_to_varnames: Mapping[str, list[str]]
) -> frozenset[str]:
    return frozenset(
        READ_CLASS_BY_VARIABLE[name]
        for subplot in subplots
        for name in subplot_to_varnames.get(subplot, ())
        if name in READ_CLASS_BY_VARIABLE
    )


def database_features(
    selected: Iterable[str], db_subplots: frozenset[str], read_subplots: frozenset[str]
) -> list[str]:
    """Drop BAM-only choices before the existing database composition path."""
    return [name for name in selected if name in db_subplots or name not in read_subplots]
