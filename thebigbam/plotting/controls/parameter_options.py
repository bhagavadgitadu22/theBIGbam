"""Pure option catalogs shared by filtering and plotting controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


def format_column(column: str) -> tuple[str, str]:
    return column, column.replace("_", " ").replace("percentage", "(%)")


@dataclass(frozen=True)
class ParameterOptionCatalog:
    columns_raw: dict[str, list[str]]
    columns: dict[str, list[tuple[str, str]]]
    metrics: dict[str, list[tuple[str, str]]]
    sample_metrics: dict[str, list[tuple[str, str]]]
    mag_metrics: dict[str, list[tuple[str, str]]]
    mag_categories: list[str]
    mag_sources: dict[str, str]
    sample_contig_categories: list[str]
    sample_mag_categories: list[str]
    sample_sources: dict[str, str]

    @classmethod
    def from_filtering_metadata(cls, metadata: Mapping[str, Mapping[str, Any]]) -> "ParameterOptionCatalog":
        columns_raw = {category: list(info.get("columns", {})) for category, info in metadata.items()}
        columns = {
            category: [format_column(column) for column in category_columns]
            for category, category_columns in columns_raw.items()
        }
        metrics = {
            category: [
                format_column(column)
                for column, info in category_info.get("columns", {}).items()
                if info.get("type") == "numeric"
            ]
            for category, category_info in metadata.items()
        }
        # Ordering supports both numeric and text columns.  Use the exact same
        # column projection as Filtering so shared categories cannot drift.
        sample_metrics = {category: list(values) for category, values in columns.items()}
        mag_metrics = {category: list(values) for category, values in columns.items()}
        mag_excluded = {"Sample", "Annotations", "MAG", "MAG coverage", "MAG misassembly", "MAG microdiversity"}
        mag_categories = [category for category in metadata if category not in mag_excluded]
        sample_contig_allowed = {
            "Coverage",
            "Misassembly",
            "Microdiversity",
            "Side misassembly",
            "Topology",
            "Termini",
        }
        sample_mag_allowed = {"MAG coverage", "MAG misassembly", "MAG microdiversity"}
        return cls(
            columns_raw=columns_raw,
            columns=columns,
            metrics=metrics,
            sample_metrics=sample_metrics,
            mag_metrics=mag_metrics,
            mag_categories=mag_categories,
            mag_sources={category: metadata[category]["source"] for category in mag_categories},
            sample_contig_categories=["Sample"]
            + [category for category in metadata if category in sample_contig_allowed],
            sample_mag_categories=["Sample"] + [category for category in metadata if category in sample_mag_allowed],
            sample_sources={category: info["source"] for category, info in metadata.items()},
        )
