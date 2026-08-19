"""Pure projections shared by the four Apply rendering modes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable


def selected_labels(group: Any) -> list[str]:
    if group is None:
        return []
    return [group.labels[index] for index in group.active]


def selected_features(groups: Iterable[Any]) -> list[str]:
    return [label for group in groups for label in selected_labels(group)]


def one_sample_features(genome_group: Any, variable_groups: Iterable[Any]) -> list[str]:
    return [*selected_labels(genome_group), *selected_features(variable_groups)]


def all_sample_features(genome_group: Any, variable_groups: Iterable[Any]) -> tuple[list[str], str | None]:
    variable = next(
        (labels[-1] for group in variable_groups if (labels := selected_labels(group))),
        None,
    )
    return selected_labels(genome_group), variable


@dataclass(frozen=True)
class TrackVisibility:
    gene_map: bool
    sequence: bool
    translation: bool


def track_visibility(
    window: float,
    *,
    max_genemap_window: int,
    max_sequence_window: int,
    sequence_requested: bool,
    translation_requested: bool,
) -> TrackVisibility:
    return TrackVisibility(
        gene_map=window <= max_genemap_window,
        sequence=sequence_requested and window <= max_sequence_window,
        translation=translation_requested and window <= max_sequence_window,
    )


def without_gene_map(features: Iterable[str], gene_map_visible: bool) -> list[str]:
    values = list(features)
    return values if gene_map_visible else [feature for feature in values if feature != "Gene map"]
