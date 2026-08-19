"""Immutable filter-expression contracts independent of Panel and DuckDB."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class FilterPredicate:
    category: str
    column: str
    operator: str
    value: Any


@dataclass(frozen=True)
class FilterSection:
    predicates: tuple[FilterPredicate, ...]
    connectors: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if len(self.connectors) != max(0, len(self.predicates) - 1):
            raise ValueError("A section needs one connector between adjacent predicates")


@dataclass(frozen=True)
class FilterExpression:
    sections: tuple[FilterSection, ...]
    connectors: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if len(self.connectors) != max(0, len(self.sections) - 1):
            raise ValueError("An expression needs one connector between adjacent sections")

    @property
    def active(self) -> bool:
        return any(section.predicates for section in self.sections)


@dataclass(frozen=True)
class CompiledFilter:
    sql: str
    parameters: tuple[Any, ...]


@dataclass(frozen=True)
class FilterResult:
    compiled: CompiledFilter
    pair_count: int
    contig_count: int
    sample_count: int
    mag_pair_count: int | None = None
