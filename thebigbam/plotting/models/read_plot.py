"""Immutable inputs and plain data for BAM read plots."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from ..shared.defaults import (
    DEFAULT_INDIVIDUAL_READ_HEIGHT,
    DEFAULT_READ_DISCREPANCY_WINDOW,
    DEFAULT_READ_PLOT_HEIGHT,
)
from .filters import FilterExpression
from .read_events import ReadEvent


@dataclass(frozen=True)
class ReadColorRule:
    field: str
    operator: str
    value: Any
    color: str


@dataclass(frozen=True)
class ReadPlotOptions:
    enabled: bool = False
    max_reads: int = 1000
    link_pairs: bool = False
    plot_height: int = DEFAULT_READ_PLOT_HEIGHT
    read_height: int = DEFAULT_INDIVIDUAL_READ_HEIGHT
    event_variables: frozenset[str] = frozenset()
    filters: FilterExpression = FilterExpression(())
    max_discrepancy_window: int = DEFAULT_READ_DISCREPANCY_WINDOW
    color_rules: tuple[ReadColorRule, ...] = ()


@dataclass(frozen=True)
class ReadPlotRequest:
    bam_path: str
    contig: str
    contig_length: int
    start: int
    end: int
    max_reads: int
    link_pairs: bool
    classes: frozenset[str]
    colors: Mapping[str, str] = field(default_factory=dict)
    plot_height: int = DEFAULT_READ_PLOT_HEIGHT
    read_height: int = DEFAULT_INDIVIDUAL_READ_HEIGHT
    members: tuple[tuple[str, int, int], ...] = ()
    sample_name: str | None = None
    profiler: Any = field(default=None, compare=False, repr=False)
    event_variables: frozenset[str] = frozenset()
    reference_source: Any = field(default=None, compare=False, repr=False)
    circular: bool = False
    filters: FilterExpression = FilterExpression(())
    max_discrepancy_window: int = DEFAULT_READ_DISCREPANCY_WINDOW
    color_rules: tuple[ReadColorRule, ...] = ()


@dataclass(frozen=True)
class ReadTrackConfig:
    options: ReadPlotOptions
    catalog: Any
    classes: frozenset[str]
    colors: Mapping[str, str]
    diagnose: bool = False
    reference_source: Any = field(default=None, compare=False, repr=False)


@dataclass(frozen=True)
class ReadAlignment:
    name: str
    start: int
    end: int
    blocks: tuple[tuple[int, int], ...]
    forward: bool
    alignment_class: str
    mapq: int
    cigar: str
    mate_start: int | None
    is_read1: bool
    is_paired: bool
    events: tuple[ReadEvent, ...] = ()
    alignment_id: tuple | None = field(default=None, compare=False, repr=False)
    color: str | None = None


@dataclass(frozen=True)
class ReadPlotData:
    alignments: tuple[ReadAlignment, ...]
    fetched: int
    message: str | None = None
    warnings: tuple[str, ...] = ()
    details_loaded: bool = True
