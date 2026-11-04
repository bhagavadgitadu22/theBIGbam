"""Immutable inputs and plain data for BAM read plots."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from ..shared.defaults import DEFAULT_INDIVIDUAL_READ_HEIGHT, DEFAULT_READ_PLOT_HEIGHT


@dataclass(frozen=True)
class ReadPlotOptions:
    enabled: bool = False
    max_reads: int = 1000
    link_pairs: bool = False
    selected_subplots: tuple[str, ...] = ()
    plot_height: int = DEFAULT_READ_PLOT_HEIGHT
    read_height: int = DEFAULT_INDIVIDUAL_READ_HEIGHT


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


@dataclass(frozen=True)
class ReadTrackConfig:
    options: ReadPlotOptions
    catalog: Any
    classes: frozenset[str]
    colors: Mapping[str, str]
    diagnose: bool = False


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


@dataclass(frozen=True)
class ReadPlotData:
    alignments: tuple[ReadAlignment, ...]
    fetched: int
    message: str | None = None
