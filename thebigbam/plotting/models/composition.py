"""Immutable requests used at the final plot-composition boundary."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class TrackSelection:
    features: tuple[str, ...] = ()
    feature_types: tuple[str, ...] = ()
    plot_isoforms: bool = True
    plot_sequence: bool = False
    plot_translation: bool = False
    label_key: str | None = None
    color_rules: tuple[Mapping[str, Any], ...] = ()


@dataclass(frozen=True)
class CompositionDisplay:
    subplot_height: int = 100
    genemap_height: int | None = None
    sequence_height: int | None = None
    translation_height: int | None = None
    same_y_scale: bool = False
    max_base_resolution: int | None = None
    max_genemap_window: int | None = None
    max_sequence_window: int | None = None
    min_relative_value: float = 0.0


@dataclass(frozen=True)
class SingleSampleCompositionRequest:
    contig_name: str
    sample_name: str | None
    xstart: int | None = None
    xend: int | None = None
    genbank_path: str | None = None
    tracks: TrackSelection = field(default_factory=TrackSelection)
    display: CompositionDisplay = field(default_factory=CompositionDisplay)
    mag_name: str | None = None
    enable_timing: bool = False
    encoding_by_feature: Mapping[str, Any] | None = None
    data_cache: Any = field(default=None, compare=False, repr=False)
    profiler: Any = field(default=None, compare=False, repr=False)


@dataclass(frozen=True)
class MagOrdering:
    source: str | None = None
    metric: str | None = None
    ascending: bool = True
    sample_name: str | None = None


@dataclass(frozen=True)
class MagCompositionRequest:
    mag_name: str
    sample_name: str | None
    xstart: int | None = None
    xend: int | None = None
    genbank_path: str | None = None
    tracks: TrackSelection = field(default_factory=TrackSelection)
    display: CompositionDisplay = field(default_factory=CompositionDisplay)
    mag_track_colors: tuple[Mapping[str, Any], ...] = ()
    max_track_dots: int | None = 1000
    is_all: bool = False
    allowed_samples: frozenset[str] | None = None
    max_samples: int | None = None
    ordering: MagOrdering = field(default_factory=MagOrdering)
    sample_ordering: MagOrdering = field(default_factory=MagOrdering)
    focus_contig: str | None = None
    enable_timing: bool = False
    encoding_by_feature: Mapping[str, Any] | None = None
    data_cache: Any = field(default=None, compare=False, repr=False)
    profiler: Any = field(default=None, compare=False, repr=False)
