"""Plain data contracts shared by plotting repositories, services and renderers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class GenomicRegion:
    start: int
    end: int
    max_base_resolution: int = 10_000

    def __post_init__(self) -> None:
        if self.start < 1 or self.end <= self.start:
            raise ValueError("Region must use 1-based coordinates with end > start")
        if self.max_base_resolution < 1:
            raise ValueError("max_base_resolution must be positive")

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class SampleOrdering:
    source: str | None = None
    column: str | None = None
    ascending: bool = True
    max_samples: int | None = None

    def __post_init__(self) -> None:
        if self.max_samples is not None and self.max_samples < 1:
            raise ValueError("max_samples must be positive")


@dataclass(frozen=True)
class DisplayOptions:
    subplot_height: int = 130
    same_y_scale: bool = False
    min_relative_value: float = 0.0

    def __post_init__(self) -> None:
        if self.subplot_height < 1:
            raise ValueError("subplot_height must be positive")
        if not 0.0 <= self.min_relative_value <= 1.0:
            raise ValueError("min_relative_value must be between 0 and 1")


@dataclass(frozen=True)
class AllSamplesPlotRequest:
    contig_name: str
    variable: str | None
    region: GenomicRegion
    allowed_samples: frozenset[str] | None = None
    ordering: SampleOrdering = field(default_factory=SampleOrdering)
    display: DisplayOptions = field(default_factory=DisplayOptions)
    genome_features: tuple[str, ...] = ()
    data_cache: Any = field(default=None, compare=False, repr=False)
    profiler: Any = field(default=None, compare=False, repr=False)


@dataclass(frozen=True)
class FeatureStyle:
    plot_type: str
    color: str
    alpha: float
    fill_alpha: float
    size: float
    title: str
    variable_name: str
    feature_table: str
    encoding: str | None = None


@dataclass(frozen=True)
class FeatureSeries:
    style: FeatureStyle
    data: Mapping[str, tuple[Any, ...]]
    is_relative_scaled: bool = False


@dataclass(frozen=True)
class SampleTrack:
    sample_id: int
    sample_name: str
    series: tuple[FeatureSeries, ...]


@dataclass(frozen=True)
class GenomeTrack:
    feature_name: str
    series: tuple[FeatureSeries, ...]


@dataclass(frozen=True)
class AllSamplesPlotData:
    contig_id: int
    contig_name: str
    contig_length: int
    region: GenomicRegion
    sample_tracks: tuple[SampleTrack, ...]
    genome_tracks: tuple[GenomeTrack, ...] = ()
    y_min: float = 0.0
    y_max: float = 0.0
