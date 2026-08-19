"""Immutable requests and results for BLOB-backed feature loading."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class FeatureRegion:
    start: int
    end: int

    def __post_init__(self):
        if self.start < 1 or self.end < self.start:
            raise ValueError("Feature regions use inclusive, 1-based coordinates")

    @property
    def length(self):
        return self.end - self.start + 1


@dataclass(frozen=True)
class FeatureLoadRequest:
    feature: str
    contig_id: int
    sample_id: int | None
    region: FeatureRegion
    max_base_resolution: int = 10_000
    minimum_relative_value: float = 0.0
    encoding_by_feature: Mapping[str, str] | None = None


@dataclass(frozen=True)
class MagMember:
    contig_id: int
    length: int
    offset: int


@dataclass(frozen=True)
class MagFeatureLoadRequest:
    feature: str
    sample_id: int | None
    region: FeatureRegion
    members: tuple[MagMember, ...]
    max_base_resolution: int = 10_000
    minimum_relative_value: float = 0.0
    encoding_by_feature: Mapping[str, str] | None = None


@dataclass(frozen=True)
class DecodedFeatureData:
    x: tuple[float, ...]
    y: tuple[float, ...]
    metadata: Mapping[str, tuple[Any, ...]]
    sparse: bool = False
