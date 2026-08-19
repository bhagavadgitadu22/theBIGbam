"""Immutable data loaded before construction of a plotting session."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterator, Mapping


@dataclass(frozen=True)
class PreloadedPlotData(Mapping[str, Any]):
    annotation_types: tuple[str, ...]
    has_mags: bool
    mag_to_contigs: Mapping[str, list[str]]
    contig_to_mag: Mapping[str, str]
    mags: tuple[str, ...]
    contigs: tuple[str, ...]
    contig_lengths: Mapping[str, int]
    contig_name_to_id: Mapping[str, int]
    contig_id_to_name: Mapping[int, str]
    sample_name_to_id: Mapping[str, int]
    sample_id_to_name: Mapping[int, str]
    samples: tuple[str, ...]
    has_samples: bool
    mag_to_sample_ids: Mapping[str, set[int]]
    module_names: tuple[str, ...]
    module_variables: tuple[tuple[str, ...], ...]
    module_helps: tuple[str, ...]
    custom_contig_subplots: tuple[str, ...]
    filtering_metadata: Mapping[str, Any]
    subplot_to_varnames: Mapping[str, list[str]]
    encoding_by_feature: Mapping[str, str]
    total_coverage_count: int
    filter_encode: Mapping[str, Any]
    mag_to_contig_offsets: dict[str, dict[str, int]] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        try:
            return getattr(self, key)
        except AttributeError as error:
            raise KeyError(key) from error

    def __iter__(self) -> Iterator[str]:
        return iter(self.__dataclass_fields__)

    def __len__(self) -> int:
        return len(self.__dataclass_fields__)
