"""Prepare feature-series data for typed plot composition requests."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any

from .composition_data import (
    get_feature_data,
    get_mag_feature_data,
    get_variable_metadata_batch,
    parse_requested_features,
    split_contig_vs_sample_features,
)


@dataclass(frozen=True)
class PreparedFeature:
    name: str
    series: list[dict[str, Any]]
    sample_name: str | None = None


@dataclass(frozen=True)
class PreparedFeatures:
    contig: tuple[PreparedFeature, ...]
    sample: tuple[PreparedFeature, ...]


class CompositionDataService:
    def __init__(self, cursor: Any, cache=None, profiler=None) -> None:
        self.cursor = cursor
        self.cache = cache
        self.profiler = profiler

    def _cached(self, key, load):
        if self.cache is None:
            return load()
        hit, value = self.cache.get(key)
        if self.profiler is not None:
            self.profiler.cache("feature_data", hit)
        if hit:
            return deepcopy(value)
        value = load()
        self.cache.put(key, deepcopy(value))
        return value

    def for_contig(self, request, contig_id: int, sample_id: int | None) -> PreparedFeatures:
        requested, metadata, contig_features, sample_features = self._classify(request)
        del requested
        contig = self._load_contig_features(request, contig_features, metadata, contig_id, None)
        sample = (
            self._load_contig_features(request, sample_features, metadata, contig_id, sample_id)
            if sample_id is not None
            else []
        )
        return PreparedFeatures(tuple(contig), tuple(sample))

    def for_mag(self, request, context, repository, sample_id: int | None) -> PreparedFeatures:
        _requested, metadata, contig_features, sample_features = self._classify(request)
        contig = self._load_mag_features(request, context, contig_features, metadata, [(None, None)])
        if request.is_all:
            samples = repository.ordered_mag_samples(request, context.mag_id)
        elif sample_id is not None:
            samples = [(sample_id, request.sample_name)]
        else:
            samples = []
        sample = self._load_mag_features(request, context, sample_features, metadata, samples)
        return PreparedFeatures(tuple(contig), tuple(sample))

    def primary_maximum(self, request, contig_id: int, sample_id: int) -> float:
        try:
            series = get_feature_data(
                self.cursor,
                "Primary alignments",
                contig_id,
                sample_id,
                request.xstart,
                request.xend,
                enable_timing=request.enable_timing,
                encoding_by_feature=request.encoding_by_feature,
            )
            return max((value for item in series for value in item.get("y", [])), default=0)
        except Exception:
            return 0

    def _classify(self, request):
        requested = parse_requested_features(request.tracks.features)
        metadata = get_variable_metadata_batch(self.cursor, requested)
        contig, sample = split_contig_vs_sample_features(metadata, requested)
        return requested, metadata, contig, sample

    def _load_contig_features(self, request, features, metadata, contig_id, sample_id):
        prepared = []
        for feature in features:
            try:
                key = (
                    "contig_feature",
                    feature,
                    contig_id,
                    sample_id,
                    request.xstart,
                    request.xend,
                    request.display.max_base_resolution,
                    request.display.min_relative_value,
                )
                series = self._cached(
                    key,
                    lambda: get_feature_data(
                        self.cursor,
                        feature,
                        contig_id,
                        sample_id,
                        request.xstart,
                        request.xend,
                        variable_metadata=metadata.get(feature),
                        max_base_resolution=request.display.max_base_resolution,
                        min_relative_value=request.display.min_relative_value,
                        enable_timing=request.enable_timing,
                        encoding_by_feature=request.encoding_by_feature,
                        profiler=self.profiler,
                    ),
                )
                prepared.append(PreparedFeature(feature, series))
            except Exception as error:
                print(f"Error processing feature '{feature}': {error}", flush=True)
        return prepared

    def _load_mag_features(self, request, context, features, metadata, samples):
        prepared = []
        for feature in features:
            for sample_id, sample_name in samples:
                try:
                    key = (
                        "mag_feature",
                        feature,
                        context.mag_id,
                        sample_id,
                        context.xstart,
                        context.xend,
                        request.display.max_base_resolution,
                        request.display.min_relative_value,
                        context.feature_members,
                    )
                    series = self._cached(
                        key,
                        lambda: get_mag_feature_data(
                            self.cursor,
                            feature,
                            context.mag_id,
                            sample_id,
                            context.total_length,
                            xstart=context.xstart,
                            xend=context.xend,
                            variable_metadata=metadata.get(feature),
                            max_base_resolution=request.display.max_base_resolution,
                            min_relative_value=request.display.min_relative_value,
                            members=context.feature_members,
                            encoding_by_feature=request.encoding_by_feature,
                            profiler=self.profiler,
                        ),
                    )
                    if sample_name is not None:
                        for item in series:
                            item["title"] = sample_name
                    prepared.append(PreparedFeature(feature, series, sample_name))
                except Exception as error:
                    suffix = f" for sample '{sample_name}'" if sample_name else ""
                    print(f"Error processing MAG feature '{feature}'{suffix}: {error}", flush=True)
        return prepared
