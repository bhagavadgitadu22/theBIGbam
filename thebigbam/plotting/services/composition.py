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
    def __init__(self, feature_repository: Any, composition_repository=None, cache=None, profiler=None) -> None:
        self.feature_repository = feature_repository
        self.composition_repository = composition_repository
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

    def for_mag(self, request, context, sample_id: int | None) -> PreparedFeatures:
        _requested, metadata, contig_features, sample_features = self._classify(request)
        contig = self._load_mag_features(request, context, contig_features, metadata, [(None, None)])
        if request.is_all:
            if self.composition_repository is None:
                raise ValueError("A composition repository is required for all-sample MAG plots")
            samples = self.ordered_mag_samples(request, context.mag_id)
        elif sample_id is not None:
            samples = [(sample_id, request.sample_name)]
        else:
            samples = []
        sample = self._load_mag_features(request, context, sample_features, metadata, samples)
        return PreparedFeatures(tuple(contig), tuple(sample))

    def ordered_mag_samples(self, request, mag_id):
        """Apply user filtering, ordering, and limits to repository sample rows."""
        membership_key = ("mag_samples", mag_id)
        if self.cache is None:
            rows = self.composition_repository.mag_samples(mag_id)
        else:
            hit, cached_rows = self.cache.get(membership_key)
            if self.profiler is not None:
                self.profiler.cache("mag_samples", hit)
            if hit:
                rows = list(cached_rows)
            else:
                rows = self.composition_repository.mag_samples(mag_id)
                self.cache.put(membership_key, tuple(rows))
        if request.allowed_samples is not None:
            rows = [(sid, name) for sid, name in rows if name in request.allowed_samples]
        ordering = request.sample_ordering
        if ordering.source and ordering.metric and ordering.metric != "Sample_name" and rows:
            try:
                values_key = (
                    "mag_sample_order_values",
                    mag_id,
                    ordering.source,
                    ordering.metric,
                )
                if self.cache is None:
                    values = self.composition_repository.mag_sample_order_values(rows, request.mag_name, ordering)
                else:
                    hit, values = self.cache.get(values_key)
                    if self.profiler is not None:
                        self.profiler.cache("mag_sample_order", hit)
                    if not hit:
                        # Load against complete membership so filter changes can
                        # reuse the same order-independent metric snapshot.
                        _, all_rows = self.cache.get(membership_key)
                        values = self.composition_repository.mag_sample_order_values(
                            list(all_rows), request.mag_name, ordering
                        )
                        self.cache.put(values_key, values)
                present = [row for row in rows if values.get(row[0]) is not None]
                missing = [row for row in rows if values.get(row[0]) is None]
                present.sort(key=lambda row: values[row[0]], reverse=not ordering.ascending)
                rows = present + missing
            except Exception as error:
                print(f"Warning: Could not order samples by '{ordering.metric}': {error}", flush=True)
                rows.sort(key=lambda row: row[1], reverse=not ordering.ascending)
        else:
            rows.sort(key=lambda row: row[1], reverse=not ordering.ascending)
        if request.max_samples is not None and len(rows) > request.max_samples:
            print(
                f"Plotting {request.max_samples}/{len(rows)} samples (limited by 'Max number of samples plotted')",
                flush=True,
            )
            rows = rows[: request.max_samples]
        return rows

    def primary_maximum(self, request, contig_id: int, sample_id: int) -> float:
        try:
            series = get_feature_data(
                self.feature_repository,
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
        metadata = get_variable_metadata_batch(self.feature_repository, requested)
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
                        self.feature_repository,
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
                            self.feature_repository,
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
