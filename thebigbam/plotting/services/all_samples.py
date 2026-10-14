"""BLOB decoding and plain-data assembly for ALL SAMPLES plots.

This module deliberately has no Bokeh or Panel imports.
"""

from __future__ import annotations

import time
from collections.abc import Mapping

import numpy as np

from ...database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, decode_zoom_standalone
from ..models.plots import AllSamplesPlotData, AllSamplesPlotRequest, FeatureSeries, GenomeTrack, SampleTrack
from ..repositories.all_samples import AllSamplesRepository
from .features import apply_minimum_frequency
from .blob_features import FeatureDataService, _format_zoom, _slice_with_anchors


def _as_tuple_mapping(data: Mapping) -> dict[str, tuple]:
    result = {}
    for key, value in data.items():
        if key == "sparse":
            continue
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, list):
            value = tuple(value)
        elif not isinstance(value, tuple):
            continue
        result[key] = value
    return result


def _format_chunks(data, plot_type, start, end, position_scale=1, position_offset=0):
    data = dict(data)
    data["x"] = data["x"] * position_scale + position_offset
    return _slice_with_anchors(data, start, end, plot_type)


def _apply_threshold(data, relative_threshold):
    if not data or relative_threshold <= 0 or not data.get("y"):
        return data
    data["y"] = apply_minimum_frequency(data["y"], relative_threshold)
    return data


class AllSamplesDataService:
    def __init__(self, repository: AllSamplesRepository, diagnostics=None) -> None:
        self.repository = repository
        self.diagnostics = diagnostics
        self.phase_seconds: dict[str, float] = {}

    def _phase(self, name, started):
        self.phase_seconds[name] = self.phase_seconds.get(name, 0.0) + time.perf_counter() - started

    def load(self, request: AllSamplesPlotRequest) -> AllSamplesPlotData:
        cache_key = ("all_samples_plot_data", request)
        if request.data_cache is not None:
            hit, cached = request.data_cache.get(cache_key)
            if request.profiler is not None:
                request.profiler.cache("all_samples_data", hit)
            if hit:
                return cached
        started = time.perf_counter()
        contig_id, contig_name, contig_length = self.repository.resolve_contig(request.contig_name)
        samples = self._resolve_samples(contig_id, request)
        self._phase("resolve_subjects", started)

        sample_series: dict[int, list[FeatureSeries]] = {sid: [] for sid, _ in samples}
        if request.variable:
            styles = self.repository.feature_styles(request.variable, "Feature_blob")
            for style in styles:
                decoded = self._load_sample_style(contig_id, [sid for sid, _ in samples], style, request)
                for sid, series in decoded.items():
                    sample_series[sid].append(series)

        tracks = tuple(SampleTrack(sid, name, tuple(sample_series[sid])) for sid, name in samples if sample_series[sid])
        genome_tracks = []
        for feature in request.genome_features:
            styles = self.repository.feature_styles(feature, "Contig_blob")
            series = tuple(filter(None, (self._load_contig_style(contig_id, style, request) for style in styles)))
            if series:
                genome_tracks.append(GenomeTrack(feature, series))

        y_values = [value for track in tracks for series in track.series for value in series.data.get("y", ())]
        result = AllSamplesPlotData(
            contig_id,
            contig_name,
            contig_length,
            request.region,
            tracks,
            tuple(genome_tracks),
            min([0.0, *y_values]),
            max([0.0, *y_values]),
        )
        if request.data_cache is not None:
            request.data_cache.put(cache_key, result)
        return result

    def _resolve_samples(self, contig_id, request):
        samples = self.repository.available_samples(contig_id)
        if request.allowed_samples is not None:
            samples = [item for item in samples if item[1] in request.allowed_samples]
        if not samples:
            raise ValueError("No samples match the current contig and filters")
        if request.ordering.column:
            samples = self.repository.order_samples(contig_id, samples, request.ordering)
        if request.ordering.max_samples is not None:
            samples = samples[: request.ordering.max_samples]
        return samples

    def _load_sample_style(self, contig_id, sample_ids, style, request):
        feature_id = self.repository.feature_id(style.variable_name)
        if feature_id is None:
            return {}
        scale, chunk_size, zoom_bins, _ = self.repository.blob_settings(style.variable_name)
        region = request.region
        started = time.perf_counter()
        raw_by_sample = {}
        use_zoom = region.length > region.max_base_resolution
        if use_zoom:
            retrieval_started = time.perf_counter()
            blobs = self.repository.sample_zoom_blobs(contig_id, sample_ids, feature_id)
            self._phase("blob_retrieval", retrieval_started)
            decoding_started = time.perf_counter()
            for sid, blob in blobs.items():
                raw_by_sample[sid] = self._decode_zoom(blob, scale, zoom_bins, style.plot_type, region)
            self._phase("blob_decoding", decoding_started)
        else:
            first = max(0, (region.start - 1) // chunk_size)
            last = (region.end - 1) // chunk_size
            retrieval_started = time.perf_counter()
            chunks = self.repository.sample_chunks(contig_id, sample_ids, feature_id, first, last)
            self._phase("blob_retrieval", retrieval_started)
            decoder = decode_raw_sparse_chunks if style.encoding == "sparse" else decode_raw_chunks
            decoding_started = time.perf_counter()
            for sid, rows in chunks.items():
                raw_by_sample[sid] = _format_chunks(
                    decoder(rows, scale, chunk_size), style.plot_type, region.start, region.end
                )
            self._phase("blob_decoding", decoding_started)
        self._phase("retrieve_and_decode", started)
        return {
            sid: FeatureSeries(style, _as_tuple_mapping(_apply_threshold(data, request.display.min_relative_value)))
            for sid, data in raw_by_sample.items()
            if data
        }

    def _load_contig_style(self, contig_id, style, request):
        feature_id = self.repository.feature_id(style.variable_name)
        if feature_id is None:
            return None
        scale, chunk_size, zoom_bins, position_scale = self.repository.blob_settings(style.variable_name)
        region = request.region
        intervals = style.variable_name in {"direct_repeat_count", "inverted_repeat_count", "direct_repeat_identity", "inverted_repeat_identity"} and self.repository.storage_settings(style.variable_name).get("intervals",False)
        use_zoom = region.length > (
            10_000_000 if style.variable_name in {"gc_content", "gc_skew"} else region.max_base_resolution
        )
        if use_zoom:
            retrieval_started = time.perf_counter()
            blob = self.repository.contig_zoom_blob(contig_id, feature_id)
            self._phase("blob_retrieval", retrieval_started)
            decoding_started = time.perf_counter()
            data = self._decode_zoom(blob, scale, zoom_bins, style.plot_type, region, intervals) if blob else None
            self._phase("blob_decoding", decoding_started)
        else:
            first = max(0, ((region.start - 1) // position_scale) // chunk_size)
            last = ((region.end - 1) // position_scale) // chunk_size
            retrieval_started = time.perf_counter()
            rows = self.repository.contig_chunks(contig_id, feature_id, first, last)
            self._phase("blob_retrieval", retrieval_started)
            decoder = decode_raw_sparse_chunks if style.encoding == "sparse" else decode_raw_chunks
            decoding_started = time.perf_counter()
            decoded = decoder(rows, scale, chunk_size) if rows else None
            if decoded and style.variable_name in {"direct_repeat_count", "inverted_repeat_count", "direct_repeat_identity", "inverted_repeat_identity"}:
                settings = self.repository.storage_settings(style.variable_name)
                if settings.get("intervals"):
                    decoded["sparse"] = False
                    decoded["intervals"] = True
                    FeatureDataService._repeat_labels(decoded)
            data = (
                _format_chunks(
                    decoded,
                    style.plot_type,
                    region.start,
                    region.end,
                    position_scale,
                    position_scale // 2 - 1 if position_scale > 1 else 0,
                )
                if rows
                else None
            )
            self._phase("blob_decoding", decoding_started)
        if data:
            FeatureDataService(self.repository)._resolve_partners(
                data, style.variable_name, ((contig_id, 0, region.end),))
        return FeatureSeries(style, _as_tuple_mapping(data)) if data else None

    @staticmethod
    def _decode_zoom(blob, scale, bins, plot_type, region, intervals=False):
        for bin_size in bins:
            if region.length // bin_size <= 10_000:
                decoded = decode_zoom_standalone(blob, bin_size, scale, bins)
                if decoded is not None:
                    return _format_zoom(decoded, plot_type, region.start, region.end, intervals)
        decoded = decode_zoom_standalone(blob, bins[-1], scale, bins)
        return _format_zoom(decoded, plot_type, region.start, region.end, intervals) if decoded is not None else None
