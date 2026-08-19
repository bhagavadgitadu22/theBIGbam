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


def _format_zoom(data, plot_type, start, end):
    starts = data["bin_start"]
    ends = data["bin_end"]
    mask = (ends >= start - 1) & (starts <= end - 1)
    starts, ends = starts[mask], ends[mask]
    if len(starts) == 0:
        return None
    midpoints = ((starts + ends) / 2.0 + 1).tolist()
    if plot_type == "bars":
        y = data.get("max", data.get("mean"))[mask].tolist()
        return {
            "x": midpoints,
            "y": y,
            "width": (ends - starts + 1).tolist(),
            "first_pos": (starts + 1).tolist(),
            "last_pos": (ends + 1).tolist(),
        }
    y = data["mean"][mask].tolist()
    if data.get("sparse") and midpoints:
        bin_size = int(ends[0] - starts[0] + 1)
        anchored_x, anchored_y = [midpoints[0] - bin_size], [0.0]
        for index, midpoint in enumerate(midpoints):
            if index and starts[index] > ends[index - 1] + 1:
                anchored_x.extend((midpoints[index - 1] + bin_size, midpoint - bin_size))
                anchored_y.extend((0.0, 0.0))
            anchored_x.append(midpoint)
            anchored_y.append(y[index])
        anchored_x.append(midpoints[-1] + bin_size)
        anchored_y.append(0.0)
        return {"x": anchored_x, "y": anchored_y}
    return {"x": midpoints, "y": y}


def _format_chunks(data, plot_type, start, end, position_scale=1, position_offset=0):
    x = data["x"] * position_scale + position_offset
    y = data["y"]
    left = max(0, int(np.searchsorted(x, start - 1, side="right")) - 1)
    right = min(len(x), int(np.searchsorted(x, end - 1, side="left")) + 1)
    if right <= left:
        return None
    original_length = len(x)
    sliced_x, sliced_y = x[left:right], y[left:right]
    result = {"x": (sliced_x + 1).tolist(), "y": sliced_y.tolist()}
    for key, value in data.items():
        if key in {"x", "y", "sparse"}:
            continue
        if hasattr(value, "__len__") and len(value) == original_length:
            result[key] = value[left:right].tolist() if hasattr(value, "tolist") else list(value[left:right])
    if plot_type == "bars":
        result.update(width=[1] * len(result["x"]), first_pos=result["x"], last_pos=result["x"])
    elif data.get("sparse") and len(sliced_x):
        metadata = {key: list(value) for key, value in result.items() if key not in {"x", "y"}}
        anchored_x, anchored_y = [], []
        anchored_metadata = {key: [] for key in metadata}

        def anchor(position):
            anchored_x.append(float(position) + 1)
            anchored_y.append(0.0)
            for values in anchored_metadata.values():
                values.append(None)

        for index, position in enumerate(sliced_x):
            if index == 0 and position > 0:
                anchor(position - 1)
            elif index and position > sliced_x[index - 1] + 1:
                anchor(sliced_x[index - 1] + 1)
                anchor(position - 1)
            anchored_x.append(float(position) + 1)
            anchored_y.append(float(sliced_y[index]))
            for key, values in metadata.items():
                anchored_metadata[key].append(values[index])
        anchor(sliced_x[-1] + 1)
        result = {"x": anchored_x, "y": anchored_y, **anchored_metadata}
    return result


def _apply_threshold(data, relative_threshold):
    if not data or relative_threshold <= 0 or not data.get("y"):
        return data
    maximum = max(data["y"])
    if maximum > 0:
        cutoff = relative_threshold * maximum
        data["y"] = [value if value >= cutoff else 0 for value in data["y"]]
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
        samples = self.repository.resolve_samples(contig_id, request.allowed_samples, request.ordering)
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
        use_zoom = region.length > (
            10_000_000 if style.variable_name in {"gc_content", "gc_skew"} else region.max_base_resolution
        )
        if use_zoom:
            retrieval_started = time.perf_counter()
            blob = self.repository.contig_zoom_blob(contig_id, feature_id)
            self._phase("blob_retrieval", retrieval_started)
            decoding_started = time.perf_counter()
            data = self._decode_zoom(blob, scale, zoom_bins, style.plot_type, region) if blob else None
            self._phase("blob_decoding", decoding_started)
        else:
            first = max(0, ((region.start - 1) // position_scale) // chunk_size)
            last = ((region.end - 1) // position_scale) // chunk_size
            retrieval_started = time.perf_counter()
            rows = self.repository.contig_chunks(contig_id, feature_id, first, last)
            self._phase("blob_retrieval", retrieval_started)
            decoder = decode_raw_sparse_chunks if style.encoding == "sparse" else decode_raw_chunks
            decoding_started = time.perf_counter()
            data = (
                _format_chunks(
                    decoder(rows, scale, chunk_size),
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
        return FeatureSeries(style, _as_tuple_mapping(data)) if data else None

    @staticmethod
    def _decode_zoom(blob, scale, bins, plot_type, region):
        for bin_size in bins:
            if region.length // bin_size <= 10_000:
                decoded = decode_zoom_standalone(blob, bin_size, scale, bins)
                if decoded is not None:
                    return _format_zoom(decoded, plot_type, region.start, region.end)
        decoded = decode_zoom_standalone(blob, bins[-1], scale, bins)
        return _format_zoom(decoded, plot_type, region.start, region.end) if decoded is not None else None
