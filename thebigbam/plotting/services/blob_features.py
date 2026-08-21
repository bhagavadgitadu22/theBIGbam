"""Decode and transform raw feature BLOBs without SQL or UI dependencies."""

from __future__ import annotations

import time

import numpy as np

from thebigbam.database.blob_decoder import (
    decode_raw_chunks,
    decode_raw_sparse_chunks,
    decode_zoom_standalone,
    is_sparse_zoom_blob,
)

from ..models.blobs import FeatureLoadRequest, MagFeatureLoadRequest
from .features import apply_minimum_frequency


def _slice_with_anchors(data, start, end, plot_type):
    x, y = np.asarray(data["x"]), np.asarray(data["y"])
    if len(x) == 0:
        return None
    left = max(0, int(np.searchsorted(x, start - 1, side="right")) - 1)
    right = min(len(x), int(np.searchsorted(x, end - 1, side="left")) + 1)
    original_length = len(x)
    x, y = x[left:right], y[left:right]
    metadata = {}
    for key, values in data.items():
        if key in {"x", "y", "sparse"} or not hasattr(values, "__len__") or len(values) != original_length:
            continue
        metadata[key] = list(values[left:right])
    if data.get("sparse") and plot_type != "bars" and len(x):
        new_x, new_y = [], []
        new_meta = {key: [] for key in metadata}

        def anchor(position):
            new_x.append(position)
            new_y.append(0.0)
            for values in new_meta.values():
                values.append(None)

        for index, position in enumerate(x):
            if index == 0 and position > 0:
                anchor(position - 1)
            elif index and position > x[index - 1] + 1:
                same_plateau = float(y[index]) == float(y[index - 1])
                previous_dense = index >= 2 and x[index - 1] == x[index - 2] + 1
                next_dense = index + 1 < len(x) and x[index + 1] == position + 1
                if not same_plateau or previous_dense or next_dense:
                    anchor(x[index - 1] + 1)
                    anchor(position - 1)
            new_x.append(position)
            new_y.append(float(y[index]))
            for key, values in metadata.items():
                new_meta[key].append(values[index])
        anchor(x[-1] + 1)
        x, y, metadata = np.asarray(new_x), np.asarray(new_y), new_meta
    result = {"x": (x + 1).tolist(), "y": y.tolist(), **metadata}
    if plot_type == "bars":
        result.update(width=[1] * len(x), first_pos=result["x"], last_pos=result["x"])
    return result


def _format_zoom(data, plot_type, start, end):
    starts, ends = np.asarray(data["bin_start"]), np.asarray(data["bin_end"])
    mask = (ends >= start - 1) & (starts <= end - 1)
    starts, ends = starts[mask], ends[mask]
    if len(starts) == 0:
        return None
    midpoints = ((starts + ends) / 2 + 1).tolist()
    if plot_type == "bars":
        values = np.asarray(data.get("max", data.get("mean", ())))[mask].tolist()
        return {
            "x": midpoints,
            "y": values,
            "width": (ends - starts + 1).tolist(),
            "first_pos": (starts + 1).tolist(),
            "last_pos": (ends + 1).tolist(),
        }
    values = np.asarray(data.get("mean", ()))[mask].tolist()
    if data.get("sparse") and midpoints:
        width, new_x, new_y = int(ends[0] - starts[0] + 1), [], []
        for index, midpoint in enumerate(midpoints):
            if index == 0:
                new_x.append(midpoint - width)
                new_y.append(0.0)
            elif starts[index] > ends[index - 1] + 1:
                new_x.extend((midpoints[index - 1] + width, midpoint - width))
                new_y.extend((0.0, 0.0))
            new_x.append(midpoint)
            new_y.append(values[index])
        new_x.append(midpoints[-1] + width)
        new_y.append(0.0)
        midpoints, values = new_x, new_y
    return {"x": midpoints, "y": values}


def _merge(parts):
    if not parts:
        return None
    keys = set.intersection(*(set(part) for part in parts))
    result = {}
    for key in keys:
        if key == "sparse":
            result[key] = any(bool(part[key]) for part in parts)
        elif hasattr(parts[0][key], "__len__"):
            result[key] = np.concatenate([np.asarray(part[key]) for part in parts])
    order = np.argsort(result["x"])
    for key, values in tuple(result.items()):
        if key != "sparse" and hasattr(values, "__len__") and len(values) == len(order):
            result[key] = values[order]
    return result


def _overlapping_members(members, start, end):
    """Return members intersecting a 1-based half-open MAG display region."""
    region_start, region_end = start - 1, end - 1
    return tuple(
        member
        for member in members
        if member.offset <= region_end and member.offset + member.length - 1 >= region_start
    )


class FeatureDataService:
    def __init__(self, repository):
        self.repository = repository
        self.phase_seconds = {}
        self.decoded_points = 0
        self.mag_members_total = 0
        self.mag_members_loaded = 0

    def _timed(self, name, function, *args, **kwargs):
        started = time.perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            self.phase_seconds[name] = self.phase_seconds.get(name, 0.0) + time.perf_counter() - started

    def _style(self, row, values):
        plot_type, color, alpha, fill_alpha, size, title, _table = row
        result = {
            "type": plot_type,
            "color": color,
            "alpha": alpha,
            "fill_alpha": fill_alpha,
            "size": size,
            "title": title,
            "is_relative_scaled": False,
            "has_stats": "mean" in values,
            "has_sequences": "sequence" in values,
        }
        result.update(values)
        self.decoded_points += len(values.get("x", ()))
        return result

    def load_contig(self, request: FeatureLoadRequest, metadata=None):
        started = time.perf_counter()
        rows = metadata if metadata is not None else self.repository.metadata(request.feature)
        output = []
        for row in rows:
            plot_type, _color, _alpha, _fill, _size, title, table = row
            if table not in {"Feature_blob", "Contig_blob"}:
                continue
            variable = self.repository.variable_name(title, table)
            feature_id = self.repository.feature_id(variable) if variable else None
            if feature_id is None:
                continue
            settings = self.repository.storage_settings(variable)
            contig_level = table == "Contig_blob"
            sample_ids = () if contig_level else (request.sample_id,)
            key = (request.contig_id, None if contig_level else request.sample_id)
            use_zoom = request.region.length - 1 > (
                10_000_000 if variable in {"gc_content", "gc_skew"} else request.max_base_resolution
            )
            zooms = self._timed(
                "blob_retrieval",
                self.repository.zoom_blobs,
                (request.contig_id,),
                sample_ids,
                feature_id,
                contig_level,
            )
            zoom = zooms.get(key)
            if use_zoom:
                values = self._timed(
                    "blob_decoding",
                    self._decode_zoom,
                    zoom,
                    settings,
                    plot_type,
                    request.region.start,
                    request.region.end,
                )
            else:
                window_size = settings["window_size"]
                first = max(0, ((request.region.start - 1) // window_size) // settings["chunk_size"])
                last = ((request.region.end - 1) // window_size) // settings["chunk_size"]
                chunks = self._timed(
                    "blob_retrieval",
                    self.repository.chunks,
                    (request.contig_id,),
                    sample_ids,
                    feature_id,
                    first,
                    last,
                    contig_level,
                ).get(key)
                sparse = request.encoding_by_feature and request.encoding_by_feature.get(variable) == "sparse"
                if sparse is None:
                    sparse = bool(zoom and is_sparse_zoom_blob(zoom))
                values = self._timed(
                    "blob_decoding",
                    self._decode_chunks,
                    chunks,
                    settings,
                    sparse,
                    plot_type,
                    request.region.start,
                    request.region.end,
                    window_size,
                )
            if values:
                transformed = time.perf_counter()
                values["y"] = apply_minimum_frequency(values["y"], request.minimum_relative_value)
                self._resolve_partners(values, variable)
                output.append(self._style(row, values))
                self.phase_seconds["transformation"] = self.phase_seconds.get("transformation", 0) + (
                    time.perf_counter() - transformed
                )
        self.phase_seconds["contig"] = self.phase_seconds.get("contig", 0) + time.perf_counter() - started
        return output

    def load_mag(self, request: MagFeatureLoadRequest, metadata=None):
        started = time.perf_counter()
        rows = metadata if metadata is not None else self.repository.metadata(request.feature)
        output, contig_ids = [], tuple(member.contig_id for member in request.members)
        self.mag_members_total = len(request.members)
        for row in rows:
            plot_type, _color, _alpha, _fill, _size, title, table = row
            if table not in {"Feature_blob", "Contig_blob"}:
                continue
            variable = self.repository.variable_name(title, table)
            feature_id = self.repository.feature_id(variable) if variable else None
            if feature_id is None:
                continue
            settings = self.repository.storage_settings(variable)
            contig_level = table == "Contig_blob"
            sample_ids = () if contig_level else (request.sample_id,)
            zooms = self._timed(
                "blob_retrieval", self.repository.zoom_blobs, contig_ids, sample_ids, feature_id, contig_level
            )
            sparse = request.encoding_by_feature and request.encoding_by_feature.get(variable) == "sparse"
            if sparse is None:
                sparse = any(is_sparse_zoom_blob(blob) for blob in zooms.values() if blob)
            use_zoom = request.region.length - 1 > request.max_base_resolution
            if use_zoom:
                self.mag_members_loaded = len(request.members)
                parts = self._timed(
                    "blob_decoding", self._mag_zoom_parts, request, zooms, settings, plot_type, contig_level
                )
                values = self._timed(
                    "transformation",
                    self._combine_zoom,
                    parts,
                    plot_type,
                    request.region.start,
                    request.region.end,
                )
            else:
                active_members = _overlapping_members(request.members, request.region.start, request.region.end)
                self.mag_members_loaded = len(active_members)
                if not active_members:
                    continue
                active_contig_ids = tuple(member.contig_id for member in active_members)
                max_local = max((member.length for member in active_members), default=0)
                chunks = self._timed(
                    "blob_retrieval",
                    self.repository.chunks,
                    active_contig_ids,
                    sample_ids,
                    feature_id,
                    0,
                    max_local // settings["chunk_size"],
                    contig_level,
                )
                parts = []
                decoding_started = time.perf_counter()
                for member in active_members:
                    key = (member.contig_id, None if contig_level else request.sample_id)
                    decoded = self._raw_chunks(chunks.get(key), settings, sparse)
                    if decoded and len(decoded["x"]):
                        decoded["x"] = decoded["x"] * settings["window_size"] + member.offset
                        parts.append(decoded)
                    else:
                        parts.append(
                            {
                                "x": np.asarray([member.offset, member.offset + member.length - 1]),
                                "y": np.asarray([0.0, 0.0]),
                                "sparse": False,
                            }
                        )
                self.phase_seconds["blob_decoding"] = self.phase_seconds.get("blob_decoding", 0) + (
                    time.perf_counter() - decoding_started
                )
                transformation_started = time.perf_counter()
                values = (
                    _slice_with_anchors(_merge(parts), request.region.start, request.region.end, plot_type)
                    if parts
                    else None
                )
                if values and settings["window_size"] > 1:
                    values["x"] = [position + settings["window_size"] // 2 - 1 for position in values["x"]]
                self.phase_seconds["transformation"] = self.phase_seconds.get("transformation", 0) + (
                    time.perf_counter() - transformation_started
                )
            if values:
                transformed = time.perf_counter()
                values["y"] = apply_minimum_frequency(values["y"], request.minimum_relative_value)
                self._resolve_partners(values, variable)
                output.append(self._style(row, values))
                self.phase_seconds["transformation"] = self.phase_seconds.get("transformation", 0) + (
                    time.perf_counter() - transformed
                )
        self.phase_seconds["mag"] = self.phase_seconds.get("mag", 0) + time.perf_counter() - started
        return output

    @staticmethod
    def _raw_chunks(chunks, settings, sparse):
        if not chunks:
            return None
        decoder = decode_raw_sparse_chunks if sparse else decode_raw_chunks
        return decoder(chunks, settings["scale"], settings["chunk_size"])

    def _decode_chunks(self, chunks, settings, sparse, plot_type, start, end, window_size=1):
        decoded = self._raw_chunks(chunks, settings, sparse)
        if not decoded:
            return None
        if window_size > 1:
            decoded["x"] = decoded["x"] * window_size
        values = _slice_with_anchors(decoded, start, end, plot_type)
        if values and window_size > 1:
            values["x"] = [position + window_size // 2 - 1 for position in values["x"]]
        return values

    @staticmethod
    def _decode_zoom(blob, settings, plot_type, start, end):
        if not blob:
            return None
        window = end - start
        target = next((size for size in settings["zoom_bins"] if window // size <= 10_000), settings["zoom_bins"][-1])
        decoded = decode_zoom_standalone(blob, target, settings["scale"], settings["zoom_bins"])
        return _format_zoom(decoded, plot_type, start, end) if decoded else None

    def _mag_zoom_parts(self, request, zooms, settings, plot_type, contig_level):
        parts = []
        for member in request.members:
            key = (member.contig_id, None if contig_level else request.sample_id)
            blob = zooms.get(key)
            if not blob:
                parts.append(
                    {
                        "bin_start": np.asarray([member.offset, member.offset + member.length - 1]),
                        "bin_end": np.asarray([member.offset, member.offset + member.length - 1]),
                        "mean": np.asarray([0.0, 0.0]),
                        "sparse": False,
                    }
                )
                continue
            window = request.region.length - 1
            target = next(
                (size for size in settings["zoom_bins"] if window // size <= 10_000), settings["zoom_bins"][-1]
            )
            decoded = decode_zoom_standalone(blob, target, settings["scale"], settings["zoom_bins"])
            if decoded:
                decoded["bin_start"] += member.offset
                decoded["bin_end"] += member.offset
                parts.append(decoded)
        return parts

    @staticmethod
    def _combine_zoom(parts, plot_type, start, end):
        if not parts:
            return None
        merged = {
            "bin_start": np.concatenate([part["bin_start"] for part in parts]),
            "bin_end": np.concatenate([part["bin_end"] for part in parts]),
            "mean": np.concatenate([part["mean"] for part in parts]),
            "sparse": any(part.get("sparse", False) for part in parts),
        }
        if all("max" in part for part in parts):
            merged["max"] = np.concatenate([part["max"] for part in parts])
        return _format_zoom(merged, plot_type, start, end)

    def _resolve_partners(self, values, variable):
        partner_ids = values.pop("partner_contig_id", None)
        if partner_ids is None:
            return
        if variable in {"direct_repeat_identity", "inverted_repeat_identity"}:
            labels = [None] * len(partner_ids)
            for index in range(0, len(partner_ids) - 1, 2):
                first, last = partner_ids[index : index + 2]
                if first is not None and last is not None and int(first) >= 0 and int(last) >= 0:
                    labels[index] = labels[index + 1] = f"{int(first) + 1}-{int(last) + 1}"
            values["repeat_positions"] = labels
        else:
            ids = {int(value) for value in partner_ids if value is not None and int(value) >= 0}
            names = self.repository.contig_names(ids)
            values["partner_contig"] = [
                names.get(int(value), "unknown") if value is not None and int(value) >= 0 else None
                for value in partner_ids
            ]
