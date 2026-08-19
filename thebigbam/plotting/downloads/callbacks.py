"""Panel download callback factories for plotting metrics."""

from __future__ import annotations

import io
import time
from collections.abc import Callable, Mapping
from typing import Any

from .data import download_contig_metrics_csv, download_mag_metrics_csv


def _safe_name(value: str) -> str:
    return "".join(character if character.isalnum() or character in "-_" else "_" for character in value)


def make_contig_metrics_callback(
    db_path: str,
    widgets: Mapping[str, Any],
    sample_scope,
    filtered_samples: Callable[[str], list[str]],
    download_widgets: Mapping[str, Any],
    report_timing: Callable[[str, float], None] | None = None,
):
    def download():
        started = time.perf_counter()
        contig = widgets["contig_select"].value
        if not contig:
            print("[plotting] Download contig metrics: no contig selected", flush=True)
            return io.StringIO("")
        if sample_scope.active == 1:
            samples = filtered_samples(contig)
            if not samples:
                print("[plotting] Download contig metrics: no samples match filters", flush=True)
                return io.StringIO("")
            filename = f"{_safe_name(contig)}_in_all_samples_contig_metrics.csv"
        else:
            sample = widgets["sample_select"].value
            if not sample:
                print("[plotting] Download contig metrics: no sample selected", flush=True)
                return io.StringIO("")
            samples = [sample]
            filename = f"{_safe_name(contig)}_in_{_safe_name(sample)}_contig_metrics.csv"
        widget = download_widgets.get("contig_metrics")
        if widget is not None:
            widget.filename = filename
        content = download_contig_metrics_csv(db_path, contig, samples)
        if report_timing:
            report_timing(f"DOWNLOAD CONTIG METRICS ({len(samples)} samples)", time.perf_counter() - started)
        return io.StringIO(content or "")

    return download


def make_mag_metrics_callback(
    db_path: str,
    widgets: Mapping[str, Any],
    download_widgets: Mapping[str, Any],
    report_timing: Callable[[str, float], None] | None = None,
):
    def download():
        started = time.perf_counter()
        is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
        if is_mag_view:
            mag = widgets["mag_select"].value
        else:
            contig = widgets["contig_select"].value
            mag = widgets["contig_to_mag"].get(contig) if contig else None
        if not mag:
            print("[plotting] Download MAG metrics: no MAG selected or found", flush=True)
            return io.StringIO("")
        sample = widgets["sample_select"].value
        samples = [sample] if sample else []
        content = download_mag_metrics_csv(db_path, mag, samples)
        if report_timing:
            report_timing("DOWNLOAD MAG METRICS", time.perf_counter() - started)
        if content:
            suffix = f"_in_{_safe_name(sample)}" if sample else ""
            widget = download_widgets.get("mag_metrics")
            if widget is not None:
                widget.filename = f"{_safe_name(mag)}{suffix}_mag_metrics.csv"
        return io.StringIO(content or "")

    return download
