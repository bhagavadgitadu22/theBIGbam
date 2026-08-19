"""Validate export requests before invoking DuckDB export repositories."""

from __future__ import annotations

from ..repositories.exports import download_contig_metrics_csv as _export_contig
from ..repositories.exports import download_mag_metrics_csv as _export_mag


def export_contig_metrics(db_path: str, contig_name: str, sample_names) -> str | None:
    if not db_path or not contig_name:
        return None
    return _export_contig(db_path, contig_name, tuple(dict.fromkeys(sample_names)))


def export_mag_metrics(db_path: str, mag_name: str, sample_names) -> str | None:
    if not db_path or not mag_name:
        return None
    return _export_mag(db_path, mag_name, tuple(dict.fromkeys(sample_names)))
