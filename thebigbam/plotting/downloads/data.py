"""Download-facing wrappers around CSV export repositories."""

from ..services.exports import export_contig_metrics, export_mag_metrics

download_contig_metrics_csv = export_contig_metrics
download_mag_metrics_csv = export_mag_metrics

__all__ = ["download_contig_metrics_csv", "download_mag_metrics_csv"]
