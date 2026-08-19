"""HTML rendering for plain summary-report models."""

from __future__ import annotations

import html
import traceback

from ..repositories.summary import SummaryRepository
from ..services.summary import (
    SummaryDataService,
    column_scales,
    decode_map,
    filter_encode,
    round_significant,
)


def _get_column_scales(conn):
    return column_scales(conn)


def _make_decode_map(conn):
    return decode_map(conn)


def _make_filter_encode(conn):
    return filter_encode(conn)


def _build_row_table_html(rows, col_names, skip_cols=(), decode_map=None):
    if not rows:
        return ""
    decoders = decode_map or {}
    included = [(index, name) for index, name in enumerate(col_names) if name not in set(skip_cols)]
    parts = ['<table style="border-collapse:collapse;width:100%;background:white">', "<thead><tr>"]
    parts.extend(f'<th style="border:1px solid #ddd;padding:8px">{html.escape(str(name))}</th>' for _, name in included)
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        for index, name in included:
            value = row[index]
            if value is not None and name in decoders:
                value = decoders[name](value)
            parts.append(f'<td style="border:1px solid #ddd;padding:8px">{html.escape(str(value or ""))}</td>')
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return "".join(parts)


def generate_summary_table_html(data, column_list):
    columns = ["Sample", *(column for column in column_list if column in data)]
    if len(columns) == 1:
        return None
    rows = tuple(zip(*(data[column] for column in columns)))
    labels = ("Sample", *(column_list[column] for column in columns[1:]))
    formatted = tuple(tuple(round_significant(value, 2) for value in row) for row in rows)
    return _build_row_table_html(formatted, labels)


def _render_metrics(sections):
    content = []
    for section in sections:
        table = generate_summary_table_html(section.data, section.labels)
        if not table:
            continue
        content.extend((f"<b>{html.escape(section.title)}:</b>", table))
        if section.note:
            content.append(f"<i>{html.escape(section.note)}</i>")
    return content


def build_summary_data(conn, contig_name, sample_names):
    service = SummaryDataService(SummaryRepository(conn))
    return _render_metrics(service.metric_sections("Contig", contig_name, sample_names))


def build_mag_summary_data(conn, mag_name, sample_names):
    service = SummaryDataService(SummaryRepository(conn))
    return _render_metrics(service.metric_sections("MAG", mag_name, sample_names))


def generate_peruse_html(conn, contig_name, sample_names, *, mag_name=None, is_mag_view=False):
    try:
        service = SummaryDataService(SummaryRepository(conn))
        report = service.report(contig_name, sample_names, mag_name=mag_name, is_mag_view=is_mag_view)
        decode = decode_map(conn)
        parts = [
            "<!DOCTYPE html><html><head><meta charset='utf-8'>",
            f"<title>theBIGbam - {html.escape(report.subject)} Summary</title>",
            "<style>body{font-family:Arial,sans-serif;margin:20px;background:#f5f5f5}h2,h3{color:#333}b{display:block;margin-top:15px;margin-bottom:5px}i{display:block;margin-bottom:15px;color:#666}</style>",
            "</head><body>",
            f"<h2>{html.escape(report.subject)}</h2>",
            _build_row_table_html(
                report.subject_table.rows, report.subject_table.columns, report.subject_table.skip, decode
            ),
        ]
        if report.parent_mag:
            parts.extend(
                (
                    f"<h3>Parent MAG ({html.escape(report.parent_mag_name or '')})</h3>",
                    _build_row_table_html(
                        report.parent_mag.rows, report.parent_mag.columns, report.parent_mag.skip, decode
                    ),
                )
            )
        if report.contigs:
            parts.extend(
                (
                    "<h3>Contigs</h3>",
                    _build_row_table_html(report.contigs.rows, report.contigs.columns, report.contigs.skip, decode),
                )
            )
        if report.samples.rows:
            parts.extend(
                (
                    "<h3>Samples</h3>",
                    _build_row_table_html(report.samples.rows, report.samples.columns, report.samples.skip),
                )
            )
        parts.extend(_render_metrics(report.metrics))
        parts.append("</body></html>")
        return "\n".join(parts)
    except Exception:
        print(f"[summary] Exception: {traceback.format_exc()}", flush=True)
        return None


def round_to_n_sigfigs(value, n=2):
    return round_significant(value, n)
