"""Show BAM matching warnings beside sample selection."""

from __future__ import annotations

from html import escape
from typing import Any, Mapping

from bokeh.models import Div

from ..services.bam_catalog import BamCatalog
from .model_events import on_model_change


def _warning_div(message: str) -> Div:
    return Div(
        text=f"<b>BAM warning:</b> {escape(message)}" if message else "",
        visible=bool(message),
        sizing_mode="stretch_width",
        styles={"color": "#9b4d00", "white-space": "normal"},
    )


def build_bam_warning_controls(
    catalog: BamCatalog,
    sample_select: Any,
    contig_select: Any,
    contig_lengths: Mapping[str, int],
) -> tuple[Div, Div]:
    """Keep startup-wide and selected-sample issues visible in the sidebar."""
    global_issues = tuple(
        warning
        for warning in catalog.warnings
        if "has no matching database sample" in warning or "No .bam files found" in warning
    )
    summary = " ".join(global_issues[:5])
    if len(global_issues) > 5:
        summary += f" And {len(global_issues) - 5} more; see the server log."
    global_status = _warning_div(summary)
    sample_status = _warning_div("")

    def update(_attr: str, _old: Any, _new: Any) -> None:
        name = sample_select.value
        message = (
            catalog.warning_for_sample(name, contig_select.value, contig_lengths.get(contig_select.value, 0))
            if name
            else None
        )
        sample_status.text = f"<b>BAM warning:</b> {escape(message)}" if message else ""
        sample_status.visible = bool(message)

    on_model_change(sample_select, "value", update)
    on_model_change(contig_select, "value", update)
    update("value", None, None)
    return global_status, sample_status
