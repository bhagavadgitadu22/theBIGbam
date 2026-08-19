"""Construction of summary and download controls for a plotting session."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import panel as pn
from bokeh.models import CustomJS, Div

from ..downloads.callbacks import make_contig_metrics_callback, make_mag_metrics_callback
from ..downloads.inspect_command import InspectCommandBindings, InspectCommandController
from .summary import SummaryBindings, SummaryController


@dataclass(frozen=True)
class OutputControls:
    peruse_button: Any
    summary_controller: SummaryController
    summary_carrier: Any
    contig_metrics_button: Any
    mag_metrics_button: Any
    data_button: Any
    command_hint: Any
    download_widgets: dict[str, Any]


def build_output_controls(
    *,
    db_path: str,
    connection: Any,
    widgets: dict[str, Any],
    sample_scope: Any,
    filtered_samples: Callable[[str], list[str]],
    combined_features: Any,
    subplot_variables: Any,
    from_position: Any,
    to_position: Any,
    stylesheet: Any,
    enable_timing: bool,
    timing: Any,
    report_timing: Callable[[str, float], None] | None,
    show_summary: Callable[..., None],
) -> OutputControls:
    peruse = pn.widgets.Button(
        name="SHOW SUMMARY",
        height=30,
        stylesheets=[stylesheet],
        css_classes=["apply-btn"],
        visible=False,
    )
    peruse.on_click(show_summary)
    carrier = Div(text="", visible=False)
    carrier.js_on_change(
        "text",
        CustomJS(
            code="""
        var b64 = cb_obj.text;
        if (!b64) return;
        var bin = atob(b64);
        var bytes = new Uint8Array(bin.length);
        for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
        var blob = new Blob([bytes], {type: 'text/html;charset=utf-8'});
        window.open(URL.createObjectURL(blob), '_blank');
        cb_obj.text = '';
    """
        ),
    )
    summary = SummaryController(
        SummaryBindings(connection, widgets, sample_scope, filtered_samples, carrier, enable_timing, timing)
    )
    downloads: dict[str, Any] = {"contig_metrics": None, "mag_metrics": None, "data": None}
    contig_button = pn.widgets.FileDownload(
        callback=make_contig_metrics_callback(
            db_path, widgets, sample_scope, filtered_samples, downloads, report_timing
        ),
        filename="contig_metrics.csv",
        label="DOWNLOAD CONTIG METRICS",
        button_type="primary",
        height=30,
        visible=False,
    )
    downloads["contig_metrics"] = contig_button
    mag_button = pn.widgets.FileDownload(
        callback=make_mag_metrics_callback(db_path, widgets, downloads, report_timing),
        filename="mag_metrics.csv",
        label="DOWNLOAD MAG METRICS",
        button_type="primary",
        height=30,
        visible=False,
    )
    downloads["mag_metrics"] = mag_button
    inspect = InspectCommandController(
        InspectCommandBindings(
            db_path,
            widgets,
            sample_scope,
            combined_features,
            subplot_variables,
            from_position,
            to_position,
            filtered_samples,
        )
    )
    downloads["data"] = inspect.button
    return OutputControls(
        peruse,
        summary,
        carrier,
        contig_button,
        mag_button,
        inspect.button,
        inspect.pane,
        downloads,
    )
