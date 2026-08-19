"""Summary/peruse generation for the active subject and sample scope."""

from __future__ import annotations

import base64
import time
import traceback
from dataclasses import dataclass
from typing import Any, Callable


@dataclass(frozen=True)
class SummaryBindings:
    connection: Any
    widgets: dict[str, Any]
    sample_scope: Any
    filtered_samples: Callable[[str], list[str]]
    carrier: Any
    enable_timing: bool
    timing: Any


class SummaryController:
    def __init__(self, bindings: SummaryBindings) -> None:
        self.bindings = bindings

    def show(self) -> None:
        conn = self.bindings.connection
        widgets = self.bindings.widgets
        views = self.bindings.sample_scope
        _get_filtered_samples_for_contig = self.bindings.filtered_samples
        summary_carrier = self.bindings.carrier
        enable_timing = self.bindings.enable_timing
        timing = self.bindings.timing

        from ..reports.summary import generate_peruse_html
        from ..repositories.summary import SummaryRepository
        from ..services.summary import SummaryDataService, decode_map

        if enable_timing:
            timing.start_phase("Peruse")
            t_peruse = time.perf_counter()

        is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
        repository = SummaryRepository(conn)
        service = SummaryDataService(repository)

        if is_mag_view:
            mag = widgets["mag_select"].value
            if not mag:
                print("[start_bokeh_server] Peruse: No MAG selected", flush=True)
                return
            sample = widgets["sample_select"].value
            sample_names = [sample] if sample else []
            try:
                report = service.report(None, sample_names, mag_name=mag, is_mag_view=True)
                html_content = generate_peruse_html(report, decode_map(repository))
            except Exception:
                print(f"[summary] Exception: {traceback.format_exc()}", flush=True)
                return
            if enable_timing:
                _step = time.perf_counter() - t_peruse
                print(f"[timing] SHOW SUMMARY (MAG view): {_step:.3f}s{timing.tag(_step)}", flush=True)
                timing.summary("Peruse")
        else:
            contig = widgets["contig_select"].value
            if not contig:
                print("[start_bokeh_server] Peruse: No contig selected", flush=True)
                return

            parent_mag = widgets["contig_to_mag"].get(contig)

            has_samples = widgets["has_samples"]
            if not has_samples:
                sample_names = []
            else:
                is_all = views.active == 1

                if is_all:
                    filtered_samples = _get_filtered_samples_for_contig(contig)

                    if not filtered_samples:
                        print("[start_bokeh_server] Peruse: No samples match filters", flush=True)
                        return

                    sample_names = filtered_samples
                else:
                    sample = widgets["sample_select"].value
                    if not sample:
                        print("[start_bokeh_server] Peruse: No sample selected", flush=True)
                        return
                    sample_names = [sample]

            try:
                report = service.report(contig, sample_names, mag_name=parent_mag, is_mag_view=False)
                html_content = generate_peruse_html(report, decode_map(repository))
            except Exception:
                print(f"[summary] Exception: {traceback.format_exc()}", flush=True)
                return
            if enable_timing:
                _step = time.perf_counter() - t_peruse
                print(
                    f"[timing] SHOW SUMMARY (contig view, {len(sample_names)} samples): {_step:.3f}s{timing.tag(_step)}",
                    flush=True,
                )
                timing.summary("Peruse")

        if not html_content:
            return

        b64 = base64.b64encode(html_content.encode("utf-8")).decode("ascii")
        summary_carrier.text = b64
