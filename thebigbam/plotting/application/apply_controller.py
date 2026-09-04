"""Thin build, dispatch, and presentation controller for Apply actions."""

from __future__ import annotations

import time
import traceback

from bokeh.io import curdoc

from ..shared.timing import rss_mb
from .apply_handlers import (
    ApplyBindings,
    ApplyRenderEngine,
    ContigAllHandler,
    ContigOneHandler,
    MagAllHandler,
    MagOneHandler,
)
from .apply_pipeline import ApplyDispatcher, ApplyMode, ApplyRequestBuilder, PlotPresenter

__all__ = ["ApplyBindings", "ApplyController"]


class ApplyController:
    """Coordinate one Apply transaction without performing plot composition."""

    def __init__(self, bindings: ApplyBindings) -> None:
        self.bindings = bindings
        self.request_builder = ApplyRequestBuilder(bindings)
        self.plot_presenter = PlotPresenter(bindings)
        engine = ApplyRenderEngine(bindings, self.plot_presenter)
        self.handlers = {
            ApplyMode.CONTIG_ONE: ContigOneHandler(engine),
            ApplyMode.CONTIG_ALL: ContigAllHandler(engine),
            ApplyMode.MAG_ONE: MagOneHandler(engine),
            ApplyMode.MAG_ALL: MagAllHandler(engine),
        }
        self.dispatcher = ApplyDispatcher({mode: handler.render for mode, handler in self.handlers.items()})

    def apply(self) -> bool:
        bindings = self.bindings
        doc = curdoc()
        doc.hold("combine")
        print(flush=True)
        started_at = time.perf_counter()
        outcome = "error"
        try:
            bindings.set_operation("apply/param_parse")
            bindings.diagnostics.next_generation()
            if bindings.enable_timing:
                bindings.timing.start_phase("APPLY")
                print(f"[timing] Memory (current RSS) at APPLY start: {rss_mb():.0f} MB", flush=True)
                started_at = time.perf_counter()
            build_result = self.request_builder.from_widgets()
            if build_result.failure is not None:
                self.plot_presenter.show_validation(build_result.failure)
                return False
            result = self.dispatcher.render(build_result.request, started_at)
            self.plot_presenter.replace(result, started_at)
            outcome = "result"
            return True
        except Exception:
            traceback_text = traceback.format_exc()
            print(f"[start_bokeh_server] Exception: {traceback_text}", flush=True)
            self.plot_presenter.show_exception(traceback_text)
            return False
        finally:
            bindings.set_operation("idle")
            if bindings.enable_timing:
                print(f"[timing] Memory (current RSS) at APPLY end: {rss_mb():.0f} MB", flush=True)
                bindings.timing.summary("APPLY")
            if callable(bindings._send_timing_ping):
                bindings._send_timing_ping(
                    f"APPLY_{outcome}", started_at if bindings.enable_timing else None
                )
            doc.unhold()
