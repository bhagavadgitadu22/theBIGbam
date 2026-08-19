"""ONE/ALL sample-scope transitions and their presentation invariants."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable

from bokeh.io import curdoc

from ..models.sample_scope import apply_changed_attributes, visibility_for_sample_scope
from .state import AvailabilitySnapshot, PlotController, PlotUiState, SampleScope, SubjectScope


@dataclass(frozen=True)
class ScopeTransitionBindings:
    widgets: dict[str, Any]
    sample_scope: Any
    sample_section: Any
    variables_one: Any
    variables_all: Any
    sample_parameters_header: Any
    sample_parameters_content: Any
    mag_sort_sample_row: Any
    mag_sort_category: Any
    interaction_lock: dict[str, bool]
    diagnostics: Any
    enable_timing: bool
    timing: Any
    send_timing_ping: Callable[..., None]
    set_operation: Callable[[str], None]
    update_completions: Callable[[Any, list[str]], None]
    refresh_availability: Callable[[], None]


class ScopeTransitionController:
    def __init__(self, bindings: ScopeTransitionBindings) -> None:
        self.bindings = bindings
        widgets = bindings.widgets
        initial = AvailabilitySnapshot(
            samples=tuple(widgets["sample_select"].options),
            contigs=tuple(widgets["contig_select"].options),
            mags=tuple(widgets["mag_select"].options) if widgets["has_mags"] else (),
        )
        self.plot_controller = PlotController(
            PlotUiState(
                sample_scope=SampleScope.ONE,
                subject_scope=(
                    SubjectScope.MAG
                    if widgets["has_mags"] and widgets["view_radio"].active == 0
                    else SubjectScope.CONTIG
                ),
                sample=widgets["sample_select"].value,
                contig=widgets["contig_select"].value,
                mag=widgets["mag_select"].value if widgets["has_mags"] else "",
            ),
            {SampleScope.ONE: initial, SampleScope.ALL: initial},
            project=self._project,
        )

    def _project(self, state: PlotUiState, snapshot: AvailabilitySnapshot) -> None:
        bindings = self.bindings
        bindings.update_completions(bindings.widgets["sample_select"], list(snapshot.samples))
        bindings.update_completions(bindings.widgets["contig_select"], list(snapshot.contigs))
        if bindings.widgets["has_mags"]:
            bindings.update_completions(bindings.widgets["mag_select"], list(snapshot.mags))

    def variable_callback(self, group: Any) -> Callable[..., None]:
        def callback(attr: str, old: list[int], new: list[int]) -> None:
            bindings = self.bindings
            if bindings.interaction_lock["locked"] or bindings.sample_scope.active != 1:
                return
            added = set(new or ()) - set(old or ())
            selected = next(iter(added)) if added else (new[-1] if new else None)
            document = curdoc()
            document.hold("combine")
            bindings.interaction_lock["locked"] = True
            try:
                for other in bindings.widgets["variables_widgets_all"]:
                    target = [selected] if other is group and selected is not None else []
                    if other.active != target:
                        other.active = target
            finally:
                bindings.interaction_lock["locked"] = False
                document.unhold()

        return callback

    def view_changed(self, attr: str, old: int, new: int) -> None:
        bindings = self.bindings
        bindings.set_operation("view_change")
        if bindings.enable_timing:
            bindings.timing.start_phase("View change")
            started = time.perf_counter()
        scope = SampleScope(new)
        visibility = visibility_for_sample_scope(scope, mag_sort_category=bindings.mag_sort_category.value)
        diagnostic_started = time.perf_counter()
        document = curdoc()
        document.hold("combine")
        bindings.interaction_lock["locked"] = True
        try:
            self.plot_controller.set_snapshot(
                SampleScope(old),
                AvailabilitySnapshot(
                    samples=tuple(bindings.widgets["sample_select"].options),
                    contigs=tuple(bindings.widgets["contig_select"].options),
                    mags=tuple(bindings.widgets["mag_select"].options) if bindings.widgets["has_mags"] else (),
                ),
            )
            changed = apply_changed_attributes(
                (
                    (bindings.sample_section, "visible", visibility.sample_section),
                    (bindings.variables_one, "visible", visibility.variables_one),
                    (bindings.variables_all, "visible", visibility.variables_all),
                    (bindings.sample_parameters_header, "visible", visibility.sample_params_header),
                    (bindings.sample_parameters_content, "visible", visibility.sample_params_content),
                    (bindings.mag_sort_sample_row, "visible", visibility.mag_sort_sample_row),
                )
            )
            self.plot_controller.switch_sample_scope(scope)
            # Stored scope snapshots may predate a filter change. Recompute
            # from the current filter after projection; sample availability is
            # refreshed first by the callback so an invalid sample cannot
            # clear an otherwise-valid MAG through reciprocal callbacks.
            bindings.refresh_availability()
            self.plot_controller.set_snapshot(
                scope,
                AvailabilitySnapshot(
                    samples=tuple(bindings.widgets["sample_select"].options),
                    contigs=tuple(bindings.widgets["contig_select"].options),
                    mags=tuple(bindings.widgets["mag_select"].options) if bindings.widgets["has_mags"] else (),
                ),
            )
            bindings.diagnostics.record(
                "view_change", diagnostic_started, query_count=0, callbacks=1, property_updates=changed
            )
            if bindings.enable_timing:
                elapsed = time.perf_counter() - started
                print(f"[timing] on_view_change (server): {elapsed:.3f}s{bindings.timing.tag(elapsed)}", flush=True)
                bindings.timing.summary("View change")
                bindings.send_timing_ping("view_change")
        finally:
            bindings.interaction_lock["locked"] = False
            document.unhold()
