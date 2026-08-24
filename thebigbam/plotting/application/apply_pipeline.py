"""Typed build, dispatch, and presentation boundaries for Apply actions."""

from __future__ import annotations

import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, ClassVar, Generic, TypeVar

import panel as pn

from ..shared.timing import estimate_grid_data_size, rss_mb
from .apply_inputs import ApplyInputs, collect_apply_inputs


class ApplyMode(Enum):
    CONTIG_ONE = "contig_one"
    CONTIG_ALL = "contig_all"
    MAG_ONE = "mag_one"
    MAG_ALL = "mag_all"


@dataclass(frozen=True)
class ApplyRequest:
    inputs: ApplyInputs
    mode: ClassVar[ApplyMode]


@dataclass(frozen=True)
class ContigOneApplyRequest(ApplyRequest):
    mode: ClassVar[ApplyMode] = ApplyMode.CONTIG_ONE


@dataclass(frozen=True)
class ContigAllApplyRequest(ApplyRequest):
    mode: ClassVar[ApplyMode] = ApplyMode.CONTIG_ALL


@dataclass(frozen=True)
class MagOneApplyRequest(ApplyRequest):
    mode: ClassVar[ApplyMode] = ApplyMode.MAG_ONE


@dataclass(frozen=True)
class MagAllApplyRequest(ApplyRequest):
    mode: ClassVar[ApplyMode] = ApplyMode.MAG_ALL


REQUEST_TYPES = {
    ApplyMode.CONTIG_ONE: ContigOneApplyRequest,
    ApplyMode.CONTIG_ALL: ContigAllApplyRequest,
    ApplyMode.MAG_ONE: MagOneApplyRequest,
    ApplyMode.MAG_ALL: MagAllApplyRequest,
}


@dataclass(frozen=True)
class ValidationFailure:
    message: str


@dataclass(frozen=True)
class RequestBuildResult:
    request: ApplyRequest | None = None
    failure: ValidationFailure | None = None

    @property
    def is_valid(self) -> bool:
        return self.request is not None


@dataclass(frozen=True)
class PlotResult:
    grid: Any
    toolbar_row: Any
    operation: str
    description: str
    total_label: str
    profiler: Any = None


@dataclass(frozen=True)
class RangeSnapshot:
    preserve: bool
    start: float | None = None
    end: float | None = None


class ApplyRequestBuilder:
    """Snapshot and validate widgets before database or rendering work starts."""

    def __init__(self, bindings: Any) -> None:
        self.bindings = bindings

    def from_widgets(self) -> RequestBuildResult:
        inputs = collect_apply_inputs(self.bindings)
        widgets = self.bindings.widgets
        is_mag = bool(widgets["has_mags"] and widgets["view_radio"].active == 0)
        mode = {
            (False, False): ApplyMode.CONTIG_ONE,
            (False, True): ApplyMode.CONTIG_ALL,
            (True, False): ApplyMode.MAG_ONE,
            (True, True): ApplyMode.MAG_ALL,
        }[(is_mag, inputs.is_all)]
        subject = widgets["mag_select"].value if is_mag else inputs.contig
        if not subject:
            label = "MAG" if is_mag else "contig"
            return RequestBuildResult(failure=ValidationFailure(f"Please select a {label}."))
        position_failure = self._validate_positions(is_mag, inputs.contig)
        if position_failure:
            return RequestBuildResult(failure=position_failure)
        selection_failure = self._validate_selection(mode, inputs)
        if selection_failure:
            return RequestBuildResult(failure=selection_failure)
        return RequestBuildResult(request=REQUEST_TYPES[mode](inputs))

    def _validate_positions(self, is_mag: bool, contig: str) -> ValidationFailure | None:
        start_text = self.bindings.from_position_input.value.strip()
        end_text = self.bindings.to_position_input.value.strip()
        try:
            start = int(start_text) if start_text else 1
            if end_text:
                end = int(end_text)
            elif is_mag:
                end = None
            else:
                end = self.bindings.widgets["contig_lengths"].get(contig, 0)
        except ValueError:
            return ValidationFailure("Invalid position range - positions must be integers.")
        if end is not None and start >= end:
            return ValidationFailure("Invalid position range - start must be less than end.")
        return None

    def _validate_selection(self, mode: ApplyMode, inputs: ApplyInputs) -> ValidationFailure | None:
        if mode not in (ApplyMode.CONTIG_ALL, ApplyMode.MAG_ALL):
            return None
        genome_features = self._selected_labels(self.bindings.combined_features_cbg)
        selected_variable = any(self._selected_labels(group) for group in inputs.active_variables_widgets)
        if mode is ApplyMode.MAG_ALL:
            genome_features = [feature for feature in genome_features if feature != "Gene map"]
        if not selected_variable and not genome_features:
            return ValidationFailure("In 'All samples' view you must select at least one variable.")
        if mode is ApplyMode.MAG_ALL:
            category = self.bindings.mag_params_category_select.value
            source = self.bindings._mag_sort_category_sources.get(category)
            if category != "Contig" and source and not self.bindings.mag_params_sort_sample_select.value:
                return ValidationFailure(
                    "Please select a sample to use for 'Order contigs by'; "
                    "the selected sort category requires per-sample values."
                )
        return None

    @staticmethod
    def _selected_labels(group: Any) -> list[str]:
        if group is None:
            return []
        return [group.labels[index] for index in group.active]


T = TypeVar("T")


class ApplyDispatcher(Generic[T]):
    """Route each Apply mode explicitly, without inspecting widgets."""

    def __init__(self, handlers: dict[ApplyMode, Callable[[ApplyRequest, float], T]]) -> None:
        missing = set(ApplyMode) - handlers.keys()
        if missing:
            raise ValueError(f"Missing Apply handlers: {sorted(mode.value for mode in missing)}")
        self._handlers = handlers

    def render(self, request: ApplyRequest, started_at: float) -> T:
        return self._handlers[request.mode](request, started_at)


class PlotPresenter:
    """Own validation/error display and replacement of the visible plot container."""

    def __init__(self, bindings: Any) -> None:
        self.bindings = bindings

    def show_validation(self, failure: ValidationFailure) -> None:
        self._hide_action_buttons()
        self.bindings.main_placeholder.objects = [pn.pane.HTML(f"<pre>Error: {failure.message}</pre>")]

    def show_exception(self, traceback_text: str) -> None:
        self._hide_action_buttons()
        self.bindings.main_placeholder.objects = [pn.pane.HTML(f"<pre>Error building plot:\n{traceback_text}</pre>")]

    def bind_range_inputs(self, x_range: Any) -> None:
        """Replace the callbacks which project the visible range into position controls."""
        if x_range is None:
            self.bindings.current_plot_state["range_callbacks"] = ()
            return

        def sync_start(attr: str, old: Any, new: Any) -> None:
            self.bindings.from_position_input.value = str(int(new))

        def sync_end(attr: str, old: Any, new: Any) -> None:
            self.bindings.to_position_input.value = str(int(new))

        x_range.on_change("start", sync_start)
        x_range.on_change("end", sync_end)
        self.bindings.current_plot_state["range_callbacks"] = (
            (x_range, "start", sync_start),
            (x_range, "end", sync_end),
        )

    def capture_range(self, subject: str, sample: str | None, is_all: bool, start: int, end: int | None):
        state = self.bindings.current_plot_state
        shared = state["shared_xrange"]
        preserve = (
            shared is not None
            and state["contig"] == subject
            and state["is_all"] == is_all
            and (is_all or state["sample"] == sample)
            and state["data_xstart"] == start
            and state["data_xend"] == end
        )
        snapshot = RangeSnapshot(preserve, shared.start if preserve else None, shared.end if preserve else None)
        # Capture the old range before detaching its callbacks and releasing it.
        # Render handlers call this exactly once, before installing the new range.
        self.bindings.plot_lifecycle.prepare_replacement()
        return snapshot

    def install_range(
        self,
        grid: Any,
        snapshot: RangeSnapshot,
        *,
        subject: str,
        sample: str | None,
        is_all: bool,
        start: int,
        end: int | None,
    ):
        x_range = self.bindings.plot_lifecycle.shared_xrange(grid)
        if snapshot.preserve and x_range is not None:
            x_range.start = snapshot.start
            x_range.end = snapshot.end
        state = self.bindings.current_plot_state
        state["contig"] = subject
        state["sample"] = sample
        state["is_all"] = is_all
        state["shared_xrange"] = x_range
        state["data_xstart"] = start
        state["data_xend"] = end
        self.bind_range_inputs(x_range)
        return x_range

    def mag_toolbar(self, dots_html: str, has_samples: bool):
        bindings = self.bindings
        bindings.peruse_button.visible = True
        bindings.download_mag_metrics_button.visible = bool(has_samples)
        bindings.download_data_button.visible = True
        bindings.command_hint_pane.visible = False
        return pn.Row(
            pn.pane.HTML(dots_html, align="center", margin=(0, 5, 0, 0)),
            pn.Spacer(sizing_mode="stretch_width"),
            bindings.peruse_button,
            bindings.download_mag_metrics_button,
            bindings.download_data_button,
            margin=(0, 0, 5, 0),
        )

    def contig_toolbar(self, has_samples: bool, has_mags: bool):
        bindings = self.bindings
        bindings.peruse_button.visible = True
        bindings.download_data_button.visible = True
        bindings.download_metrics_button.visible = bool(has_samples)
        bindings.download_mag_metrics_button.visible = bool(has_samples and has_mags)
        bindings.command_hint_pane.visible = False
        return pn.Row(
            pn.Spacer(sizing_mode="stretch_width"),
            bindings.peruse_button,
            bindings.download_mag_metrics_button,
            bindings.download_metrics_button,
            bindings.download_data_button,
            margin=(0, 0, 5, 0),
        )

    def replace(self, result: PlotResult, started_at: float) -> None:
        bindings = self.bindings
        presentation_started = time.perf_counter()
        bindings.set_operation(result.operation)
        payload = None
        if bindings.enable_timing or bindings.diagnostics.enabled:
            payload = estimate_grid_data_size(result.grid)
        if bindings.enable_timing:
            print(f"[timing] RSS after plot generation: {rss_mb():.0f} MB", flush=True)
            print(f"[timing] Sending: {result.description}", flush=True)
            data_bytes, n_sources, n_models = payload
            print(
                f"[timing] Data to frontend: {data_bytes / 1024 / 1024:.1f} MB approx "
                f"({n_models} models, {n_sources} data sources)",
                flush=True,
            )
            sent_at = time.perf_counter()
        bindings.main_placeholder.objects = [
            pn.Column(result.toolbar_row, bindings.command_hint_pane, result.grid, sizing_mode="stretch_both")
        ]
        if bindings.diagnostics.enabled:
            data_bytes, n_sources, n_models = payload
            bindings.diagnostics.record(
                "apply",
                started_at,
                models=n_models,
                data_sources=n_sources,
                payload_bytes=data_bytes,
                callbacks=len(bindings.current_plot_state["range_callbacks"]),
            )
        if bindings.enable_timing:
            elapsed = time.perf_counter() - sent_at
            print(f"[timing] Panel objects update: {elapsed:.3f}s{bindings.timing.tag(elapsed)}", flush=True)
            elapsed = time.perf_counter() - started_at
            print(f"[timing] Total APPLY ({result.total_label}): {elapsed:.3f}s{bindings.timing.tag(0)}", flush=True)
        if result.profiler is not None:
            result.profiler.seconds["presentation"] = time.perf_counter() - presentation_started
            result.profiler.report(bindings.data_cache.stats())

    def _hide_action_buttons(self) -> None:
        self.bindings.peruse_button.visible = False
        self.bindings.download_metrics_button.visible = False
        self.bindings.download_mag_metrics_button.visible = False
        self.bindings.download_data_button.visible = False
