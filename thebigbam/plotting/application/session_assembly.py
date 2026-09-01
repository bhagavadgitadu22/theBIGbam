"""Final layout, Apply controller, and inspectable session assembly."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..models.session import PlotSessionContext
from .apply_controller import ApplyController
from .layout import LayoutParts, assemble_layout
from .wiring import make_apply_bindings


@dataclass(frozen=True)
class FinalizedSession:
    layout: Any
    placeholder: Any
    controls: Any
    history: Any
    apply_controller: ApplyController


def finalize_session(
    *,
    parts: LayoutParts,
    has_samples: bool,
    summary_carrier: Any,
    stylesheet: Any,
    timing_models: tuple[Any, ...],
    apply_arguments: dict[str, Any],
    history_drawer: Any,
) -> FinalizedSession:
    assembled = assemble_layout(
        parts,
        has_samples=has_samples,
        summary_carrier=summary_carrier,
        stylesheet=stylesheet,
        timing_models=timing_models,
        history_drawer=history_drawer,
    )
    apply_arguments = {**apply_arguments, "placeholder": assembled.placeholder}
    controller = ApplyController(make_apply_bindings(**apply_arguments))
    return FinalizedSession(assembled.layout, assembled.placeholder, assembled.controls, assembled.history, controller)


def inspectable_application(
    application_type,
    *,
    session: FinalizedSession,
    connection: Any,
    preloaded: Any,
    widgets: dict[str, Any],
    diagnostics: Any,
    plot_state: Any,
    plot_controller: Any,
    apply: Any,
):
    exposed = {**widgets, "sample_scope": widgets["sample_scope"]}
    context = PlotSessionContext(connection, preloaded, exposed, diagnostics, plot_state)
    return application_type(session.layout, context, plot_controller, exposed, diagnostics, apply)
