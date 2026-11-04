"""Lifecycle management for replacing interactive Bokeh plot trees."""

from __future__ import annotations

from typing import Any

from ..models.session import CurrentPlotState


class PlotLifecycle:
    """Own callback cleanup and shared-range discovery for the active plot."""

    def __init__(self, state: CurrentPlotState | None = None) -> None:
        self.state = state if state is not None else CurrentPlotState()

    def prepare_replacement(self) -> None:
        """Detach callbacks and release the range retained by the previous plot."""
        for model, attribute, callback in self.state.get("range_callbacks", ()):
            try:
                model.remove_on_change(attribute, callback)
            except (ValueError, RuntimeError):
                pass
        self.state["range_callbacks"] = ()
        self.state["shared_xrange"] = None

    @staticmethod
    def shared_xrange(grid: Any) -> Any | None:
        """Find the first figure range, including figures inside nested rows."""
        pending = [grid]
        visited = set()
        while pending:
            model = pending.pop(0)
            if id(model) in visited:
                continue
            visited.add(id(model))
            if hasattr(model, "x_range"):
                return model.x_range
            for child_spec in getattr(model, "children", ()):
                child = child_spec[0] if isinstance(child_spec, tuple) else child_spec
                pending.append(child)
        return None
