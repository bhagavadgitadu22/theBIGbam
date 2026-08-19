"""Lifecycle management for replacing interactive Bokeh plot trees."""

from __future__ import annotations

from typing import Any

from ..models.session import CurrentPlotState


class PlotLifecycle:
    """Own callback cleanup and shared-range discovery for the active plot."""

    def __init__(self, state: CurrentPlotState | None = None) -> None:
        self.state = state or CurrentPlotState()

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
        """Return the first figure range in a Bokeh grid, when present."""
        for child_spec in getattr(grid, "children", ()):
            child = child_spec[0] if isinstance(child_spec, tuple) else child_spec
            if hasattr(child, "x_range"):
                return child.x_range
        return None
