"""Observe both Bokeh widgets and Panel custom components."""

from __future__ import annotations

from typing import Any, Callable


def on_model_change(model: Any, attribute: str, callback: Callable[[str, Any, Any], None]) -> None:
    if hasattr(model, "on_change"):
        model.on_change(attribute, callback)
    else:
        model.param.watch(lambda event: callback(attribute, event.old, event.new), attribute)
