"""Stylesheet adapters for the Bokeh and Panel component boundaries."""

from __future__ import annotations

from typing import Any


def panel_stylesheet(stylesheet: Any) -> Any:
    """Return inline CSS for Bokeh stylesheet models passed to Panel widgets."""
    css = getattr(stylesheet, "css", None)
    return css if isinstance(css, str) else stylesheet
