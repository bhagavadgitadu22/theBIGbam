"""Stylesheet adapters for the Bokeh and Panel component boundaries."""

from __future__ import annotations

from typing import Any

import panel as pn
from bokeh.layouts import row
from bokeh.models import Div, Spacer, Tooltip

CONTROL_GAP_PX = 4
ACTION_CONTROL_MARGIN = (0, 0, 0, CONTROL_GAP_PX)


def right_panel_tooltip(content: str) -> Tooltip:
    """Place right-sidebar help inward so it remains inside the viewport."""
    return Tooltip(content=content, position="left")


def panel_stylesheet(stylesheet: Any) -> Any:
    """Return inline CSS for Bokeh stylesheet models passed to Panel widgets."""
    css = getattr(stylesheet, "css", None)
    return css if isinstance(css, str) else stylesheet


def bokeh_control_row(*children: Any, **params: Any) -> Any:
    """Build a Bokeh row with the application's standard sibling gap."""
    classes = list(params.pop("css_classes", ()))
    params["css_classes"] = [*classes, "control-row"]
    params.setdefault("spacing", CONTROL_GAP_PX)
    return row(*children, **params)


def panel_control_row(*children: Any, **params: Any) -> pn.Row:
    """Build a Panel row whose sibling gap comes from the shared stylesheet."""
    gap = params.pop("control_gap", CONTROL_GAP_PX)
    classes = list(params.pop("css_classes", ()))
    params["css_classes"] = [*classes, "control-row"]
    styles = dict(params.pop("styles", {}))
    styles.setdefault("column-gap", f"{gap}px")
    styles.setdefault("gap", f"{gap}px")
    params["styles"] = styles
    return pn.Row(*children, **params)


def section_header(title: str | Div, *, toggle: Any | None = None, **params: Any) -> Any:
    """Build a consistently centered sidebar section title.

    Collapsible headers reserve an equal-width slot opposite the toggle so
    the title remains centered across the full sidebar rather than merely in
    the space remaining after the toggle.
    """
    title_model = (
        title
        if isinstance(title, Div)
        else Div(text=f"<span style='font-size: 1.2em;'><b>{title}</b></span>", align="center")
    )
    title_model.css_classes = [*title_model.css_classes, "section-title"]
    title_model.sizing_mode = "stretch_width"
    children: list[Any] = [title_model]
    if toggle is not None:
        toggle_width = toggle.width or 20
        children = [toggle, title_model, Spacer(width=toggle_width, height=toggle.height)]
    classes = list(params.pop("css_classes", ()))
    return bokeh_control_row(
        *children,
        sizing_mode="stretch_width",
        align="center",
        css_classes=[*classes, "section-header"],
        **params,
    )
