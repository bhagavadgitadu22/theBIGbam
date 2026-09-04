"""Assemble the stable outer Panel layout for the plotting application."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import panel as pn
from bokeh.models import Div

from ..controls.panel_resizer import PanelResizer
from ..shared.styles import panel_stylesheet


def separator() -> Div:
    return Div(
        text="",
        height=2,
        min_height=2,
        max_height=2,
        sizing_mode="stretch_width",
        margin=(10, 0),
        styles={"background-color": "#333", "box-sizing": "border-box"},
        css_classes=["section-separator"],
    )


@dataclass(frozen=True)
class LayoutParts:
    logo: Any
    sample_scope: Any
    filtering_header: Any
    filtering_content: Any
    mag_separator: Any
    mag_header: Any
    view_radio: Any
    mag_select: Any
    contig_header: Any
    contig_select: Any
    below_contig_content: Any
    sample_section: Any
    variables_one: Any
    variables_all: Any
    plotting_separator: Any
    plotting_header: Any
    plotting_content: Any
    buttons: Any


@dataclass(frozen=True)
class AssembledLayout:
    layout: Any
    placeholder: Any
    controls: Any
    history: Any
    left_toggle: Any
    right_toggle: Any
    left_resizer: Any
    right_resizer: Any
    left_rail: Any
    right_rail: Any
    left_shell: Any
    right_shell: Any


def assemble_layout(
    parts: LayoutParts,
    *,
    has_samples: bool,
    summary_carrier: Any,
    stylesheet: Any,
    history_drawer: Any,
    timing_models: tuple[Any, ...] = (),
) -> AssembledLayout:
    filtering_separator = separator()
    contig_separator = separator()
    variable_separator = separator()
    if has_samples:
        children = [
            parts.logo,
            parts.sample_scope,
            filtering_separator,
            parts.filtering_header,
            parts.filtering_content,
            parts.mag_separator,
            parts.mag_header,
            parts.view_radio,
            parts.mag_select,
            contig_separator,
            parts.contig_header,
            parts.contig_select,
            parts.below_contig_content,
            parts.sample_section,
            variable_separator,
            parts.variables_one,
            parts.variables_all,
            parts.plotting_separator,
            parts.plotting_header,
            parts.plotting_content,
            parts.buttons,
        ]
        placeholder_text = (
            '<i>No plot yet. Select one sample, one contig and at least one variable in "One sample" mode '
            'or one contig and one variable in "All samples" mode and click Apply.</i>'
        )
    else:
        children = [
            parts.logo,
            parts.filtering_header,
            parts.filtering_content,
            parts.mag_separator,
            parts.mag_header,
            parts.view_radio,
            parts.mag_select,
            contig_separator,
            parts.contig_header,
            parts.contig_select,
            parts.below_contig_content,
            parts.plotting_separator,
            parts.plotting_header,
            parts.plotting_content,
            parts.buttons,
        ]
        placeholder_text = "<i>No plot yet. Select one contig and click Apply to view the genome annotation.</i>"
    children.extend(timing_models)
    children.append(pn.Spacer(height=20, sizing_mode="fixed", css_classes=["sidebar-bottom-space"]))
    sidebar_content = pn.Column(*children, sizing_mode="stretch_width", css_classes=["sidebar-content"])
    controls = pn.Column(sidebar_content, sizing_mode="stretch_height", width=400, css_classes=["left-col"])
    sidebar_content.stylesheets = [stylesheet]
    controls.stylesheets = [stylesheet]
    placeholder = pn.Column(
        pn.pane.HTML(placeholder_text),
        sizing_mode="stretch_both",
        css_classes=["main-right", "plot-area"],
    )
    left_toggle = pn.widgets.Button(
        name="◀",
        width=28,
        height=38,
        css_classes=["drawer-toggle", "drawer-toggle-left"],
        margin=0,
        stylesheets=[panel_stylesheet(stylesheet)],
    )
    history = history_drawer
    history.visible = False
    right_toggle = pn.widgets.Button(
        name="◀",
        width=28,
        height=38,
        css_classes=["drawer-toggle", "drawer-toggle-right"],
        margin=0,
        stylesheets=[panel_stylesheet(stylesheet)],
    )

    def toggle_left() -> None:
        controls.visible = not controls.visible
        left_resizer.enabled = controls.visible
        left_shell.width = left_resizer.value if controls.visible else 0
        left_toggle.name = "◀" if controls.visible else "▶"

    def toggle_right() -> None:
        history.visible = not history.visible
        right_resizer.enabled = history.visible
        right_shell.width = right_resizer.value if history.visible else 0
        right_toggle.name = "▶" if history.visible else "◀"

    left_toggle.on_click(lambda _event: toggle_left())
    right_toggle.on_click(lambda _event: toggle_right())
    left_resizer = PanelResizer(
        side="left",
        storage_key="thebigbam-left-panel-width",
        minimum=100,
        maximum=700,
        default_width=400,
        value=400,
        width=8,
        sizing_mode="stretch_height",
        align="center",
        margin=0,
    )
    right_resizer = PanelResizer(
        side="right",
        storage_key="thebigbam-right-panel-width",
        minimum=100,
        maximum=700,
        default_width=250,
        value=250,
        enabled=False,
        width=8,
        sizing_mode="stretch_height",
        align="center",
        margin=0,
    )

    def resize_left(event) -> None:
        controls.width = event.new
        if controls.visible:
            left_shell.width = event.new

    def resize_right(event) -> None:
        history.width = event.new
        if history.visible:
            right_shell.width = event.new

    left_resizer.param.watch(resize_left, "value")
    right_resizer.param.watch(resize_right, "value")
    left_rail = pn.Column(
        left_resizer,
        left_toggle,
        width=8,
        sizing_mode="stretch_height",
        align="center",
        margin=0,
        css_classes=["panel-resize-rail", "panel-resize-rail-left"],
        styles={"bottom": "0", "position": "absolute", "right": "-8px", "top": "0"},
    )
    right_rail = pn.Column(
        right_resizer,
        right_toggle,
        width=8,
        sizing_mode="stretch_height",
        align="center",
        margin=0,
        css_classes=["panel-resize-rail", "panel-resize-rail-right"],
        styles={"bottom": "0", "left": "-8px", "position": "absolute", "top": "0"},
    )
    left_shell = pn.Column(
        controls,
        left_rail,
        width=400,
        sizing_mode="stretch_height",
        margin=0,
        css_classes=["panel-shell", "panel-shell-left"],
        stylesheets=[stylesheet],
    )
    right_shell = pn.Column(
        history,
        right_rail,
        width=0,
        sizing_mode="stretch_height",
        margin=0,
        css_classes=["panel-shell", "panel-shell-right"],
        stylesheets=[stylesheet],
    )
    layout = pn.Row(
        left_shell,
        placeholder,
        pn.pane.Bokeh(summary_carrier),
        right_shell,
        sizing_mode="stretch_both",
        css_classes=["main-layout"],
    )
    layout.stylesheets = [stylesheet]
    return AssembledLayout(
        layout,
        placeholder,
        controls,
        history,
        left_toggle,
        right_toggle,
        left_resizer,
        right_resizer,
        left_rail,
        right_rail,
        left_shell,
        right_shell,
    )
