"""ONE/ALL variable panel construction and module synchronization."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

from bokeh.layouts import column, row
from bokeh.models import Div, Tooltip
from bokeh.models.widgets import Button, HelpButton


@dataclass(frozen=True)
class VariableControls:
    one_sample: Any
    all_samples: Any
    genome_checkbox: Any | None
    genome_index: int | None


def build_variable_controls(
    widgets: Mapping[str, Any],
    toggle_stylesheet: Any,
    make_toggle_callback: Callable[[Any, Any], Callable[..., None]],
) -> VariableControls:
    ## Build Variables section - TWO SEPARATE SECTIONS for each view
    variables_title_one = Div(text="<span style='font-size: 1.2em;'><b>Variables</b></span>")
    variables_title_all = Div(text="<span style='font-size: 1.2em;'><b>Variables</b></span>")
    genome_cbg_one = None  # Will store reference to Genome module's CheckboxButtonGroup (shared between views)

    # Build "One Sample" view variables section
    # Has module checkboxes + collapsible variable groups
    # NOTE: Genome module is built separately and placed in Contigs section
    controls_variables_one = []
    module_toggles_one = []
    module_contents_one = []
    genome_index_one = None  # Track Genome module index for separate handling

    for i, module_widget in enumerate(widgets["module_widgets_one"]):
        module_name = widgets["module_names"][i]
        help_tooltip = widgets["helps_widgets"][i]

        # Skip Genome module - will be added to Contigs section
        if module_name == "Genome":
            genome_index_one = i
            genome_cbg_one = widgets["variables_widgets_one"][i]
            continue

        # Create toggle button for collapsible section
        toggle_btn = Button(
            label="▶",
            width=20,
            height=20,
            button_type="primary",
            align="center",
            margin=0,
            stylesheets=[toggle_stylesheet],
        )
        toggle_btn.styles = {"padding": "0px", "line-height": "20px"}
        module_toggles_one.append(toggle_btn)

        # Build header with checkbox (module_widget has module name as label)
        if help_tooltip is not None:
            help_btn = HelpButton(
                tooltip=help_tooltip,
                width=20,
                height=20,
                align="center",
                button_type="light",
                stylesheets=[toggle_stylesheet],
            )
            hdr = row(
                toggle_btn,
                module_widget,
                help_btn,
                sizing_mode="stretch_width",
                align="center",
                styles={"overflow": "hidden"},
            )
        else:
            hdr = row(toggle_btn, module_widget, sizing_mode="stretch_width", align="center")

        controls_variables_one.append(hdr)

        # Add the module's CheckboxButtonGroup for variables (collapsible, starts folded)
        cbg = widgets["variables_widgets_one"][i]
        cbg.visible = False
        module_contents_one.append(cbg)
        controls_variables_one.append(cbg)

    # Add toggle callbacks for One Sample view
    for i, toggle_btn in enumerate(module_toggles_one):
        content = module_contents_one[i]
        toggle_btn.on_click(make_toggle_callback(toggle_btn, content))

    # Build "All Samples" view variables section
    # NOTE: Genome module is built separately and placed in Contigs section
    # Other modules have title-only headers (no checkbox)
    controls_variables_all = []
    module_toggles_all = []
    module_contents_all = []

    for i in range(len(widgets["module_names"])):
        module_name = widgets["module_names"][i]
        help_tooltip = widgets["helps_widgets"][i]

        # Skip Genome module - it's in the Contigs section (shared between views)
        if module_name == "Genome":
            continue

        # Create toggle button for collapsible section
        toggle_btn = Button(
            label="▶",
            width=20,
            height=20,
            button_type="primary",
            align="center",
            margin=0,
            stylesheets=[toggle_stylesheet],
        )
        toggle_btn.styles = {"padding": "0px", "line-height": "20px"}
        module_toggles_all.append(toggle_btn)

        # Build header with title only (no checkbox) for non-Genome modules
        module_title_div = Div(text=f"{module_name}", align="center")
        if help_tooltip is not None:
            # Create new tooltip instance to avoid "already in doc" error
            help_text = help_tooltip.content
            tooltip_all = Tooltip(content=help_text, position="right")
            help_btn = HelpButton(
                tooltip=tooltip_all,
                width=20,
                height=20,
                align="center",
                button_type="light",
                stylesheets=[toggle_stylesheet],
            )
            hdr = row(
                toggle_btn,
                module_title_div,
                help_btn,
                sizing_mode="stretch_width",
                align="center",
                styles={"overflow": "hidden"},
            )
        else:
            hdr = row(toggle_btn, module_title_div, sizing_mode="stretch_width", align="center")

        controls_variables_all.append(hdr)

        # Add the module's CheckboxButtonGroup for variables (collapsible, starts folded)
        cbg = widgets["variables_widgets_all"][i]
        cbg.visible = False
        module_contents_all.append(cbg)
        controls_variables_all.append(cbg)

    # Add toggle callbacks for All Samples view
    for i, toggle_btn in enumerate(module_toggles_all):
        content = module_contents_all[i]
        toggle_btn.on_click(make_toggle_callback(toggle_btn, content))

    # Create the two variables section containers
    variables_section_one = column(
        variables_title_one, *controls_variables_one, visible=True, sizing_mode="stretch_width"
    )
    variables_section_all = column(
        variables_title_all, *controls_variables_all, visible=False, sizing_mode="stretch_width"
    )

    return VariableControls(
        one_sample=variables_section_one,
        all_samples=variables_section_all,
        genome_checkbox=genome_cbg_one,
        genome_index=genome_index_one,
    )


def attach_module_synchronization(widgets: Mapping[str, Any], interaction_lock: Mapping[str, bool]) -> None:
    for module_checkbox, variables in zip(widgets["module_widgets_one"], widgets["variables_widgets_one"]):
        lock = {"locked": False}

        def module_changed(attr, old, new, module_checkbox=module_checkbox, variables=variables, lock=lock):
            if lock["locked"] or interaction_lock.get("locked", False):
                return
            lock["locked"] = True
            try:
                variables.active = list(range(len(variables.labels))) if 0 in module_checkbox.active else []
            finally:
                lock["locked"] = False

        def variables_changed(attr, old, new, module_checkbox=module_checkbox, variables=variables, lock=lock):
            if lock["locked"] or interaction_lock.get("locked", False):
                return
            lock["locked"] = True
            try:
                total = len(variables.labels)
                module_checkbox.active = [0] if total > 0 and len(variables.active) == total else []
            finally:
                lock["locked"] = False

        module_checkbox.on_change("active", module_changed)
        variables.on_change("active", variables_changed)
