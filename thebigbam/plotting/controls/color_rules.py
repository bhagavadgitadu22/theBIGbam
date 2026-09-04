"""Dynamic annotation and MAG color-rule controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import panel as pn
from bokeh.models import Div
from bokeh.models.widgets import ColorPicker, Select, Spinner

from ..shared.styles import panel_control_row, panel_stylesheet
from .distinct_value_select import build_distinct_value_select
from .searchable_select import SearchableSelect


@dataclass(frozen=True)
class ColorRuleControls:
    template_select: Any | None
    qualifier_options: list[str]
    custom_rows: list[dict[str, Any]]
    custom_column: Any
    custom_title: Any
    mag_rows: list[dict[str, Any]]
    mag_column: Any
    mag_title: Any
    feature_label_select: Any | None
    create_row: Any
    rebuild_custom: Any
    rebuild_mag: Any


def build_color_rule_controls(
    metadata_service: Any,
    filtering_metadata: Mapping[str, Any],
    color_templates: Mapping[str, list[dict[str, Any]]],
    stylesheet: Any,
    enable_timing: bool,
) -> ColorRuleControls:
    template_select = None
    if color_templates:
        template_select = Select(
            title="Use a color template:",
            value="(none)",
            options=["(none)", *color_templates],
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["sidebar-field"],
        )

    # Build qualifier key options for custom color rows — reuse the exact same
    # list the Filtering section iterates over (text + numeric alike).
    annotation_meta = filtering_metadata.get("Annotations", {}).get("columns", {})
    color_qualifier_options = list(annotation_meta.keys())

    # Custom color row system (Panel widgets, same pattern as filtering rows)
    custom_color_rows = []
    custom_title = Div(text="Customise genomic annotations (0 rules)", align="center")
    add_color_btn = pn.widgets.Button(
        name="+ Add coloring rule",
        margin=(2, 0, 2, 0),
        button_type="success",
        stylesheets=[panel_stylesheet(stylesheet)],
        css_classes=["action-add", "benchmark-add-annotation-rule"],
    )
    custom_rows_column = pn.Column(
        sizing_mode="stretch_width",
        css_classes=["color-rule-list"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )
    custom_color_column = pn.Column(
        add_color_btn,
        custom_rows_column,
        sizing_mode="stretch_width",
        css_classes=["nested-section"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )

    TEXT_OPS = ["=", "!=", "has", "has not", "Use random colors"]
    NUMERIC_OPS = ["=", ">", "<", "!=", "Use random colors"]

    def _build_color_value_widget(is_text, operator, *, category="Annotations", column="", name):
        """Return the value widget that matches (type, operator).

        Widgets live directly inside the row (no wrapper container), so each
        one carries its own stretch_width sizing and the 2px right margin
        that separates it from the color picker.
        """
        if is_text:
            # Substring operators retain server suggestions while accepting a
            # literal value that isn't one of the complete distinct values.
            return build_distinct_value_select(
                metadata_service,
                category,
                column,
                name=name,
                allow_custom=operator in ("has", "has not"),
                sizing_mode="stretch_width",
                margin=0,
                css_classes=["responsive-field", "color-value"],
            ), True
        # Numeric columns always get a Spinner regardless of operator.
        return Spinner(
            name=name,
            value=0,
            placeholder="Value...",
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["responsive-field", "color-value"],
        ), False

    def create_color_row(target_rows=None, rebuild_fn=None):
        """Create a single custom color row with qualifier / operator / value / color / remove widgets.

        target_rows / rebuild_fn: the list and rebuild function to use for the
        remove button. Default to custom_color_rows / rebuild_color_rows so
        existing callers that omit these args continue to work.

        Mirrors create_query_row(): text columns get SearchableSelect and
        numeric columns get Spinner. Both types
        gain a 'Use random colors' operator that hides the value container and
        the color picker.
        """
        target = target_rows if target_rows is not None else custom_color_rows
        rule_kind = "mag" if target is mag_track_color_rows else "annotation"
        rule_index = len(target)
        initial_key = color_qualifier_options[0] if color_qualifier_options else ""
        initial_info = annotation_meta.get(initial_key, {})
        initial_is_text = initial_info.get("type") == "text"

        qualifier_select = Select(
            name=f"benchmark-{rule_kind}-rule-{rule_index}-qualifier",
            options=[(k, k.replace("_", " ").replace("percentage", "(%)")) for k in color_qualifier_options],
            value=initial_key,
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["responsive-field", "color-qualifier"],
        )

        operator_select = Select(
            name=f"benchmark-{rule_kind}-rule-{rule_index}-operator",
            options=TEXT_OPS if initial_is_text else NUMERIC_OPS,
            value="=",
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["responsive-field", "color-operator"],
        )

        # Dynamic value widget sits directly at index 2 of the Row — swapped
        # in place when the qualifier type or operator changes. No wrapper
        # pn.Column, because stretch_width on a wrapper introduces extra
        # vertical/horizontal padding that pushed the widget out of line.
        input_name = f"benchmark-{rule_kind}-rule-{rule_index}-value"
        initial_input, initial_is_panel = _build_color_value_widget(
            initial_is_text, "=", column=initial_key, name=input_name
        )
        current_input_ref = {"widget": initial_input, "is_panel": initial_is_panel}

        color_picker = ColorPicker(
            color="#cccccc",
            width=42,
            height=30,
            margin=0,
            css_classes=["responsive-action", "color-picker-field"],
        )

        minus_btn = pn.widgets.Button(
            name="\u2212",
            width=30,
            height=30,
            margin=0,
            stylesheets=[panel_stylesheet(stylesheet)],
            css_classes=["responsive-action", "color-remove"],
        )

        row_widget = panel_control_row(
            qualifier_select,
            operator_select,
            initial_input,
            color_picker,
            minus_btn,
            sizing_mode="stretch_width",
            margin=(2, 0, 2, 0),
            css_classes=[
                "responsive-control-row",
                "coloring-row",
            ],
            stylesheets=[panel_stylesheet(stylesheet)],
        )

        VALUE_IDX = 2  # position of the value widget inside row_widget

        row_data = {
            "qualifier_select": qualifier_select,
            "operator_select": operator_select,
            "input_ref": current_input_ref,
            "color_picker": color_picker,
            "minus_btn": minus_btn,
            "row_widget": row_widget,
        }

        def _swap_value_widget(new_widget, is_panel):
            row_widget[VALUE_IDX] = new_widget
            current_input_ref["widget"] = new_widget
            current_input_ref["is_panel"] = is_panel

        def _apply_random_visibility():
            use_random = operator_select.value == "Use random colors"
            current_input_ref["widget"].visible = not use_random
            color_picker.visible = not use_random
            # Responsive flex sizing lets the remaining fields absorb the
            # space released by hidden random-color controls.

        def update_color_input_widget(col_name):
            """Rebuild operator options and the value widget for the new column type."""
            col_info = annotation_meta.get(col_name, {})
            is_text = col_info.get("type") == "text"

            # Swap operator options for the new type; preserve current operator
            # if still valid, otherwise default back to "=".
            new_ops = TEXT_OPS if is_text else NUMERIC_OPS
            operator_select.options = new_ops
            if operator_select.value not in new_ops:
                operator_select.value = "="

            new_widget, is_panel = _build_color_value_widget(
                is_text, operator_select.value, column=col_name, name=input_name
            )
            _swap_value_widget(new_widget, is_panel)

            # Respect random-mode visibility even after a qualifier swap.
            _apply_random_visibility()

        def on_qualifier_change(attr, old, new):
            update_color_input_widget(new)

        qualifier_select.on_change("value", on_qualifier_change)

        def on_operator_change(attr, old, new):
            col_info = annotation_meta.get(qualifier_select.value, {})
            is_text = col_info.get("type") == "text"
            use_random = new == "Use random colors"

            # Rebuild only when the selection/free-text mode changes, and keep
            # the current value across that mode transition.
            if is_text and not use_random:
                allow_custom = new in ("has", "has not")
                cur = current_input_ref["widget"]
                if not isinstance(cur, SearchableSelect) or cur.allow_custom != allow_custom:
                    new_input, is_panel = _build_color_value_widget(
                        True, new, column=qualifier_select.value, name=input_name
                    )
                    new_input.value = cur.value or ""
                    _swap_value_widget(new_input, is_panel)

            _apply_random_visibility()

        operator_select.on_change("value", on_operator_change)

        _target = target_rows if target_rows is not None else custom_color_rows
        _rebuild = rebuild_fn if rebuild_fn is not None else rebuild_color_rows

        def remove_row_callback(event):
            if row_data in _target:
                _target.remove(row_data)
                _rebuild()

        minus_btn.on_click(remove_row_callback)
        return row_data

    def rebuild_color_rows():
        children = [rd["row_widget"] for rd in custom_color_rows]
        custom_rows_column.objects = children
        count = len(children)
        custom_title.text = f"Customise genomic annotations ({count} rule{'s' if count != 1 else ''})"

    def add_color_callback(event):
        new_row = create_color_row()
        custom_color_rows.append(new_row)
        rebuild_color_rows()

    add_color_btn.on_click(add_color_callback)

    # --- MAG track coloring rules (same widget pattern, separate list) ---
    mag_track_color_rows = []
    mag_title = Div(text="Customise MAG track (0 rules)", align="center")
    add_mag_track_btn = pn.widgets.Button(
        name="+ Add coloring rule",
        margin=(2, 0, 2, 0),
        button_type="success",
        stylesheets=[panel_stylesheet(stylesheet)],
        css_classes=["action-add", "benchmark-add-mag-rule"],
    )
    mag_rows_column = pn.Column(
        sizing_mode="stretch_width",
        css_classes=["color-rule-list"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )
    mag_track_color_column = pn.Column(
        add_mag_track_btn,
        mag_rows_column,
        sizing_mode="stretch_width",
        css_classes=["nested-section"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )

    def rebuild_mag_track_color_rows():
        children = [rd["row_widget"] for rd in mag_track_color_rows]
        mag_rows_column.objects = children
        count = len(children)
        mag_title.text = f"Customise MAG track ({count} rule{'s' if count != 1 else ''})"

    def add_mag_track_color_callback(event):
        new_row = create_color_row(mag_track_color_rows, rebuild_mag_track_color_rows)
        mag_track_color_rows.append(new_row)
        rebuild_mag_track_color_rows()

    add_mag_track_btn.on_click(add_mag_track_color_callback)

    # Template selection callback - populates color rows from template rules
    if template_select is not None:

        def on_template_change(attr, old, new):
            custom_color_rows.clear()
            if new != "(none)" and new in color_templates:
                for rule in color_templates[new]:
                    row_data = create_color_row()
                    # Order matters: qualifier first (triggers widget rebuild),
                    # then operator (may switch autocomplete free-text mode),
                    # then the value on whichever widget is now current.
                    row_data["qualifier_select"].value = rule["qualifier_name"]
                    row_data["operator_select"].value = rule["operator"]
                    if rule["operator"] != "Use random colors":
                        widget = row_data["input_ref"]["widget"]
                        if isinstance(widget, SearchableSelect):
                            # Seed values (e.g. "DNA, RNA and nucleotide metabolism")
                            # may not be in the current distinct list; append so
                            # the Select can display them.
                            opts = list(widget.options)
                            if rule["qualifier_value"] not in opts:
                                widget.options = opts + [rule["qualifier_value"]]
                        widget.value = rule["qualifier_value"]
                        row_data["color_picker"].color = rule["color"]
                    custom_color_rows.append(row_data)
            rebuild_color_rows()

        template_select.on_change("value", on_template_change)

    # Feature label dropdown — same qualifier list the Filtering Annotations
    # category and the coloring rules use, so every annotation attribute (KV
    # keys + direct Contig_annotation columns) can be picked as the tooltip.
    feature_label_select = None
    label_keys = color_qualifier_options
    if label_keys:
        feature_label_select = Select(
            title="Label annotations with:",
            value="product" if "product" in label_keys else label_keys[0],
            options=label_keys,
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["sidebar-field"],
        )

    return ColorRuleControls(
        template_select=template_select,
        qualifier_options=color_qualifier_options,
        custom_rows=custom_color_rows,
        custom_column=custom_color_column,
        custom_title=custom_title,
        mag_rows=mag_track_color_rows,
        mag_column=mag_track_color_column,
        mag_title=mag_title,
        feature_label_select=feature_label_select,
        create_row=create_color_row,
        rebuild_custom=rebuild_color_rows,
        rebuild_mag=rebuild_mag_track_color_rows,
    )
