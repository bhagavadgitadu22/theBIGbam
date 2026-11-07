"""Color templates and rule rows for BAM read alignments."""

from __future__ import annotations

import panel as pn
from bokeh.models import Div, TextInput
from bokeh.models.widgets import ColorPicker, Select, Spinner

from ..models.read_filters import (
    READ_FILTER_CHOICES,
    READ_FILTER_FIELDS,
    READ_FILTER_LABELS,
    read_color_rule_operators,
)
from ..models.read_plot import ReadColorRule
from ..shared.read_colors import (
    DEFAULT_READ_GRADIENT_END,
    DEFAULT_READ_GRADIENT_START,
    read_color_templates,
)
from ..shared.styles import FIELD_LABEL_MARGIN, panel_control_row, panel_stylesheet


def _color_value_widget(field, operator, *, name):
    if field in READ_FILTER_CHOICES:
        if operator in ("has", "has not"):
            return TextInput(name=name, placeholder="Text...", sizing_mode="stretch_width", margin=0)
        choices = READ_FILTER_CHOICES[field]
        return Select(name=name, options=list(choices), value=choices[0], sizing_mode="stretch_width", margin=0)
    if field in ("Non-inward pair", "Has mate mapped", "Circular"):
        return Select(name=name, options=["TRUE", "FALSE"], value="TRUE", sizing_mode="stretch_width", margin=0)
    percentage = field == "Read % aligned"
    return Spinner(
        name=name, value=None, step=0.1 if percentage else 1,
        low=0 if percentage else (1 if field == "Alignment includes" else None),
        high=100 if percentage else None, placeholder="Value...",
        sizing_mode="stretch_width", margin=0,
    )


class ReadColorRuleControls:
    """Manage ordered read-color rules and their built-in templates."""

    def __init__(self, stylesheet="", *, circular_available=False, strand_colors=None):
        self.stylesheet = stylesheet
        self.fields = [field for field in READ_FILTER_FIELDS if field != "Circular" or circular_available]
        self.templates = read_color_templates(strand_colors)
        self.template_select = Select(
            title="Use a color template for reads:",
            options=list(self.templates),
            value="Alignment class",
            sizing_mode="stretch_width",
            margin=FIELD_LABEL_MARGIN,
            css_classes=["sidebar-field"],
        )
        self.color_rules_label = Div(text="Color annotations by:", margin=FIELD_LABEL_MARGIN)
        self.title = Div(text="", align="center")
        self.rows = []
        self.rows_column = pn.Column(
            sizing_mode="stretch_width",
            css_classes=["color-rule-list"],
            stylesheets=[panel_stylesheet(stylesheet)],
        )
        self.add_button = pn.widgets.Button(
            name="+ Add coloring rule",
            margin=(2, 0, 2, 0),
            button_type="success",
            stylesheets=[panel_stylesheet(stylesheet)],
            css_classes=["action-add", "read-color-rule-add"],
        )
        self.rules_content = pn.Column(
            self.add_button,
            self.rows_column,
            sizing_mode="stretch_width",
            css_classes=["nested-section"],
            stylesheets=[panel_stylesheet(stylesheet)],
        )
        self.content = pn.Column(
            self.template_select,
            self.color_rules_label,
            self.rules_content,
            sizing_mode="stretch_width",
            margin=0,
            stylesheets=[panel_stylesheet(stylesheet)],
        )
        self.add_button.on_click(self._add_rule)
        self.template_select.on_change("value", self._template_changed)
        self._replace_rules(self.templates["Alignment class"])

    def _create_rule(self, initial=None):
        initial = dict(initial or {})
        rule_index = len(self.rows)
        initial_field = initial.get("field", self.fields[0])
        if initial_field not in self.fields:
            initial_field = self.fields[0]
        field = Select(
            name=f"benchmark-read-color-rule-{rule_index}-field",
            options=[(value, READ_FILTER_LABELS.get(value, value)) for value in self.fields],
            value=initial_field,
            sizing_mode="stretch_width",
            margin=0,
            styles={"grid-column": "1"},
            css_classes=["responsive-field", "color-qualifier"],
        )
        operator = Select(
            name=f"benchmark-read-color-rule-{rule_index}-operator",
            options=read_color_rule_operators(initial_field),
            value="=",
            sizing_mode="stretch_width",
            margin=0,
            styles={"grid-column": "2"},
            css_classes=["responsive-field", "color-operator"],
        )
        if initial.get("operator") in operator.options:
            operator.value = initial["operator"]
        value = _color_value_widget(
            initial_field, operator.value, name=f"benchmark-read-color-rule-{rule_index}-value"
        )
        value.visible = operator.value not in ("Use random colors", "Use color gradient")
        if initial.get("value") not in (None, ""):
            value.value = initial["value"]
        color = ColorPicker(
            color=initial.get(
                "color",
                DEFAULT_READ_GRADIENT_START if operator.value == "Use color gradient" else "#cccccc",
            ),
            width=42,
            height=30,
            margin=0,
            css_classes=["responsive-action", "color-picker-field"],
        )
        gradient_end = ColorPicker(
            color=initial.get("gradient_end_color", DEFAULT_READ_GRADIENT_END),
            width=42,
            height=30,
            margin=0,
            css_classes=["responsive-action", "color-picker-field", "gradient-color-picker"],
        )
        color_slot = pn.Row(
            color,
            gradient_end,
            sizing_mode="fixed",
            width=88,
            height=30,
            margin=0,
            css_classes=["gradient-color-pair"],
            stylesheets=[panel_stylesheet(self.stylesheet)],
            styles={
                "grid-column": "4", "column-gap": "4px", "gap": "4px", "align-items": "center",
            },
        )
        minus = pn.widgets.Button(
            name="−",
            width=30,
            height=30,
            margin=0,
            stylesheets=[panel_stylesheet(self.stylesheet)],
            css_classes=["responsive-action", "color-remove"],
            styles={"grid-column": "5"},
        )
        value_slot = pn.Row(
            value,
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["color-value-slot"],
            stylesheets=[panel_stylesheet(self.stylesheet)],
            styles={"grid-column": "3"},
        )
        row_widget = panel_control_row(
            field,
            operator,
            value_slot,
            color_slot,
            minus,
            sizing_mode="stretch_width",
            margin=(2, 0, 2, 0),
            css_classes=["responsive-control-row", "coloring-row", "read-coloring-row"],
            stylesheets=[panel_stylesheet(self.stylesheet)],
        )
        row = {"field": field, "operator": operator, "value": value, "value_slot": value_slot,
               "color": color, "gradient_end": gradient_end, "color_slot": color_slot,
               "gradient_initialized": operator.value == "Use color gradient",
               "minus": minus, "row_widget": row_widget}

        def update_mode_visibility(mode):
            random = mode == "Use random colors"
            gradient = mode == "Use color gradient"
            value_slot.visible = not (random or gradient)
            color.visible = not random
            gradient_end.visible = gradient
            color_slot.visible = not random
            color_slot.width = 88 if gradient else 42
            classes = [name for name in row_widget.css_classes if name not in {"read-color-random", "read-color-gradient"}]
            if random:
                classes.append("read-color-random")
            if gradient:
                classes.append("read-color-gradient")
            row_widget.css_classes = classes
            styles = dict(row_widget.styles)
            styles["grid-template-columns"] = (
                "minmax(0, 3fr) minmax(0, 2fr) 0 88px 30px"
                if gradient else
                "minmax(0, 3fr) minmax(0, 5.7fr) 0 0 30px"
                if random else
                "minmax(0, 3fr) minmax(0, 2fr) minmax(0, 2.5fr) 42px 30px"
            )
            row_widget.styles = styles

        update_mode_visibility(operator.value)

        def update_field(_attr, _old, new):
            choices = read_color_rule_operators(new)
            operator.options = choices
            if operator.value not in choices:
                operator.value = choices[0]
            new_value = _color_value_widget(
                new, operator.value, name=f"benchmark-read-color-rule-{rule_index}-value"
            )
            row["value"] = new_value
            value_slot.objects = [new_value]
            update_mode_visibility(operator.value)

        field.on_change("value", update_field)

        def update_operator(_attr, old, new):
            current_value = getattr(row["value"], "value", None)
            new_value = _color_value_widget(
                field.value, new, name=f"benchmark-read-color-rule-{rule_index}-value"
            )
            if new == "Use color gradient" and old != new and not row["gradient_initialized"]:
                color.color = DEFAULT_READ_GRADIENT_START
                gradient_end.color = DEFAULT_READ_GRADIENT_END
                row["gradient_initialized"] = True
            if new not in ("Use random colors", "Use color gradient") and current_value not in (None, ""):
                new_value.value = current_value
            row["value"] = new_value
            value_slot.objects = [new_value]
            update_mode_visibility(new)

        operator.on_change("value", update_operator)
        minus.on_click(lambda _event: self._remove_rule(row))
        return row

    def _replace_rules(self, rules):
        self.rows = []
        for rule in rules:
            self.rows.append(self._create_rule(rule))
        self._rebuild()

    def _template_changed(self, _attr, _old, new):
        self._replace_rules(self.templates[new])

    def _add_rule(self, _event):
        self.rows.append(self._create_rule())
        self._rebuild()

    def _remove_rule(self, row):
        if row in self.rows:
            self.rows.remove(row)
            self._rebuild()

    def _rebuild(self):
        self.rows_column.objects = [row["row_widget"] for row in self.rows]
        count = len(self.rows)
        self.title.text = f"Customise reads ({count} rule{'s' if count != 1 else ''})"

    def rules(self):
        rules = []
        for row in self.rows:
            operator = row["operator"].value
            value = None if operator in ("Use random colors", "Use color gradient") else row["value"].value
            if operator not in ("Use random colors", "Use color gradient") and (value is None or value == ""):
                continue
            rules.append(ReadColorRule(
                row["field"].value,
                operator,
                value,
                row["color"].color,
                row["gradient_end"].color if operator == "Use color gradient" else None,
            ))
        return tuple(rules)

    def serialize(self):
        serialized = []
        for row in self.rows:
            operator = row["operator"].value
            rule = {
                "field": row["field"].value,
                "operator": operator,
                "value": None if operator in ("Use random colors", "Use color gradient") else row["value"].value,
                "color": row["color"].color,
            }
            if operator == "Use color gradient":
                rule["gradient_end_color"] = row["gradient_end"].color
            serialized.append(rule)
        return serialized

    def restore(self, template, rules=None):
        template = template if template in self.templates else "Alignment class"
        self.template_select.value = template
        if rules is None:
            return
        valid = []
        for rule in rules:
            field = rule.get("field", rule.get("qualifier"))
            if field not in self.fields:
                continue
            operator = rule.get("operator", "=")
            value = rule.get("value")
            if operator not in read_color_rule_operators(field):
                continue
            if field == "Strand" and value in (1, "1", "+1", "+"):
                value = "Forward"
            elif field == "Strand" and value in (-1, "-1", "-"):
                value = "Reverse"
            valid.append({"field": field, "operator": operator, "value": value,
                          "color": rule.get("color", "#cccccc"),
                          "gradient_end_color": rule.get("gradient_end_color", DEFAULT_READ_GRADIENT_END)})
        self._replace_rules(valid)
