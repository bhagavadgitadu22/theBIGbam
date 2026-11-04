"""Three-field read filter rows using the standard Filtering action rows."""

from types import SimpleNamespace

import panel as pn
from bokeh.models import Select, Spinner

from ..models.read_filters import BOOLEAN_READ_FILTERS, READ_FILTER_LABELS, read_filter_operators
from ..shared.styles import panel_control_row, panel_stylesheet


class ReadFilterRowFactory:
    def __init__(self, fields, stylesheet):
        self.fields = list(fields)
        self.stylesheet = stylesheet
        self.controller = None

    def create_row(self, section, initial_column=None, initial_operator=None):
        field = Select(options=[(key, READ_FILTER_LABELS.get(key, key)) for key in self.fields],
                       value=initial_column or self.fields[0], sizing_mode="stretch_width",
                       margin=0, css_classes=["responsive-field", "filter-metric"])
        operator = Select(options=read_filter_operators(field.value), sizing_mode="stretch_width", margin=0,
                          css_classes=["responsive-field", "filter-operator"])
        operator.value = initial_operator if initial_operator in operator.options else operator.options[0]
        container = pn.Column(sizing_mode="stretch_width", margin=0, css_classes=["responsive-field", "filter-value"])
        controls = panel_control_row(field, operator, container, sizing_mode="stretch_width", margin=(3, 0, 3, 0),
                                      css_classes=["filtering-row", "read-filter-row"],
                                      stylesheets=[panel_stylesheet(self.stylesheet)])
        minus = pn.widgets.Button(name="−", width=30, height=30, margin=0,
                                  stylesheets=[panel_stylesheet(self.stylesheet)], css_classes=["responsive-action", "filter-remove"])
        row = {"category_select": SimpleNamespace(value="Reads"), "subcategory_select": field,
               "comparison_select": operator, "input_ref": {"widget": None, "is_panel": False},
               "row_wrapper": controls, "minus_btn": minus, "and_div": None}

        def update(_attr, _old, _new):
            boolean = field.value in BOOLEAN_READ_FILTERS
            operator.options = read_filter_operators(field.value)
            if operator.value not in operator.options:
                operator.value = operator.options[0]
            operator.visible = not boolean
            controls.css_classes = ["control-row", "filtering-row", "read-filter-row", *(["read-filter-boolean"] if boolean else [])]
            if boolean:
                widget = Select(options=["TRUE", "FALSE"], value="TRUE", sizing_mode="stretch_width", margin=0)
            else:
                percentage = field.value == "Alignment percentage"
                widget = Spinner(value=None, step=0.1 if percentage else 1, low=0 if percentage else None,
                                 high=100 if percentage else None, placeholder="Value...", sizing_mode="stretch_width", margin=0)
            container.objects = [widget]
            row["input_ref"]["widget"] = widget

        field.on_change("value", update)
        minus.on_click(lambda event: self.controller.remove_row(row))
        update("initial", None, None)
        return row
