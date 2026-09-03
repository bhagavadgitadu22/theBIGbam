import json

from bokeh.models import InlineStyleSheet, Spinner

from thebigbam.plotting.controls.color_rules import build_color_rule_controls
from thebigbam.plotting.controls.searchable_select import SearchableSelect
from thebigbam.plotting.repositories.color_templates import ColorTemplateRepository


class MetadataService:
    def distinct_values(self, category, column):
        return ["Defense"] if column == "activity" else []

    def search_distinct_values(self, category, column, query, limit=100):
        return ["Defense"] if column == "activity" and query.lower() in "defense" else []


class Connection:
    def execute(self, sql):
        assert "Color_templates" in sql
        return self

    def fetchall(self):
        return [("default", "product", "has", "kinase", "#123456")]


def test_color_template_repository_groups_rules():
    templates = ColorTemplateRepository(Connection()).load()
    assert templates["default"] == [
        {
            "qualifier_name": "product",
            "operator": "has",
            "qualifier_value": "kinase",
            "color": "#123456",
        }
    ]


def test_color_controls_expose_independent_annotation_and_mag_rule_lists():
    stylesheet = InlineStyleSheet(css=":host { color: red; }")
    controls = build_color_rule_controls(
        MetadataService(),
        {"Annotations": {"columns": {}}},
        {"default": []},
        stylesheet,
        False,
    )

    assert controls.template_select.options == ["(none)", "default"]
    assert controls.custom_rows == []
    assert controls.mag_rows == []
    assert controls.custom_rows is not controls.mag_rows
    assert controls.custom_column.objects[0].name == "+ Add coloring rule"
    assert controls.custom_column.objects[0].stylesheets == [stylesheet.css]
    assert controls.custom_column.objects[1].css_classes == ["color-rule-list"]
    assert controls.custom_title.text == "Customise genomic annotations (0 rules)"
    assert controls.mag_title.text == "Customise MAG track (0 rules)"


def test_color_rule_count_updates_without_moving_add_button():
    controls = build_color_rule_controls(
        MetadataService(),
        {"Annotations": {"columns": {"Start": {"type": "numeric"}}}},
        {"preset": []},
        InlineStyleSheet(css=""),
        False,
    )

    add_button, rows = controls.custom_column.objects
    add_button.clicks += 1
    add_button.clicks += 1

    assert controls.custom_column.objects[0] is add_button
    assert len(rows.objects) == 2
    assert controls.custom_title.text == "Customise genomic annotations (2 rules)"


def test_color_rule_switches_from_numeric_to_text_qualifier_values():
    controls = build_color_rule_controls(
        MetadataService(),
        {
            "Annotations": {
                "columns": {
                    "Start": {"type": "numeric"},
                    "activity": {"type": "text", "source": "Annotation_qualifier"},
                }
            }
        },
        {},
        InlineStyleSheet(css=""),
        False,
    )

    row = controls.create_row()
    assert isinstance(row["input_ref"]["widget"], Spinner)

    row["qualifier_select"].value = "activity"

    value_widget = row["input_ref"]["widget"]
    assert isinstance(value_widget, SearchableSelect)
    assert value_widget.options == ["Defense"]
    assert value_widget.server_search
    assert value_widget.min_search_chars == 0

    value_widget.search_query = "De"
    value_widget.search_request = json.dumps({"nonce": 1, "query": "De"})
    assert value_widget.options == ["Defense"]


def test_color_rule_fields_use_shared_proportional_layout():
    controls = build_color_rule_controls(
        MetadataService(),
        {"Annotations": {"columns": {"Start": {"type": "numeric"}}}},
        {"preset": []},
        InlineStyleSheet(css=""),
        False,
    )

    row = controls.create_row()

    assert row["qualifier_select"].css_classes == ["responsive-field", "color-qualifier"]
    assert row["operator_select"].css_classes == ["responsive-field", "color-operator"]
    assert row["input_ref"]["widget"].css_classes == ["responsive-field", "color-value"]
    assert row["color_picker"].css_classes == ["responsive-action", "color-picker-field"]
    assert row["minus_btn"].css_classes == ["responsive-action", "color-remove"]
    assert row["row_widget"].css_classes == [
        "control-row",
        "responsive-control-row",
        "coloring-row",
    ]
    assert controls.template_select.css_classes == ["sidebar-field"]
    assert controls.template_select.margin == 0
    assert controls.feature_label_select.css_classes == ["sidebar-field"]
    assert controls.feature_label_select.margin == 0


def test_color_has_rule_uses_free_text_autocomplete_with_suggestions():
    controls = build_color_rule_controls(
        MetadataService(),
        {"Annotations": {"columns": {"activity": {"type": "text"}}}},
        {},
        InlineStyleSheet(css=""),
        False,
    )
    row = controls.create_row()

    row["operator_select"].value = "has"
    widget = row["input_ref"]["widget"]

    assert isinstance(widget, SearchableSelect)
    assert widget.allow_custom
    assert widget.options == ["Defense"]
