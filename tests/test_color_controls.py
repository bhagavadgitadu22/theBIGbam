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
    assert controls.custom_column.objects[-1].name == "+ Add coloring rule"
    assert controls.custom_column.objects[-1].stylesheets == [stylesheet.css]


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
    assert value_widget.options == []
    assert value_widget.server_search
    assert value_widget.min_search_chars == 2

    value_widget.search_query = "De"
    value_widget.search_nonce += 1
    assert value_widget.options == ["Defense"]
