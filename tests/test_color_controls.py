from thebigbam.plotting.controls.color_rules import build_color_rule_controls
from thebigbam.plotting.repositories.color_templates import ColorTemplateRepository


class MetadataService:
    def distinct_values(self, category, column):
        return []


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
    controls = build_color_rule_controls(
        MetadataService(),
        {"Annotations": {"columns": {}}},
        {"default": []},
        "",
        False,
    )

    assert controls.template_select.options == ["(none)", "default"]
    assert controls.custom_rows == []
    assert controls.mag_rows == []
    assert controls.custom_rows is not controls.mag_rows
    assert controls.custom_column.objects[-1].name == "+ Add coloring rule"
