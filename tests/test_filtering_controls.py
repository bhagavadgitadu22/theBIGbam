import panel as pn
from bokeh.models import InlineStyleSheet

from thebigbam.plotting.controls.filtering import FilterSectionController


def _row(section):
    return {
        "row_wrapper": pn.Row(),
        "minus_btn": pn.widgets.Button(name="-"),
        "and_div": None,
    }


def test_filter_section_controller_owns_sections_and_connectors():
    changes = []
    controller = FilterSectionController(_row, lambda: changes.append(True), "", "")

    assert controller.count_rows() == 1
    second = controller.create_section()
    controller.replace_sections([controller.sections[0], second])

    assert controller.count_rows() == 2
    assert len(controller.inter_section_selects) == 1
    assert controller.inter_section_selects[0].value == "AND"
    assert controller.content.objects[-1].name == "+ Add AND/OR"


def test_filter_section_controller_add_row_refreshes_once():
    changes = []
    controller = FilterSectionController(_row, lambda: changes.append(True), "", "")
    section = controller.sections[0]
    first_wrapper = section["rows"][0]["row_wrapper"]

    section["add_and_btn"].param.trigger("clicks")

    assert controller.count_rows() == 2
    assert changes == [True]
    assert section["rows"][1]["and_div"].value == "AND"
    assert section["column"].objects[0] is first_wrapper


def test_filter_section_converts_bokeh_stylesheets_for_panel_widgets():
    stylesheet = InlineStyleSheet(css=":host { color: red; }")

    controller = FilterSectionController(_row, lambda: None, stylesheet, stylesheet)

    assert controller._global_add.stylesheets == [stylesheet.css]
    assert controller.sections[0]["add_and_btn"].stylesheets == [stylesheet.css]
    assert controller._global_add.button_type == "success"
    assert controller.sections[0]["add_and_btn"].button_type == "success"
    assert controller._global_add.css_classes == ["action-add"]
    assert controller.sections[0]["column"].css_classes == ["nested-section"]


def test_filter_action_controls_share_vertical_geometry():
    controller = FilterSectionController(_row, lambda: None, "", "")
    section = controller.sections[0]
    action_row = section["column"].objects[1]

    assert action_row.css_classes == ["control-row", "nested-control-row", "filter-action-row"]
    assert action_row.height == 30
    assert section["add_and_btn"].height == 30
    assert section["add_and_btn"].margin == 0

    section["add_and_btn"].param.trigger("clicks")
    connector = section["rows"][1]["and_div"]
    assert connector.height == 30
    assert connector.margin == 0
