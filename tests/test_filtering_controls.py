import panel as pn

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

    section["add_and_btn"].param.trigger("clicks")

    assert controller.count_rows() == 2
    assert changes == [True]
    assert section["rows"][1]["and_div"].value == "AND"
