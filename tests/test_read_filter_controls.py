"""Read filtering reuses the existing block editor and settings schema."""

import duckdb
from bokeh.events import ButtonClick
from bokeh.models import Spinner

from thebigbam.plotting.controls.filtering import FilterSectionController
from thebigbam.plotting.controls.read_filters import ReadFilterControls
from thebigbam.plotting.controls.read_plotting import build_read_plot_controls
from thebigbam.plotting.repositories.genome_tracks import GenomeTrackRepository


def _set(row, field, value):
    row["subcategory_select"].value = field
    row["input_ref"]["widget"].value = value


def test_read_filter_rows_switch_numeric_position_and_boolean_inputs():
    controls = ReadFilterControls(circular_available=True)
    assert isinstance(controls.controller, FilterSectionController)
    assert not controls.expression().active
    row = controls.controller.sections[0]["rows"][0]
    assert ("Read % aligned", "Read % aligned") in row["subcategory_select"].options
    _set(row, "MAPQ", 30)
    assert row["comparison_select"].visible
    assert row["comparison_select"].options == [">", "<", "=", "!="]
    _set(row, "Read % aligned", 95.5)
    spinner = row["input_ref"]["widget"]
    assert spinner.step == 0.1 and spinner.low == 0 and spinner.high == 100
    _set(row, "Alignment includes", 123)
    assert row["comparison_select"].options == ["=", "!="]
    _set(row, "Has mate mapped", "FALSE")
    assert not row["comparison_select"].visible
    assert row["input_ref"]["widget"].options == [("", "Choose..."), "TRUE", "FALSE"]
    expression = controls.expression()
    assert expression.sections[0].predicates[0].value == "FALSE"
    assert expression.sections[0].predicates[0].operator == "="


def test_shared_add_remove_controls_and_settings_round_trip():
    controls = ReadFilterControls(circular_available=True)
    controller = controls.controller
    section = controller.sections[0]
    _set(section["rows"][0], "MAPQ", 30)
    section["add_and_btn"].clicks += 1
    _set(section["rows"][1], "Read length", 1000)
    section["rows"][1]["and_div"].value = "OR"
    controller._global_add.clicks += 1
    _set(controller.sections[1]["rows"][0], "Circular", "TRUE")
    controller.inter_section_selects[0].value = "OR"
    saved = controls.serialize()
    restored = ReadFilterControls(circular_available=True)
    restored.restore(saved)
    assert restored.serialize() == saved
    assert restored.expression() == controls.expression()
    restored.controller.sections[1]["rows"][0]["minus_btn"].clicks += 1
    assert len(restored.controller.sections) == 1
    assert restored.controller.count_rows() == 2
    assert restored.controller.inter_section_selects == []


def test_circular_field_only_appears_when_available():
    assert "Circular" not in ReadFilterControls().factory.fields
    assert "Circular" in ReadFilterControls(circular_available=True).factory.fields
    conn = duckdb.connect(":memory:")
    try:
        conn.execute("CREATE TABLE Sample (Sample_name VARCHAR, Circular_mapping BOOLEAN)")
        conn.execute("INSERT INTO Sample VALUES ('linear', false), ('circular', true), ('not_supplied', true)")
        repository = GenomeTrackRepository(conn)
        assert repository.circular_samples(("linear",)) == frozenset()
        assert repository.circular_samples(("linear", "circular")) == frozenset({"circular"})
    finally:
        conn.close()


def test_restoring_unavailable_circular_filter_does_not_create_hidden_predicate(capsys):
    controls = ReadFilterControls()
    controls.restore([{"rows": [{"category": "Reads", "column": "Circular", "operator": "=", "value": "TRUE"}]}])
    assert not controls.expression().active
    assert "unavailable read filter" in capsys.readouterr().out


def test_read_filter_collapse_button_matches_annotation_toggle_and_preserves_rules():
    controls = build_read_plot_controls(
        max_reads=Spinner(value=1000),
        plot_height=Spinner(value=400),
        read_height=Spinner(value=14),
        toggle_stylesheet=":host { color: red; }",
    )
    button, title = controls.filter_header.children
    assert title.text == "Filter alignments to show"
    assert button.width == button.height == 20
    assert button.button_type == "primary" and button.label == "▼"
    assert button.stylesheets == [":host { color: red; }"]
    assert controls.filter_content.visible
    assert not controls.filters.expression().active
    _set(controls.filters.controller.sections[0]["rows"][0], "MAPQ", 30)
    expression = controls.filters.expression()
    button._trigger_event(ButtonClick(model=button))
    assert not controls.filter_content.visible and button.label == "▶"
    button._trigger_event(ButtonClick(model=button))
    assert controls.filter_content.visible and button.label == "▼"
    assert controls.filters.expression() == expression
