from types import SimpleNamespace

from thebigbam.plotting.application.apply_inputs import collect_apply_inputs


def _value(value):
    return SimpleNamespace(value=value)


def _color_row(operator, value, *, color="#123456"):
    return {
        "qualifier_select": _value("product"),
        "operator_select": _value(operator),
        "input_ref": {"widget": _value(value)},
        "color_picker": SimpleNamespace(color=color),
    }


def test_collect_apply_inputs_parses_rules_and_stable_display_values():
    widgets = {
        "has_samples": True,
        "contig_select": _value("c1"),
        "sample_select": _value("s1"),
        "variables_widgets_one": ["one"],
        "variables_widgets_all": ["all"],
    }
    bindings = SimpleNamespace(
        widgets=widgets,
        views=SimpleNamespace(active=1),
        feature_type_multichoice=_value(["CDS"]),
        db_path="db.duckdb",
        custom_color_rows=[_color_row("Use random colors", ""), _color_row(">", 3)],
        mag_track_color_rows=[_color_row("Use random colors", ""), _color_row("=", "capsid")],
        apply_annotation_rules_cbg=SimpleNamespace(active=[0]),
        plot_isoforms_cbg=SimpleNamespace(active=[0]),
        feature_label_select=_value("product"),
        mag_track_max_dots_input=_value(200),
        max_genemap_window_input=_value(1000),
        same_y_scale_cbg=SimpleNamespace(active=[0]),
        subplot_height_input=_value(90),
        genemap_height_input=_value(80),
        sequence_height_input=_value(40),
        translated_sequence_height_input=_value(30),
        max_binning_window_input=_value(5000),
        min_coverage_freq_input=_value(0.1),
    )

    result = collect_apply_inputs(bindings)

    assert result.is_all is True
    assert result.active_variables_widgets == ["all"]
    assert result.custom_colors[0]["match_mode"] == "random"
    assert result.custom_colors[1]["value"] == 3.0
    assert len(result.mag_track_colors) == 2
    assert result.genbank_path == "db.duckdb"
