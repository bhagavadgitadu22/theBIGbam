from types import SimpleNamespace

from thebigbam.plotting.settings.collection import SettingsCollector, SettingsCollectorBindings
from thebigbam.plotting.settings.restoration import SettingsRestoreBindings, restore_settings


def _control(*, value=None, active=None, options=None, labels=None, color=None):
    return SimpleNamespace(
        value=value,
        active=[] if active is None else active,
        options=[] if options is None else options,
        labels=[] if labels is None else labels,
        color=color,
    )


def _filter_row(category, column, operator, value, connector=None):
    return {
        "category_select": _control(value=category),
        "subcategory_select": _control(value=column),
        "comparison_select": _control(value=operator),
        "input_ref": {"widget": _control(value=value)},
        "and_div": None if connector is None else _control(value=connector),
    }


def _color_row(qualifier, operator, value, color):
    return {
        "qualifier_select": _control(value=qualifier),
        "operator_select": _control(value=operator, options=["Equals", "Contains", "Use random colors"]),
        "input_ref": {"widget": _control(value=value)},
        "color_picker": _control(color=color),
    }


def _state(*, populated):
    suffix = "saved" if populated else "default"
    widgets = {
        "view_radio": _control(active=0 if populated else 1),
        "sample_select": _control(value=f"sample-{suffix}", options=["sample-default", "sample-saved"]),
        "contig_select": _control(value=f"contig-{suffix}", options=["contig-default", "contig-saved"]),
        "mag_select": _control(value=f"mag-{suffix}", options=["mag-default", "mag-saved"]),
        "module_names": ["Coverage", "Variants"],
        "module_widgets_one": [_control(active=[0] if populated else []), _control(active=[0])],
        "variables_widgets_one": [
            _control(labels=["depth", "breadth"], active=[1] if populated else [0]),
            _control(labels=["snps", "indels"], active=[0, 1] if populated else []),
        ],
        "variables_widgets_all": [
            _control(labels=["depth", "breadth"], active=[0, 1] if populated else []),
            _control(labels=["snps", "indels"], active=[1] if populated else [0]),
        ],
    }
    controls = {
        "views": _control(active=1 if populated else 0),
        "from_position_input": _control(value="12" if populated else "1"),
        "to_position_input": _control(value="345" if populated else ""),
        "feature_type_multichoice": _control(
            value=["CDS", "rRNA"] if populated else [], options=["CDS", "rRNA"]
        ),
        "plot_isoforms_cbg": _control(active=[0] if populated else []),
        "combined_features_cbg": _control(
            labels=["Genes", "Domains"], active=[0, 1] if populated else []
        ),
        "sequence_cbg": _control(active=[0] if populated else []),
        "translated_sequence_cbg": _control(active=[0] if populated else []),
        "genome_master_cbg": _control(active=[0] if populated else []),
        "feature_label_select": _control(
            value="product" if populated else "locus_tag", options=["locus_tag", "product"]
        ),
        "apply_annotation_rules_cbg": _control(active=[0] if populated else []),
        "min_coverage_freq_input": _control(value=0.37 if populated else 0.0),
        "max_genemap_window_input": _control(value=1111 if populated else 1000),
        "max_sequence_window_input": _control(value=2222 if populated else 2000),
        "max_binning_window_input": _control(value=3333 if populated else 3000),
        "genemap_height_input": _control(value=144 if populated else 100),
        "sequence_height_input": _control(value=155 if populated else 100),
        "translated_sequence_height_input": _control(value=166 if populated else 100),
        "subplot_height_input": _control(value=177 if populated else 100),
        "mag_params_category_select": _control(
            value="MAG" if populated else "Sample", options=["Sample", "MAG"]
        ),
        "mag_params_metric_select": _control(
            value="quality" if populated else "name", options=["name", "quality"]
        ),
        "mag_params_direction": _control(active=1 if populated else 0),
        "mag_params_sort_sample_select": _control(
            value="sample-saved" if populated else "", options=["", "sample-saved"]
        ),
        "mag_track_max_dots_input": _control(value=4321 if populated else 1000),
        "max_samples_input": _control(value=42 if populated else 10),
        "sample_order_category_select": _control(
            value="Coverage" if populated else "Sample", options=["Sample", "Coverage"]
        ),
        "sample_order_metric_select": _control(
            value="mean" if populated else "name", options=["name", "mean"]
        ),
        "sample_order_direction": _control(active=1 if populated else 0),
        "same_y_scale_cbg": _control(active=[0] if populated else []),
    }
    if populated:
        sections = [
            {"rows": [_filter_row("Contig", "Length", ">", 1000), _filter_row("Contig", "GC", "<", 0.7, "OR")]},
            {"rows": [_filter_row("Sample", "Site", "=", "ocean")]},
        ]
        section_connectors = [_control(value="AND")]
        custom_colors = [_color_row("product", "Contains", "kinase", "#112233")]
        mag_colors = [_color_row("taxonomy", "Equals", "Bacteria", "#445566")]
    else:
        sections = [{"rows": []}]
        section_connectors = []
        custom_colors = []
        mag_colors = []
    return widgets, controls, sections, section_connectors, custom_colors, mag_colors


def _collector(widgets, controls, sections, connectors, custom_colors, mag_colors):
    return SettingsCollector(
        SettingsCollectorBindings(
            db_path="example.db",
            widgets=widgets,
            or_sections=sections,
            inter_section_selects=connectors,
            custom_color_rows=custom_colors,
            mag_track_color_rows=mag_colors,
            **controls,
        )
    )


def test_every_persisted_left_panel_setting_round_trips_through_restore():
    source = _state(populated=True)
    expected = _collector(*source).collect()
    target = _state(populated=False)
    widgets, controls, sections, connectors, custom_colors, mag_colors = target

    def create_section():
        return {"rows": []}

    def create_query_row(section, initial_category, initial_column, initial_operator):
        return _filter_row(
            initial_category,
            initial_column,
            initial_operator,
            None,
            "AND" if section["rows"] else None,
        )

    def rebuild_filtering():
        connectors[:] = [_control(value="AND") for _ in range(max(0, len(sections) - 1))]

    def create_color_row(_rows, _rebuild):
        return _color_row("product", "Equals", "", "#000000")

    restore_settings(
        expected,
        SettingsRestoreBindings(
            widgets=widgets,
            or_sections=sections,
            inter_section_selects=connectors,
            custom_color_rows=custom_colors,
            mag_track_color_rows=mag_colors,
            filtering_metadata={
                "Contig": {"columns": {"Length": {}, "GC": {}}},
                "Sample": {"columns": {"Site": {}}},
            },
            color_qualifier_options=["product", "taxonomy"],
            create_or_section=create_section,
            create_query_row=create_query_row,
            rebuild_section=lambda _section: None,
            rebuild_filtering_content=rebuild_filtering,
            create_color_row=create_color_row,
            rebuild_color_rows=lambda: None,
            rebuild_mag_track_color_rows=lambda: None,
            **controls,
        ),
    )

    actual = _collector(*target).collect()
    expected["_meta"].pop("saved_at")
    actual["_meta"].pop("saved_at")
    assert actual == expected
