"""Tolerant restoration of saved plotting settings into live controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from ..controls.searchable_select import SearchableSelect


@dataclass(frozen=True)
class SettingsRestoreBindings:
    apply_annotation_rules_cbg: Any
    color_qualifier_options: Any
    combined_features_cbg: Any
    create_color_row: Any
    create_or_section: Any
    create_query_row: Any
    custom_color_rows: Any
    feature_label_select: Any
    feature_type_multichoice: Any
    filtering_metadata: Any
    from_position_input: Any
    genemap_height_input: Any
    genome_master_cbg: Any
    inter_section_selects: Any
    mag_params_category_select: Any
    mag_params_direction: Any
    mag_params_metric_select: Any
    mag_params_sort_sample_select: Any
    mag_track_color_rows: Any
    mag_track_max_dots_input: Any
    max_binning_window_input: Any
    max_genemap_window_input: Any
    max_samples_input: Any
    max_sequence_window_input: Any
    min_coverage_freq_input: Any
    or_sections: Any
    plot_isoforms_cbg: Any
    rebuild_color_rows: Any
    rebuild_filtering_content: Any
    rebuild_mag_track_color_rows: Any
    rebuild_section: Any
    same_y_scale_cbg: Any
    sample_order_category_select: Any
    sample_order_direction: Any
    sample_order_metric_select: Any
    sequence_cbg: Any
    sequence_height_input: Any
    subplot_height_input: Any
    to_position_input: Any
    translated_sequence_cbg: Any
    translated_sequence_height_input: Any
    views: Any
    widgets: Any


def _safe_restore(label, function) -> None:
    try:
        function()
    except Exception as error:
        print(f"[settings] Warning: could not restore '{label}': {error}", flush=True)


def _warn_skipped(label, reason) -> None:
    print(f"[settings] Warning: skipped '{label}' — {reason}", flush=True)


def _restore_genome_features(contig_settings, bindings) -> None:
    saved = contig_settings.get("feature_widgets") or {}
    feature_types = bindings.feature_type_multichoice
    if feature_types is not None and saved.get("feature_types") is not None:
        valid = [value for value in saved["feature_types"] if value in feature_types.options]
        for value in (value for value in saved["feature_types"] if value not in feature_types.options):
            _warn_skipped("contig.feature_widgets.feature_types", f"'{value}' not found in this database")
        _safe_restore("contig.feature_widgets.feature_types", lambda: setattr(feature_types, "value", valid))
    toggles = (
        ("plot_isoforms", bindings.plot_isoforms_cbg),
        ("sequence", bindings.sequence_cbg),
        ("translated_sequence", bindings.translated_sequence_cbg),
        ("genome_master", bindings.genome_master_cbg),
    )
    for key, control in toggles:
        if control is not None and saved.get(key) is not None:
            _safe_restore(
                f"contig.feature_widgets.{key}",
                lambda control=control, key=key: setattr(control, "active", [0] if saved[key] else []),
            )
    combined = bindings.combined_features_cbg
    if combined is not None and saved.get("combined_features") is not None:
        indices = [combined.labels.index(label) for label in saved["combined_features"] if label in combined.labels]
        for label in (label for label in saved["combined_features"] if label not in combined.labels):
            _warn_skipped("contig.feature_widgets.combined_features", f"'{label}' not found in this database")
        _safe_restore("contig.feature_widgets.combined_features", lambda: setattr(combined, "active", indices))
    label_control = bindings.feature_label_select
    label = saved.get("feature_label")
    if label_control is not None and label:
        if label in label_control.options:
            _safe_restore("contig.feature_widgets.feature_label", lambda: setattr(label_control, "value", label))
        else:
            _warn_skipped("contig.feature_widgets.feature_label", f"'{label}' not found in this database")


def restore_settings(settings: Mapping[str, Any], bindings: SettingsRestoreBindings) -> None:
    """Restore compatible values and warn about fields unavailable in this database."""
    apply_annotation_rules_cbg = bindings.apply_annotation_rules_cbg
    color_qualifier_options = bindings.color_qualifier_options
    create_color_row = bindings.create_color_row
    create_or_section = bindings.create_or_section
    create_query_row = bindings.create_query_row
    custom_color_rows = bindings.custom_color_rows
    filtering_metadata = bindings.filtering_metadata
    from_position_input = bindings.from_position_input
    genemap_height_input = bindings.genemap_height_input
    inter_section_selects = bindings.inter_section_selects
    mag_params_category_select = bindings.mag_params_category_select
    mag_params_direction = bindings.mag_params_direction
    mag_params_metric_select = bindings.mag_params_metric_select
    mag_params_sort_sample_select = bindings.mag_params_sort_sample_select
    mag_track_color_rows = bindings.mag_track_color_rows
    mag_track_max_dots_input = bindings.mag_track_max_dots_input
    max_binning_window_input = bindings.max_binning_window_input
    max_genemap_window_input = bindings.max_genemap_window_input
    max_samples_input = bindings.max_samples_input
    max_sequence_window_input = bindings.max_sequence_window_input
    min_coverage_freq_input = bindings.min_coverage_freq_input
    or_sections = bindings.or_sections
    rebuild_color_rows = bindings.rebuild_color_rows
    rebuild_filtering_content = bindings.rebuild_filtering_content
    rebuild_mag_track_color_rows = bindings.rebuild_mag_track_color_rows
    rebuild_section = bindings.rebuild_section
    same_y_scale_cbg = bindings.same_y_scale_cbg
    sample_order_category_select = bindings.sample_order_category_select
    sample_order_direction = bindings.sample_order_direction
    sample_order_metric_select = bindings.sample_order_metric_select
    sequence_height_input = bindings.sequence_height_input
    subplot_height_input = bindings.subplot_height_input
    to_position_input = bindings.to_position_input
    translated_sequence_height_input = bindings.translated_sequence_height_input
    views = bindings.views
    widgets = bindings.widgets

    _restore = _safe_restore
    _warn = _warn_skipped

    # 1. View mode
    view_mode = settings.get("view_mode") or {}
    if "mag_or_contig" in view_mode:
        _restore(
            "view_mode.mag_or_contig", lambda: setattr(widgets["view_radio"], "active", view_mode["mag_or_contig"])
        )
    if "one_or_all_samples" in view_mode:
        _restore("view_mode.one_or_all_samples", lambda: setattr(views, "active", view_mode["one_or_all_samples"]))

    # 2. Sample/contig/MAG selection
    selection = settings.get("selection") or {}
    for key, widget in (
        ("sample", widgets["sample_select"]),
        ("contig", widgets["contig_select"]),
        ("mag", widgets["mag_select"]),
    ):
        if key not in selection:
            continue
        val = selection[key]
        if not val:
            _restore(f"selection.{key}", lambda widget=widget: setattr(widget, "value", ""))
            continue
        if val in widget.options:
            _restore(f"selection.{key}", lambda widget=widget, val=val: setattr(widget, "value", val))
        else:
            _warn(f"selection.{key}", f"'{val}' not found in this database")

    # 3. Position range
    contig_settings = settings.get("contig") or {}
    position_range = contig_settings.get("position_range") or {}
    if "from" in position_range:
        _restore("contig.position_range.from", lambda: setattr(from_position_input, "value", position_range["from"]))
    if "to" in position_range:
        _restore("contig.position_range.to", lambda: setattr(to_position_input, "value", position_range["to"]))

    # 4. Module / variable checkboxes (matched by name, not index)
    for mod_name, mod_settings in (settings.get("variables") or {}).items():
        if mod_name not in widgets["module_names"]:
            _warn(f"variables.{mod_name}", "module not present in this database")
            continue
        i = widgets["module_names"].index(mod_name)
        module_cbg = widgets["module_widgets_one"][i]
        if "module_enabled" in mod_settings:
            _restore(
                f"variables.{mod_name}.module_enabled",
                lambda module_cbg=module_cbg, enabled=bool(mod_settings["module_enabled"]): setattr(
                    module_cbg, "active", [0] if enabled else []
                ),
            )
        for scope_key, cbg in (
            ("selected_one", widgets["variables_widgets_one"][i]),
            ("selected_all", widgets["variables_widgets_all"][i]),
        ):
            if scope_key not in mod_settings:
                continue
            saved_labels = mod_settings[scope_key] or []
            indices = [cbg.labels.index(lbl) for lbl in saved_labels if lbl in cbg.labels]
            missing = [lbl for lbl in saved_labels if lbl not in cbg.labels]
            for lbl in missing:
                _warn(f"variables.{mod_name}.{scope_key}", f"variable '{lbl}' not found")
            _restore(
                f"variables.{mod_name}.{scope_key}",
                lambda cbg=cbg, indices=indices: setattr(cbg, "active", indices),
            )

    # 5. Filtering query builder
    saved_filtering = settings.get("filtering") or []
    if saved_filtering:
        new_sections = []
        for section_i, saved_section in enumerate(saved_filtering):
            section_data = create_or_section()
            section_data["rows"].clear()
            valid_saved_rows = []
            for saved_row in saved_section.get("rows", []):
                category = saved_row.get("category")
                column = saved_row.get("column")
                col_info = filtering_metadata.get(category, {}).get("columns", {}).get(column)
                if col_info is None:
                    _warn(f"filtering[{section_i}] row", f"column '{category}.{column}' not found in this database")
                    continue
                operator = saved_row.get("operator")
                row_data = create_query_row(
                    section_data, initial_category=category, initial_column=column, initial_operator=operator
                )
                section_data["rows"].append(row_data)
                valid_saved_rows.append(saved_row)
            if not section_data["rows"]:
                continue
            rebuild_section(section_data)
            for ri, (row_data, saved_row) in enumerate(zip(section_data["rows"], valid_saved_rows)):
                value = saved_row.get("value")
                widget = row_data["input_ref"]["widget"]

                def _set_row_value(widget=widget, value=value, row_data=row_data):
                    if isinstance(widget, SearchableSelect) and value and value not in widget.options:
                        widget.options = list(widget.options) + [value]
                    widget.value = value

                _restore(f"filtering[{section_i}] row value", _set_row_value)
                row_and_or = saved_row.get("row_and_or")
                if ri > 0 and row_and_or and row_data["and_div"] is not None:
                    _restore(
                        f"filtering[{section_i}] row_and_or",
                        lambda row_data=row_data, row_and_or=row_and_or: setattr(
                            row_data["and_div"], "value", row_and_or
                        ),
                    )
            new_sections.append((section_data, saved_section.get("section_and_or")))
        if new_sections:
            or_sections.clear()
            for section_data, _ in new_sections:
                or_sections.append(section_data)
            rebuild_filtering_content()
            for i, (_, section_and_or) in enumerate(new_sections):
                if i > 0 and section_and_or and (i - 1) < len(inter_section_selects):
                    _restore(
                        "filtering section_and_or",
                        lambda i=i, section_and_or=section_and_or: setattr(
                            inter_section_selects[i - 1], "value", section_and_or
                        ),
                    )
        else:
            _warn("filtering", "no valid rows could be restored — keeping the default empty filter")

    # 6. Coloring rules (custom + MAG track), mirroring on_template_change's restore pattern
    def _restore_color_rows(saved_rows, target_rows, rebuild_fn, label):
        target_rows.clear()
        for saved_rule in saved_rows or []:
            qualifier = saved_rule.get("qualifier")
            if qualifier not in color_qualifier_options:
                _warn(label, f"qualifier '{qualifier}' not found in this database")
                continue
            row_data = create_color_row(target_rows, rebuild_fn)
            row_data["qualifier_select"].value = qualifier
            operator = saved_rule.get("operator")
            if operator in row_data["operator_select"].options:
                row_data["operator_select"].value = operator
            if operator != "Use random colors":
                widget = row_data["input_ref"]["widget"]
                value = saved_rule.get("value")

                def _set_color_value(widget=widget, value=value):
                    if isinstance(widget, SearchableSelect) and value and value not in widget.options:
                        widget.options = list(widget.options) + [value]
                    widget.value = value

                _restore(f"{label} value", _set_color_value)
                color = saved_rule.get("color")
                if color:
                    _restore(
                        f"{label} color",
                        lambda row_data=row_data, color=color: setattr(row_data["color_picker"], "color", color),
                    )
            target_rows.append(row_data)
        rebuild_fn()

    coloring = contig_settings.get("coloring") or {}
    if coloring.get("custom_color_rows") is not None:
        _restore_color_rows(
            coloring.get("custom_color_rows"),
            custom_color_rows,
            rebuild_color_rows,
            "contig.coloring.custom_color_rows",
        )
    if coloring.get("mag_track_color_rows") is not None:
        _restore_color_rows(
            coloring.get("mag_track_color_rows"),
            mag_track_color_rows,
            rebuild_mag_track_color_rows,
            "contig.coloring.mag_track_color_rows",
        )
    if "apply_annotation_rules_to_mag_track" in coloring and apply_annotation_rules_cbg is not None:
        _restore(
            "contig.coloring.apply_annotation_rules_to_mag_track",
            lambda: setattr(
                apply_annotation_rules_cbg,
                "active",
                [0] if coloring["apply_annotation_rules_to_mag_track"] else [],
            ),
        )

    # 7. Genomic feature widgets
    _restore_genome_features(contig_settings, bindings)

    # 8. Plotting parameters
    pp = settings.get("plotting_params") or {}
    for key, widget in (
        ("min_coverage_freq", min_coverage_freq_input),
        ("max_genemap_window", max_genemap_window_input),
        ("max_sequence_window", max_sequence_window_input),
        ("max_binning_window", max_binning_window_input),
        ("genemap_height", genemap_height_input),
        ("sequence_height", sequence_height_input),
        ("translated_sequence_height", translated_sequence_height_input),
        ("subplot_height", subplot_height_input),
    ):
        if key in pp and pp[key] is not None:
            _restore(f"plotting_params.{key}", lambda widget=widget, value=pp[key]: setattr(widget, "value", value))

    mag_params = pp.get("mag_params") or {}
    if mag_params.get("category") in mag_params_category_select.options:
        _restore(
            "plotting_params.mag_params.category",
            lambda: setattr(mag_params_category_select, "value", mag_params["category"]),
        )
    elif mag_params.get("category"):
        _warn("plotting_params.mag_params.category", f"'{mag_params['category']}' not available in this database")
    _mag_metric_values = (
        [v for v, _l in mag_params_metric_select.options]
        if all(isinstance(o, tuple) for o in mag_params_metric_select.options)
        else list(mag_params_metric_select.options)
    )
    if mag_params.get("metric") in _mag_metric_values:
        _restore(
            "plotting_params.mag_params.metric",
            lambda: setattr(mag_params_metric_select, "value", mag_params["metric"]),
        )
    elif mag_params.get("metric"):
        _warn("plotting_params.mag_params.metric", f"'{mag_params['metric']}' not available for this category")
    if "direction" in mag_params:
        _restore(
            "plotting_params.mag_params.direction",
            lambda: setattr(mag_params_direction, "active", mag_params["direction"]),
        )
    if "sort_sample" in mag_params and (
        not mag_params["sort_sample"] or mag_params["sort_sample"] in mag_params_sort_sample_select.options
    ):
        _restore(
            "plotting_params.mag_params.sort_sample",
            lambda: setattr(mag_params_sort_sample_select, "value", mag_params["sort_sample"]),
        )
    if "max_dots" in mag_params and mag_params["max_dots"] is not None:
        _restore(
            "plotting_params.mag_params.max_dots",
            lambda: setattr(mag_track_max_dots_input, "value", mag_params["max_dots"]),
        )

    sample_params = pp.get("sample_params") or {}
    if "max_samples" in sample_params and sample_params["max_samples"] is not None:
        _restore(
            "plotting_params.sample_params.max_samples",
            lambda: setattr(max_samples_input, "value", sample_params["max_samples"]),
        )
    if sample_params.get("order_category") in sample_order_category_select.options:
        _restore(
            "plotting_params.sample_params.order_category",
            lambda: setattr(sample_order_category_select, "value", sample_params["order_category"]),
        )
    elif sample_params.get("order_category"):
        _warn(
            "plotting_params.sample_params.order_category",
            f"'{sample_params['order_category']}' not available in this database",
        )
    _metric_values = (
        [v for v, _l in sample_order_metric_select.options]
        if all(isinstance(o, tuple) for o in sample_order_metric_select.options)
        else sample_order_metric_select.options
    )
    if sample_params.get("order_metric") in _metric_values:
        _restore(
            "plotting_params.sample_params.order_metric",
            lambda: setattr(sample_order_metric_select, "value", sample_params["order_metric"]),
        )
    elif sample_params.get("order_metric"):
        _warn(
            "plotting_params.sample_params.order_metric",
            f"'{sample_params['order_metric']}' not available for this category",
        )
    if "order_direction" in sample_params:
        _restore(
            "plotting_params.sample_params.order_direction",
            lambda: setattr(sample_order_direction, "active", sample_params["order_direction"]),
        )
    if sample_params.get("same_y_scale") is not None:
        _restore(
            "plotting_params.sample_params.same_y_scale",
            lambda: setattr(same_y_scale_cbg, "active", [0] if sample_params["same_y_scale"] else []),
        )
