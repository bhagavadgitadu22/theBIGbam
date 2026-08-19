"""Projection of live plotting controls into a persistent settings document."""

from __future__ import annotations

import datetime
import os
from dataclasses import dataclass
from typing import Any

from .persistence import serialize_color_rows, serialize_filter_sections, serialize_variable_selection


@dataclass(frozen=True)
class SettingsCollectorBindings:
    apply_annotation_rules_cbg: Any
    combined_features_cbg: Any
    custom_color_rows: Any
    db_path: Any
    feature_label_select: Any
    feature_type_multichoice: Any
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


class SettingsCollector:
    def __init__(self, bindings: SettingsCollectorBindings) -> None:
        self.bindings = bindings

    def collect(self) -> dict[str, Any]:
        bindings = self.bindings
        apply_annotation_rules_cbg = bindings.apply_annotation_rules_cbg
        combined_features_cbg = bindings.combined_features_cbg
        custom_color_rows = bindings.custom_color_rows
        db_path = bindings.db_path
        feature_label_select = bindings.feature_label_select
        feature_type_multichoice = bindings.feature_type_multichoice
        from_position_input = bindings.from_position_input
        genemap_height_input = bindings.genemap_height_input
        genome_master_cbg = bindings.genome_master_cbg
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
        plot_isoforms_cbg = bindings.plot_isoforms_cbg
        same_y_scale_cbg = bindings.same_y_scale_cbg
        sample_order_category_select = bindings.sample_order_category_select
        sample_order_direction = bindings.sample_order_direction
        sample_order_metric_select = bindings.sample_order_metric_select
        sequence_cbg = bindings.sequence_cbg
        sequence_height_input = bindings.sequence_height_input
        subplot_height_input = bindings.subplot_height_input
        to_position_input = bindings.to_position_input
        translated_sequence_cbg = bindings.translated_sequence_cbg
        translated_sequence_height_input = bindings.translated_sequence_height_input
        views = bindings.views
        widgets = bindings.widgets

        settings = {
            "_meta": {
                "source_db": os.path.basename(db_path),
                "saved_at": datetime.datetime.now().isoformat(timespec="seconds"),
            },
            "view_mode": {
                "mag_or_contig": widgets["view_radio"].active,
                "one_or_all_samples": views.active,
            },
            "selection": {
                "sample": widgets["sample_select"].value,
                "contig": widgets["contig_select"].value,
                "mag": widgets["mag_select"].value,
            },
            "contig": {
                "position_range": {
                    "from": from_position_input.value,
                    "to": to_position_input.value,
                },
                "feature_widgets": {
                    "feature_types": list(feature_type_multichoice.value)
                    if feature_type_multichoice is not None
                    else None,
                    "plot_isoforms": bool(plot_isoforms_cbg.active) if plot_isoforms_cbg is not None else None,
                    "combined_features": [combined_features_cbg.labels[i] for i in combined_features_cbg.active]
                    if combined_features_cbg is not None
                    else None,
                    "sequence": bool(sequence_cbg.active) if sequence_cbg is not None else None,
                    "translated_sequence": bool(translated_sequence_cbg.active)
                    if translated_sequence_cbg is not None
                    else None,
                    "genome_master": bool(genome_master_cbg.active) if genome_master_cbg is not None else None,
                    "feature_label": feature_label_select.value if feature_label_select is not None else None,
                },
                "coloring": {
                    "custom_color_rows": serialize_color_rows(custom_color_rows),
                    "mag_track_color_rows": serialize_color_rows(mag_track_color_rows),
                    "apply_annotation_rules_to_mag_track": bool(apply_annotation_rules_cbg.active)
                    if apply_annotation_rules_cbg is not None
                    else None,
                },
            },
            "variables": serialize_variable_selection(widgets),
            "filtering": serialize_filter_sections(or_sections, inter_section_selects),
            "plotting_params": {
                "min_coverage_freq": min_coverage_freq_input.value,
                "max_genemap_window": max_genemap_window_input.value,
                "max_sequence_window": max_sequence_window_input.value,
                "max_binning_window": max_binning_window_input.value,
                "genemap_height": genemap_height_input.value,
                "sequence_height": sequence_height_input.value,
                "translated_sequence_height": translated_sequence_height_input.value,
                "subplot_height": subplot_height_input.value,
                "mag_params": {
                    "category": mag_params_category_select.value,
                    "metric": mag_params_metric_select.value,
                    "direction": mag_params_direction.active,
                    "sort_sample": mag_params_sort_sample_select.value,
                    "max_dots": mag_track_max_dots_input.value,
                },
                "sample_params": {
                    "max_samples": max_samples_input.value,
                    "order_category": sample_order_category_select.value,
                    "order_metric": sample_order_metric_select.value,
                    "order_direction": sample_order_direction.active,
                    "same_y_scale": bool(same_y_scale_cbg.active),
                },
            },
        }
        return settings
