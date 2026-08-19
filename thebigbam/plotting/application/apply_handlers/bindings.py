"""Dependencies supplied to Apply orchestration and mode handlers."""

from dataclasses import dataclass
from typing import Any, Callable


@dataclass(frozen=True)
class ApplyBindings:
    _encoding_by_feature: Any
    _get_filtered_samples_for_contig: Any
    _get_filtered_samples_for_mag: Any
    _mag_sort_category_sources: Any
    _sample_sort_category_sources: Any
    _send_timing_ping: Any
    _sync_from_to_for_selected_contig: Any
    apply_annotation_rules_cbg: Any
    combined_features_cbg: Any
    command_hint_pane: Any
    conn: Any
    current_plot_state: Any
    data_cache: Any
    custom_color_rows: Any
    db_path: Any
    diagnostics: Any
    download_data_button: Any
    download_mag_metrics_button: Any
    download_metrics_button: Any
    enable_timing: Any
    feature_label_select: Any
    feature_type_multichoice: Any
    from_position_input: Any
    genemap_height_input: Any
    mag_params_category_select: Any
    mag_params_direction: Any
    mag_params_metric_select: Any
    mag_params_sort_sample_select: Any
    mag_track_color_rows: Any
    mag_track_max_dots_input: Any
    main_placeholder: Any
    max_binning_window_input: Any
    max_genemap_window_input: Any
    max_samples_input: Any
    max_sequence_window_input: Any
    min_coverage_freq_input: Any
    peruse_button: Any
    plot_isoforms_cbg: Any
    plot_lifecycle: Any
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
    timing: Any
    set_operation: Callable[[str], None]
