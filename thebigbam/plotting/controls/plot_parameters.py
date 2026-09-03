"""Construction of plotting-parameter controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import panel as pn
from bokeh.layouts import row
from bokeh.models import Div
from bokeh.models.widgets import Button, CheckboxGroup, RadioButtonGroup, Select, Spinner

from ..shared.defaults import (
    AUTOCOMPLETE_LIMIT,
    DEFAULT_GENEMAP_WINDOW,
    DEFAULT_MAX_BASE_RESOLUTION,
    DEFAULT_SEQUENCE_WINDOW,
)
from .parameter_options import ParameterOptionCatalog
from .searchable_select import SearchableSelect, decode_search_request


@dataclass(frozen=True)
class PlotParameterControls:
    separator: Any
    header: Any
    content: Any
    min_coverage_freq: Any
    max_genemap_window: Any
    max_sequence_window: Any
    max_binning_window: Any
    genemap_height: Any
    sequence_height: Any
    translated_sequence_height: Any
    subplot_height: Any
    mag_header: Any
    mag_content: Any
    mag_category: Any
    mag_metric: Any
    mag_direction: Any
    mag_sort_sample: Any
    mag_sort_sample_row: Any
    mag_max_dots: Any
    mag_category_sources: Mapping[str, str]
    sample_header: Any
    sample_content: Any
    max_samples: Any
    sample_order_category: Any
    sample_order_metric: Any
    sample_order_direction: Any
    same_y_scale: Any
    sample_category_sources: Mapping[str, str]
    sample_contig_categories: list[str]
    sample_mag_categories: list[str]
    sample_current_categories: list[str]


def _collapsible_section(title, children, toggle_stylesheet, make_toggle_callback):
    toggle = Button(
        label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet]
    )
    toggle.styles = {"padding": "0px", "line-height": "20px"}
    header = row(
        toggle,
        Div(text=title, align="center"),
        sizing_mode="stretch_width",
        align="center",
        margin=(5, 0, 0, 0),
    )
    content = pn.Column(*children, sizing_mode="stretch_width", visible=False)
    toggle.on_click(make_toggle_callback(toggle, content))
    return header, content


def build_plot_parameter_controls(
    widgets: Mapping[str, Any],
    sample_scope: Any,
    original_samples: list[str],
    parameter_options: ParameterOptionCatalog,
    toggle_stylesheet: Any,
    stylesheet: Any,
    make_toggle_callback: Callable[[Any, Any], Callable[..., None]],
    compute_sample_completions: Callable[[str], list[str]],
    push_search_completions: Callable[[Any, list[str]], None],
    sort_sample_completions: Callable[[list[str]], list[str]],
    interaction_lock: Mapping[str, bool],
) -> PlotParameterControls:
    ## Plotting parameters section
    separator_plotting_params = Div(
        text="",
        height=2,
        sizing_mode="stretch_width",
        styles={"background-color": "#333", "margin-top": "10px", "margin-bottom": "10px"},
    )
    plotting_params_title = Div(
        text="<span style='font-size: 1.2em;'><b>Plotting parameters</b></span>", align="center"
    )
    plotting_params_header = row(plotting_params_title, sizing_mode="stretch_width", align="center")

    ## Plotting parameters useful in both views
    min_coverage_freq_input = Spinner(value=0.0, low=0.0, high=1.0, step=0.01, width=100, margin=(0, 2, 0, 0))
    min_coverage_freq_label = Div(text="Minimum frequency for coverage-related features", margin=(5, 0, 5, 5))
    min_coverage_freq_row = row(
        min_coverage_freq_input, min_coverage_freq_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0)
    )

    # Subsection: Max window sizes for plotting
    max_genemap_window_input = Spinner(
        value=DEFAULT_GENEMAP_WINDOW,
        low=10,
        high=1_000_000_000,
        step=1000,
        width=100,
        margin=(0, 2, 0, 0),
    )
    max_genemap_window_label = Div(text="Gene map (bp)", margin=(5, 0, 5, 5))
    max_genemap_window_row = row(
        max_genemap_window_input, max_genemap_window_label, sizing_mode="stretch_width", margin=(5, 0, 5, 0)
    )

    max_sequence_window_input = Spinner(
        value=DEFAULT_SEQUENCE_WINDOW, low=10, high=1000000, step=100, width=100, margin=(0, 2, 0, 0)
    )
    max_sequence_window_label = Div(text="Sequence plots (bp)", margin=(5, 0, 5, 5))
    max_sequence_window_row = row(
        max_sequence_window_input, max_sequence_window_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0)
    )

    max_binning_window_input = Spinner(
        value=DEFAULT_MAX_BASE_RESOLUTION, low=10, high=1000000, step=1000, width=100, margin=(0, 2, 0, 0)
    )
    max_binning_window_label = Div(text="Feature plots without binning (bp)", margin=(5, 0, 5, 5))
    max_binning_window_row = row(
        max_binning_window_input, max_binning_window_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0)
    )

    max_window_header, max_window_content = _collapsible_section(
        "Max window size for plotting",
        [max_genemap_window_row, max_sequence_window_row, max_binning_window_row],
        toggle_stylesheet,
        make_toggle_callback,
    )

    # Subsection: Plot heights
    genemap_height_input = Spinner(value=100, low=10, high=1000, step=10, width=100, margin=(0, 2, 0, 0))
    genemap_height_label = Div(text="Of gene map (px)", margin=(5, 0, 5, 5))
    genemap_height_row = row(
        genemap_height_input, genemap_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0)
    )

    sequence_height_input = Spinner(value=50, low=10, high=1000, step=10, width=100, margin=(0, 2, 0, 0))
    sequence_height_label = Div(text="Of nucleotide sequence (px)", margin=(5, 0, 5, 5))
    sequence_height_row = row(
        sequence_height_input, sequence_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0)
    )

    translated_sequence_height_input = Spinner(value=50, low=10, high=1000, step=10, width=100, margin=(0, 2, 0, 0))
    translated_sequence_height_label = Div(text="Of translated sequence (px)", margin=(5, 0, 5, 5))
    translated_sequence_height_row = row(
        translated_sequence_height_input,
        translated_sequence_height_label,
        sizing_mode="stretch_width",
        margin=(0, 0, 5, 0),
    )

    subplot_height_input = Spinner(value=100, low=10, high=1000, step=10, width=100, margin=(0, 2, 0, 0))
    subplot_height_label = Div(text="Per feature plot (px)", margin=(5, 0, 5, 5))
    subplot_height_row = row(subplot_height_input, subplot_height_label, sizing_mode="stretch_width")

    plot_heights_header, plot_heights_content = _collapsible_section(
        "Plot heights",
        [genemap_height_row, sequence_height_row, translated_sequence_height_row, subplot_height_row],
        toggle_stylesheet,
        make_toggle_callback,
    )

    # MAG parameters (contig ordering)
    _mag_sort_categories = parameter_options.mag_categories
    _mag_sort_category_sources = parameter_options.mag_sources
    _metrics_by_category = parameter_options.metrics
    _sample_metrics_by_category = parameter_options.sample_metrics
    _mag_metrics_by_category = parameter_options.mag_metrics
    _initial_metrics = _mag_metrics_by_category.get("Contig", []) if "Contig" in _mag_sort_categories else []

    mag_params_toggle_btn = Button(
        label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet]
    )
    mag_params_toggle_btn.styles = {"padding": "0px", "line-height": "20px"}
    mag_params_title = Div(text="MAG parameters", align="center")
    mag_params_header = row(
        mag_params_toggle_btn, mag_params_title, sizing_mode="stretch_width", align="center", margin=(5, 0, 0, 0)
    )
    mag_params_header.visible = bool(widgets["has_mags"])

    mag_params_order_label = Div(text="Order contigs by:", margin=(5, 0, 2, 0))
    mag_params_category_select = Select(
        name="benchmark-mag-order-category",
        value=_mag_sort_categories[0] if _mag_sort_categories else "",
        options=_mag_sort_categories,
        sizing_mode="stretch_width",
        margin=(0, 2, 0, 0),
    )
    _initial_metric_values = [v for v, _l in _initial_metrics]
    mag_params_metric_select = Select(
        name="benchmark-mag-order-metric",
        value="Contig_length"
        if "Contig_length" in _initial_metric_values
        else (_initial_metric_values[0] if _initial_metric_values else ""),
        options=_initial_metrics if _initial_metrics else [""],
        sizing_mode="stretch_width",
        margin=(0, 2, 0, 2),
    )
    mag_params_direction = RadioButtonGroup(labels=["↑", "↓"], active=1, width=60, margin=(0, 0, 0, 2))
    mag_params_controls_row = row(
        mag_params_category_select,
        mag_params_metric_select,
        mag_params_direction,
        sizing_mode="stretch_width",
        margin=(0, 10, 0, 0),
    )

    mag_params_sort_sample_label = Div(text="Using values from sample", margin=(5, 5, 5, 5))
    mag_params_sort_sample_select = SearchableSelect(
        value=original_samples[0] if original_samples else "",
        options=list(original_samples[:AUTOCOMPLETE_LIMIT]),
        placeholder="Type to search samples...",
        server_search=True,
        sizing_mode="stretch_width",
        margin=(0, 0, 0, 0),
    )
    # A Panel component can't be a child of bokeh's native row() (it expects
    # Bokeh model children), so this uses pn.Row instead — same pattern as
    # sample_section above, which mixes bokeh Divs with a SearchableSelect.
    mag_params_sort_sample_row = pn.Row(
        mag_params_sort_sample_label, mag_params_sort_sample_select, sizing_mode="stretch_width", margin=(5, 10, 0, 0)
    )
    mag_params_sort_sample_row.visible = False

    def _on_sort_sample_search(event):
        if interaction_lock["locked"]:
            return
        request_nonce, query = decode_search_request(event.new)
        base = compute_sample_completions(query)
        push_search_completions(mag_params_sort_sample_select, sort_sample_completions(base))
        mag_params_sort_sample_select.search_result_query = query
        mag_params_sort_sample_select.search_result_nonce = request_nonce

    mag_params_sort_sample_select.param.watch(_on_sort_sample_search, "search_request")

    mag_track_max_dots_input = Spinner(value=1000, low=1, step=100, width=100, margin=(0, 2, 0, 0))
    mag_track_max_dots_label = Div(text="Maximum number of points on MAG track", margin=(5, 0, 5, 5))
    mag_track_max_dots_row = row(
        mag_track_max_dots_input, mag_track_max_dots_label, sizing_mode="stretch_width", margin=(5, 10, 0, 0)
    )

    mag_params_content = pn.Column(
        mag_params_order_label,
        mag_params_controls_row,
        mag_params_sort_sample_row,
        mag_track_max_dots_row,
        sizing_mode="stretch_width",
        visible=False,
    )
    mag_params_toggle_btn.on_click(make_toggle_callback(mag_params_toggle_btn, mag_params_content))

    def _on_mag_sort_category_change(attr, old, new):
        # Cache lookup — O(1), no re-derivation from filtering_metadata.
        metrics = _mag_metrics_by_category.get(new, [])
        new_options = metrics if metrics else [""]
        if mag_params_metric_select.options != new_options:
            mag_params_metric_select.options = new_options
        # Clamp the selected value to the first valid option for the new category.
        valid_values = {v for v, _l in metrics} if metrics else set()
        if mag_params_metric_select.value not in valid_values:
            mag_params_metric_select.value = metrics[0][0] if metrics else ""
        is_sample_dep = new != "Contig"
        is_all = sample_scope.active == 1
        mag_params_sort_sample_row.visible = is_sample_dep and is_all

    mag_params_category_select.on_change("value", _on_mag_sort_category_change)

    # Sample parameters (only useful in All Samples view)
    sample_params_toggle_btn = Button(
        label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet]
    )
    sample_params_toggle_btn.styles = {"padding": "0px", "line-height": "20px"}
    sample_params_title = Div(text="Sample parameters", align="center")
    sample_params_header = row(
        sample_params_toggle_btn, sample_params_title, sizing_mode="stretch_width", align="center", margin=(5, 0, 0, 0)
    )
    sample_params_header.visible = False  # Only shown in All Samples mode

    max_samples_input = Spinner(value=20, low=1, high=500, step=5, width=100, margin=(0, 2, 0, 0))
    max_samples_input.name = "benchmark-max-samples"
    max_samples_label = Div(text="Max number of samples plotted", margin=(5, 0, 5, 5))
    max_samples_row = row(max_samples_input, max_samples_label, sizing_mode="stretch_width", margin=(5, 0, 0, 0))

    # Sample ordering: two-row layout matching "Order contigs by:"
    _sample_contig_categories = parameter_options.sample_contig_categories
    _sample_mag_categories = parameter_options.sample_mag_categories
    _sample_sort_category_sources = parameter_options.sample_sources
    _sample_sort_current_categories = list(_sample_contig_categories)

    def _sample_order_columns_for(category: str) -> list[tuple[str, str]]:
        # Cache lookup — reads from _sample_metrics_by_category built at startup.
        return _sample_metrics_by_category.get(category, [])

    _initial_sample_metrics = _sample_order_columns_for("Sample")
    _initial_sample_metric_values = [value for value, _label in _initial_sample_metrics]
    sample_order_label = Div(text="Order samples by:", margin=(5, 0, 2, 0))
    sample_order_category_select = Select(
        name="benchmark-sample-order-category",
        value=_sample_sort_current_categories[0] if _sample_sort_current_categories else "",
        options=_sample_sort_current_categories,
        sizing_mode="stretch_width",
        margin=(0, 2, 0, 0),
    )
    sample_order_metric_select = Select(
        name="benchmark-sample-order-metric",
        value="Sample_name"
        if "Sample_name" in _initial_sample_metric_values
        else (_initial_sample_metric_values[0] if _initial_sample_metric_values else ""),
        options=_initial_sample_metrics if _initial_sample_metrics else [""],
        sizing_mode="stretch_width",
        margin=(0, 2, 0, 2),
    )
    sample_order_direction = RadioButtonGroup(labels=["↑", "↓"], active=0, width=60, margin=(0, 0, 0, 2))
    sample_order_controls_row = row(
        sample_order_category_select,
        sample_order_metric_select,
        sample_order_direction,
        sizing_mode="stretch_width",
        margin=(0, 10, 5, 0),
    )

    def _on_sample_order_category_change(attr, old, new):
        # Cache lookup — O(1), no re-derivation from filtering_metadata.
        metrics = _sample_order_columns_for(new)
        new_options = metrics if metrics else [""]
        if sample_order_metric_select.options != new_options:
            sample_order_metric_select.options = new_options
        # Clamp the selected value to the first valid option for the new category.
        valid_values = {v for v, _l in metrics} if metrics else set()
        if sample_order_metric_select.value not in valid_values:
            sample_order_metric_select.value = metrics[0][0] if metrics else ""

    sample_order_category_select.on_change("value", _on_sample_order_category_change)

    same_y_scale_cbg = CheckboxGroup(
        name="benchmark-same-y-scale",
        labels=["Use same y scale for all samples"],
        active=[],
        margin=(5, 0, 5, 0),
    )
    same_y_scale_row = row(same_y_scale_cbg, sizing_mode="stretch_width")

    sample_params_content = pn.Column(
        max_samples_row,
        sample_order_label,
        sample_order_controls_row,
        same_y_scale_row,
        sizing_mode="stretch_width",
        visible=False,
    )
    sample_params_toggle_btn.on_click(make_toggle_callback(sample_params_toggle_btn, sample_params_content))

    plotting_params_content = pn.Column(
        min_coverage_freq_row,
        max_window_header,
        max_window_content,
        plot_heights_header,
        plot_heights_content,
        mag_params_header,
        mag_params_content,
        sample_params_header,
        sample_params_content,
        sizing_mode="stretch_width",
    )

    # Hide samples-only parameters when no BAM files were provided
    if not widgets["has_samples"]:
        min_coverage_freq_row.visible = False
        # sample_params_header is already visible=False by default (All Samples mode only)

    return PlotParameterControls(
        separator=separator_plotting_params,
        header=plotting_params_header,
        content=plotting_params_content,
        min_coverage_freq=min_coverage_freq_input,
        max_genemap_window=max_genemap_window_input,
        max_sequence_window=max_sequence_window_input,
        max_binning_window=max_binning_window_input,
        genemap_height=genemap_height_input,
        sequence_height=sequence_height_input,
        translated_sequence_height=translated_sequence_height_input,
        subplot_height=subplot_height_input,
        mag_header=mag_params_header,
        mag_content=mag_params_content,
        mag_category=mag_params_category_select,
        mag_metric=mag_params_metric_select,
        mag_direction=mag_params_direction,
        mag_sort_sample=mag_params_sort_sample_select,
        mag_sort_sample_row=mag_params_sort_sample_row,
        mag_max_dots=mag_track_max_dots_input,
        mag_category_sources=_mag_sort_category_sources,
        sample_header=sample_params_header,
        sample_content=sample_params_content,
        max_samples=max_samples_input,
        sample_order_category=sample_order_category_select,
        sample_order_metric=sample_order_metric_select,
        sample_order_direction=sample_order_direction,
        same_y_scale=same_y_scale_cbg,
        sample_category_sources=_sample_sort_category_sources,
        sample_contig_categories=_sample_contig_categories,
        sample_mag_categories=_sample_mag_categories,
        sample_current_categories=_sample_sort_current_categories,
    )
