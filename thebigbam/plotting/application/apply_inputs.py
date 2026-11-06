"""Parse live plotting controls into an immutable Apply request snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..models.filters import FilterExpression
from ..models.read_events import READ_DISCREPANCY_CHOICES
from ..models.read_plot import ReadPlotOptions
from ..shared.defaults import DEFAULT_READ_DISCREPANCY_WINDOW


@dataclass(frozen=True)
class ApplyInputs:
    contig: str
    has_samples: bool
    is_all: bool
    sample: str | None
    selected_feature_types: Any
    genbank_path: str | None
    plot_isoforms: bool
    feature_label_key: str | None
    custom_colors: list[dict[str, Any]]
    mag_track_colors: list[dict[str, Any]]
    max_track_dots: int
    active_variables_widgets: Any
    max_genemap_window: int
    same_y_scale: bool
    show_subplot_titles: bool
    freeze_gene_map: bool
    subplot_size: int
    genemap_size: int
    sequence_size: int
    translated_sequence_size: int
    max_binning: int
    min_coverage_freq: float
    read_plot: ReadPlotOptions = ReadPlotOptions()


_OPERATOR_MODES = {
    "=": "exact",
    "!=": "not_equal",
    "has": "has",
    "has not": "has_not",
    "<": "lt",
    ">": "gt",
    "Use random colors": "random",
}


def _color_rules(rows, *, allow_random: bool) -> list[dict[str, Any]]:
    rules = []
    for row in rows:
        key = row["qualifier_select"].value
        if not key:
            continue
        mode = _OPERATOR_MODES.get(row["operator_select"].value, "exact")
        if mode == "random":
            if allow_random:
                rules.append({"qualifier_key": key, "match_mode": "random"})
            continue
        raw_value = row["input_ref"]["widget"].value
        color = row["color_picker"].color
        if raw_value is None or raw_value == "" or not color:
            continue
        value = float(raw_value) if isinstance(raw_value, (int, float)) else raw_value
        rules.append({"qualifier_key": key, "value": value, "color": color, "match_mode": mode})
    return rules


def collect_apply_inputs(bindings) -> ApplyInputs:
    """Read all controls once so downstream branches operate on stable values."""
    widgets = bindings.widgets
    has_samples = widgets["has_samples"]
    is_all = bindings.views.active == 1 if has_samples or widgets.get("bam_only_samples") else False
    feature_types = bindings.feature_type_multichoice.value if bindings.feature_type_multichoice is not None else None
    genbank_path = bindings.db_path if feature_types else None
    custom_colors = _color_rules(bindings.custom_color_rows, allow_random=True)
    mag_colors = _color_rules(bindings.mag_track_color_rows, allow_random=False)
    active_variables = widgets["variables_widgets_all"] if is_all else widgets["variables_widgets_one"]
    read_controls = widgets.get("read_plot_controls")
    read_options = ReadPlotOptions()
    if read_controls is not None:
        discrepancy_window = getattr(read_controls, "max_discrepancy_window", None)
        read_options = ReadPlotOptions(
            enabled=0 in (read_controls.all_reads if is_all else read_controls.add_reads).active,
            max_reads=int(read_controls.max_reads.value or 0),
            link_pairs=0 in read_controls.link_pairs.active,
            plot_height=int(read_controls.plot_height.value),
            read_height=int(read_controls.read_height.value),
            event_variables=frozenset(
                event
                for choice in read_controls.discrepancies.value
                for event in READ_DISCREPANCY_CHOICES[choice]
            ),
            max_discrepancy_window=int(discrepancy_window.value) if discrepancy_window is not None and discrepancy_window.value is not None else DEFAULT_READ_DISCREPANCY_WINDOW,
            filters=(read_controls.filters.expression() if getattr(read_controls, "filters", None) is not None else FilterExpression(())),
            color_rules=read_controls.color_rules.rules(),
        )
    if (
        bindings.apply_annotation_rules_cbg is not None
        and 0 in bindings.apply_annotation_rules_cbg.active
        and custom_colors
    ):
        mag_colors.extend(rule for rule in custom_colors if rule.get("match_mode") != "random")
    return ApplyInputs(
        contig=widgets["contig_select"].value,
        has_samples=has_samples,
        is_all=is_all,
        sample=(widgets["sample_select"].value if has_samples or widgets.get("bam_only_samples") else None),
        selected_feature_types=feature_types,
        genbank_path=genbank_path,
        plot_isoforms=(
            0 in bindings.plot_isoforms_cbg.active if bindings.plot_isoforms_cbg is not None and genbank_path else True
        ),
        feature_label_key=(bindings.feature_label_select.value if bindings.feature_label_select is not None else None),
        custom_colors=custom_colors,
        mag_track_colors=mag_colors,
        max_track_dots=int(bindings.mag_track_max_dots_input.value),
        active_variables_widgets=active_variables,
        max_genemap_window=int(bindings.max_genemap_window_input.value),
        same_y_scale=0 in bindings.same_y_scale_cbg.active,
        show_subplot_titles=0 in bindings.show_subplot_titles_cbg.active,
        freeze_gene_map=0 in bindings.freeze_gene_map_cbg.active,
        subplot_size=int(bindings.subplot_height_input.value),
        genemap_size=int(bindings.genemap_height_input.value),
        sequence_size=int(bindings.sequence_height_input.value),
        translated_sequence_size=int(bindings.translated_sequence_height_input.value),
        max_binning=int(bindings.max_binning_window_input.value),
        min_coverage_freq=float(bindings.min_coverage_freq_input.value),
        read_plot=read_options,
    )
