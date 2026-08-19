"""Parse live plotting controls into an immutable Apply request snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


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
    subplot_size: int
    genemap_size: int
    sequence_size: int
    translated_sequence_size: int
    max_binning: int
    min_coverage_freq: float


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
    is_all = bindings.views.active == 1 if has_samples else False
    feature_types = bindings.feature_type_multichoice.value if bindings.feature_type_multichoice is not None else None
    genbank_path = bindings.db_path if feature_types else None
    custom_colors = _color_rules(bindings.custom_color_rows, allow_random=True)
    mag_colors = _color_rules(bindings.mag_track_color_rows, allow_random=False)
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
        sample=widgets["sample_select"].value if has_samples else None,
        selected_feature_types=feature_types,
        genbank_path=genbank_path,
        plot_isoforms=(
            0 in bindings.plot_isoforms_cbg.active if bindings.plot_isoforms_cbg is not None and genbank_path else True
        ),
        feature_label_key=(bindings.feature_label_select.value if bindings.feature_label_select is not None else None),
        custom_colors=custom_colors,
        mag_track_colors=mag_colors,
        max_track_dots=int(bindings.mag_track_max_dots_input.value),
        active_variables_widgets=(widgets["variables_widgets_all"] if is_all else widgets["variables_widgets_one"]),
        max_genemap_window=int(bindings.max_genemap_window_input.value),
        same_y_scale=0 in bindings.same_y_scale_cbg.active,
        subplot_size=int(bindings.subplot_height_input.value),
        genemap_size=int(bindings.genemap_height_input.value),
        sequence_size=int(bindings.sequence_height_input.value),
        translated_sequence_size=int(bindings.translated_sequence_height_input.value),
        max_binning=int(bindings.max_binning_window_input.value),
        min_coverage_freq=float(bindings.min_coverage_freq_input.value),
    )
