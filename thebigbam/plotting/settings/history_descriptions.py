"""Human-facing descriptions derived from canonical plotting settings."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

from ..shared.defaults import (
    DEFAULT_GENEMAP_HEIGHT,
    DEFAULT_GENEMAP_WINDOW,
    DEFAULT_MAG_MAX_DOTS,
    DEFAULT_MAG_ORDER_DIRECTION,
    DEFAULT_MAX_BASE_RESOLUTION,
    DEFAULT_MAX_SAMPLES,
    DEFAULT_MIN_COVERAGE_FREQUENCY,
    DEFAULT_SAMPLE_ORDER_DIRECTION,
    DEFAULT_SEQUENCE_HEIGHT,
    DEFAULT_SEQUENCE_WINDOW,
    DEFAULT_SUBPLOT_HEIGHT,
    DEFAULT_TRANSLATED_SEQUENCE_HEIGHT,
)


@dataclass(frozen=True)
class HistoryDescriptionLine:
    """One stable, human-readable element in a history snapshot or diff."""

    key: str
    text: str
    removed: bool = False
    default: bool = False


@dataclass(frozen=True)
class HistoryDescriptionContext:
    """Database-dependent defaults used only for concise presentation."""

    mag_order: tuple[Any, Any] = (None, None)
    sample_order: tuple[Any, Any] = ("Sample", "Sample_name")


_SETTING_LABELS = {
    "view_mode": "View mode",
    "mag_or_contig": "MAG/contig view",
    "one_or_all_samples": "Sample mode",
    "selection": "Selection",
    "sample": "Sample",
    "contig": "Contig",
    "mag": "MAG",
    "position_range": "Position range",
    "plotting_params": "Plotting parameters",
    "mag_params": "MAG params",
    "sample_params": "Sample params",
    "feature_widgets": "Genome features",
    "coloring": "Coloring",
    "apply_annotation_rules_to_mag_track": "Apply annotation rules to MAG track",
}
_ENUM_LABELS = {
    ("view_mode", "mag_or_contig"): {0: "MAG", 1: "Contig"},
    ("view_mode", "one_or_all_samples"): {0: "One sample", 1: "All samples"},
}
_SELECTION_ORDER = {"mag": 0, "contig": 1, "sample": 2}
_DEFAULT_VALUES = {
    ("plotting_params", "min_coverage_freq"): DEFAULT_MIN_COVERAGE_FREQUENCY,
    ("plotting_params", "max_genemap_window"): DEFAULT_GENEMAP_WINDOW,
    ("plotting_params", "max_sequence_window"): DEFAULT_SEQUENCE_WINDOW,
    ("plotting_params", "max_binning_window"): DEFAULT_MAX_BASE_RESOLUTION,
    ("plotting_params", "genemap_height"): DEFAULT_GENEMAP_HEIGHT,
    ("plotting_params", "sequence_height"): DEFAULT_SEQUENCE_HEIGHT,
    ("plotting_params", "translated_sequence_height"): DEFAULT_TRANSLATED_SEQUENCE_HEIGHT,
    ("plotting_params", "subplot_height"): DEFAULT_SUBPLOT_HEIGHT,
    ("plotting_params", "mag_params", "max_dots"): DEFAULT_MAG_MAX_DOTS,
    ("plotting_params", "sample_params", "max_samples"): DEFAULT_MAX_SAMPLES,
    ("plotting_params", "sample_params", "same_y_scale"): False,
}
_MISSING = object()


def _display_value(value: Any) -> str:
    if isinstance(value, bool):
        return "Enabled" if value else "Disabled"
    if value is None or value == "":
        return "None"
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(", ", ": "))


def _display_setting_value(path: tuple[str, ...], value: Any) -> str:
    return _ENUM_LABELS.get(path, {}).get(value, _display_value(value))


def _path_label(path: tuple[str, ...]) -> str:
    return " · ".join(
        token if token.startswith("#") else _SETTING_LABELS.get(token, token.replace("_", " ").capitalize())
        for token in path
    )


def _is_described_value(value: Any) -> bool:
    return value is not None and value is not False and value != ""


def _position_range(value: Mapping[str, Any], path: tuple[str, ...]) -> list[HistoryDescriptionLine]:
    start, end = value.get("from"), value.get("to")
    available = [item for item in (start, end) if _is_described_value(item)]
    if not available:
        return []
    displayed = "–".join(_display_value(item) for item in available)
    return [HistoryDescriptionLine(".".join(path), f"{_path_label(path)}: {displayed}")]


def _color_rules(value: Sequence[Any], path: tuple[str, ...]) -> list[HistoryDescriptionLine]:
    prefix = "Contig · Coloring ·" if path[-1] == "custom_color_rows" else "Contig · MAG coloring ·"
    lines = []
    for index, rule in enumerate(value, start=1):
        if not isinstance(rule, Mapping) or not all(
            _is_described_value(rule.get(field)) for field in ("qualifier", "value", "color")
        ):
            continue
        lines.append(
            HistoryDescriptionLine(
                f"{'.'.join(path)}.#{index}",
                f"{prefix} #{index}: {rule['qualifier']} {rule.get('operator') or '='} "
                f"{rule['value']} · {rule['color']}",
            )
        )
    return lines


_MAPPING_FORMATTERS: dict[
    tuple[str, ...], Callable[[Mapping[str, Any], tuple[str, ...]], list[HistoryDescriptionLine]]
] = {("contig", "position_range"): _position_range}
_LIST_FORMATTERS = {
    "custom_color_rows": _color_rules,
    "mag_track_color_rows": _color_rules,
}


def _variable_lines(
    value: Mapping[str, Any], *, all_samples: bool
) -> list[HistoryDescriptionLine]:
    lines = []
    selection_key = "selected_all" if all_samples else "selected_one"
    for module, state in value.items():
        if not isinstance(state, Mapping):
            continue
        selected = [str(item) for item in state.get(selection_key) or []]
        if not selected:
            continue
        displayed = "All" if not all_samples and state.get("module_enabled") is True else ", ".join(selected)
        lines.append(
            HistoryDescriptionLine(
                f"variables.{module}.{selection_key}",
                f"Variables · {module}: {displayed}",
            )
        )
    return lines


def _ordering_group(
    value: Mapping[str, Any],
    path: tuple[str, ...],
    context: HistoryDescriptionContext,
    *,
    all_samples: bool,
) -> list[HistoryDescriptionLine]:
    is_sample = path[-1] == "sample_params"
    category_key = "order_category" if is_sample else "category"
    metric_key = "order_metric" if is_sample else "metric"
    category, metric = value.get(category_key), value.get(metric_key)
    default_pair = context.sample_order if is_sample else context.mag_order
    direction = value.get("order_direction" if is_sample else "direction")
    default_direction = DEFAULT_SAMPLE_ORDER_DIRECTION if is_sample else DEFAULT_MAG_ORDER_DIRECTION
    label = "Sample order" if is_sample else "MAG order"
    lines = []
    if _is_described_value(category) or _is_described_value(metric):
        pair = (category, metric)
        displayed = " · ".join(str(item) for item in pair if _is_described_value(item))
        if direction in (0, 1):
            displayed += f" · {'↑' if direction == 0 else '↓'}"
        lines.append(
            HistoryDescriptionLine(
                f"{'.'.join(path)}.order",
                f"Plotting parameters · {label}: {displayed}",
                default=pair == default_pair and direction == default_direction,
            )
        )
    for key, child in value.items():
        if key == "sort_sample" and not all_samples:
            continue
        if key not in {category_key, metric_key, "order_direction", "direction"}:
            lines.extend(_flatten_settings(child, (*path, str(key)), context=context))
    return lines


def _flatten_settings(
    value: Any,
    path: tuple[str, ...] = (),
    *,
    context: HistoryDescriptionContext,
    all_samples: bool = False,
    mag_view: bool = False,
) -> list[HistoryDescriptionLine]:
    if isinstance(value, Mapping):
        if path == ("variables",):
            return _variable_lines(value, all_samples=all_samples)
        if path == ("plotting_params", "sample_params"):
            return _ordering_group(value, path, context, all_samples=True) if all_samples else []
        if path == ("plotting_params", "mag_params"):
            return (
                _ordering_group(value, path, context, all_samples=all_samples)
                if mag_view
                else []
            )
        formatter = _MAPPING_FORMATTERS.get(path)
        if formatter is not None:
            return formatter(value, path)
        items = value.items()
        if path == ("selection",):
            items = sorted(items, key=lambda item: _SELECTION_ORDER.get(str(item[0]), 3))
        lines = []
        for key, child in items:
            if key == "_meta":
                continue
            if not path and key == "filtering":
                filter_lines = _filter_description_lines({"filtering": child})
                if not (len(filter_lines) == 1 and filter_lines[0].key == "filtering"):
                    lines.extend(filter_lines)
                continue
            lines.extend(
                _flatten_settings(
                    child,
                    (*path, str(key)),
                    context=context,
                    all_samples=all_samples,
                    mag_view=mag_view,
                )
            )
        return lines
    if isinstance(value, list):
        if not value:
            return []
        formatter = _LIST_FORMATTERS.get(path[-1]) if path else None
        if formatter is not None:
            return formatter(value, path)
        if all(not isinstance(item, (Mapping, list)) for item in value):
            shown = ", ".join(_display_value(item) for item in value if _is_described_value(item))
            return (
                [HistoryDescriptionLine(".".join(path), f"{_path_label(path)}: {shown}")]
                if shown
                else []
            )
        lines = []
        for index, child in enumerate(value, start=1):
            lines.extend(
                _flatten_settings(
                    child,
                    (*path, f"#{index}"),
                    context=context,
                    all_samples=all_samples,
                    mag_view=mag_view,
                )
            )
        return lines
    default_value = _DEFAULT_VALUES.get(path, _MISSING)
    if default_value is _MISSING and not _is_described_value(value):
        return []
    return [
        HistoryDescriptionLine(
            ".".join(path),
            f"{_path_label(path)}: {_display_setting_value(path, value)}",
            default=default_value is not _MISSING and value == default_value,
        )
    ]


def _filter_description_lines(settings: Mapping[str, Any]) -> list[HistoryDescriptionLine]:
    lines = []
    for section_index, section in enumerate(settings.get("filtering") or [], start=1):
        if section_index > 1:
            connector = section.get("section_and_or") or "AND"
            lines.append(
                HistoryDescriptionLine(
                    f"filtering.section.{section_index}.connector",
                    f"Filter group #{section_index}: {connector}",
                )
            )
        for row_index, row in enumerate(section.get("rows") or [], start=1):
            value = row.get("value")
            if value in (None, ""):
                continue
            connector = f"{row.get('row_and_or') or 'AND'} · " if row_index > 1 else ""
            lines.append(
                HistoryDescriptionLine(
                    f"filtering.section.{section_index}.row.{row_index}",
                    f"Filter #{section_index}.{row_index}: {connector}{row.get('category', '?')} · "
                    f"{row.get('column', '?')} {row.get('operator', '=')} {_display_value(value)}",
                )
            )
    return lines or [HistoryDescriptionLine("filtering", "Filters: None")]


def canonical_history_description_lines(
    entry: Any, context: HistoryDescriptionContext | None = None
) -> tuple[HistoryDescriptionLine, ...]:
    """Describe an entry while retaining defaults needed for later diffs."""
    context = context or HistoryDescriptionContext()
    view_mode = entry.settings.get("view_mode") or {}
    lines = (
        _filter_description_lines(entry.settings)
        if entry.action == "apply_filters"
        else _flatten_settings(
            entry.settings,
            context=context,
            all_samples=view_mode.get("one_or_all_samples") == 1,
            mag_view=view_mode.get("mag_or_contig") == 0,
        )
    )
    return tuple(lines)


def history_description_lines(
    entry: Any, context: HistoryDescriptionContext | None = None
) -> tuple[HistoryDescriptionLine, ...]:
    """Describe a complete entry without routine default-valued controls."""
    return tuple(line for line in canonical_history_description_lines(entry, context) if not line.default)


def diff_description_lines(
    previous: Sequence[HistoryDescriptionLine], current: Sequence[HistoryDescriptionLine]
) -> tuple[HistoryDescriptionLine, ...]:
    """Compare preformatted entries without traversing either settings tree again."""
    before = {line.key: line for line in previous}
    after = {line.key: line for line in current}
    result = []
    for line in current:
        old = before.get(line.key)
        if old is None:
            if not line.default:
                result.append(line)
        elif old.text != line.text:
            result.append(line)
    result.extend(
        HistoryDescriptionLine(key, line.text, removed=True)
        for key, line in before.items()
        if key not in after
    )
    return tuple(result or [HistoryDescriptionLine("unchanged", "No changes")])
