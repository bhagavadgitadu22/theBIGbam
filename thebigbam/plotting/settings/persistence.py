"""Serialization and persistence helpers for plotting settings."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


def serialize_variable_selection(widgets: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Return variable selections keyed by stable module names."""
    result = {}
    for index, module_name in enumerate(widgets["module_names"]):
        one = widgets["variables_widgets_one"][index]
        all_samples = widgets["variables_widgets_all"][index]
        result[module_name] = {
            "module_enabled": bool(widgets["module_widgets_one"][index].active),
            "selected_one": [one.labels[item] for item in one.active],
            "selected_all": [all_samples.labels[item] for item in all_samples.active],
        }
    return result


def serialize_filter_sections(
    sections: Sequence[Mapping[str, Any]], inter_section_selects: Sequence[Any]
) -> list[dict[str, Any]]:
    """Serialize the dynamic filtering UI without retaining Bokeh models."""
    result = []
    for section_index, section in enumerate(sections):
        section_operator = (
            inter_section_selects[section_index - 1].value if 0 < section_index <= len(inter_section_selects) else None
        )
        rows = []
        for row_index, row in enumerate(section["rows"]):
            rows.append(
                {
                    "category": row["category_select"].value,
                    "column": row["subcategory_select"].value,
                    "operator": row["comparison_select"].value,
                    "value": row["input_ref"]["widget"].value,
                    "row_and_or": (row["and_div"].value if row_index > 0 and row["and_div"] is not None else None),
                }
            )
        result.append({"section_and_or": section_operator, "rows": rows})
    return result


def serialize_color_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Serialize annotation coloring rules."""
    return [
        {
            "qualifier": row["qualifier_select"].value,
            "operator": row["operator_select"].value,
            "value": row["input_ref"]["widget"].value,
            "color": row["color_picker"].color,
        }
        for row in rows
    ]


def save_settings_document(
    settings: Mapping[str, Any],
    db_path: str | Path,
    directory: str | Path | None = None,
    now: Callable[[], dt.datetime] = dt.datetime.now,
) -> Path:
    """Write a settings document and return its generated path."""
    output_directory = Path.cwd() if directory is None else Path(directory)
    timestamp = now().strftime("%Y%m%d_%H%M%S")
    output_path = output_directory / f"{Path(db_path).stem}_{timestamp}.json"
    with output_path.open("w", encoding="utf-8") as stream:
        json.dump(settings, stream, indent=2)
    return output_path


def load_settings_document(path: str | Path) -> dict[str, Any]:
    """Load and minimally validate a plotting settings document."""
    with Path(path).open(encoding="utf-8") as stream:
        settings = json.load(stream)
    if not isinstance(settings, dict):
        raise ValueError("settings JSON must contain an object at its root")
    return settings
