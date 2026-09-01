"""Per-session history of successful user apply actions."""

from __future__ import annotations

import datetime as dt
import uuid
from dataclasses import asdict, dataclass
from typing import Any, Mapping

HISTORY_FORMAT = "thebigbam-session-history"
HISTORY_VERSION = 1
VALID_ACTIONS = {"apply_filters", "apply_plot"}


@dataclass(frozen=True)
class HistoryEntry:
    id: str
    sequence: int
    action: str
    created_at: str
    settings: dict[str, Any]


def _filter_predicates(settings: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    return [
        row
        for section in settings.get("filtering") or []
        for row in section.get("rows") or []
        if row.get("value") not in (None, "")
    ]


def _selected_variables(settings: Mapping[str, Any], *, all_samples: bool) -> list[str]:
    selection_key = "selected_all" if all_samples else "selected_one"
    selected = []
    for module in (settings.get("variables") or {}).values():
        selected.extend(str(label) for label in module.get(selection_key) or [])
    return selected


def _limited(values: list[str], *, limit: int = 5) -> str:
    visible = values[:limit]
    text = ", ".join(visible)
    remaining = len(values) - len(visible)
    return f"{text} (+{remaining} more)" if remaining else text


def describe_history_entry(entry: HistoryEntry) -> str:
    """Build a concise, display-independent description of a history event."""
    settings = entry.settings
    predicates = _filter_predicates(settings)
    if entry.action == "apply_filters":
        lines = []
        if not predicates:
            lines.append("No active filters")
        else:
            lines.append(f"{len(predicates)} condition{'s' if len(predicates) != 1 else ''}:")
            lines.extend(
                f"• {row.get('category', '?')} · {row.get('column', '?')} "
                f"{row.get('operator', '=')} {row.get('value')}"
                for row in predicates[:5]
            )
            if len(predicates) > 5:
                lines.append(f"• +{len(predicates) - 5} more")
    else:
        view_mode = settings.get("view_mode") or {}
        all_samples = view_mode.get("one_or_all_samples") == 1
        mag_view = view_mode.get("mag_or_contig") == 0
        selection = settings.get("selection") or {}
        lines = [f"Mode: {'All samples' if all_samples else 'One sample'}"]
        subject_key = "mag" if mag_view else "contig"
        subject = selection.get(subject_key)
        if subject:
            lines.append(f"{'MAG' if mag_view else 'Contig'}: {subject}")
        if not all_samples and selection.get("sample"):
            lines.append(f"Sample: {selection['sample']}")
        position = ((settings.get("contig") or {}).get("position_range") or {})
        if position.get("from") is not None and position.get("to") is not None:
            lines.append(f"Range: {position['from']}–{position['to']}")
        variables = _selected_variables(settings, all_samples=all_samples)
        if variables:
            lines.append(f"Variable{'s' if len(variables) != 1 else ''}: {_limited(variables)}")
        lines.append(
            f"Filters: {len(predicates)} active condition{'s' if len(predicates) != 1 else ''}"
            if predicates
            else "Filters: No active filters"
        )
    return "\n".join(lines)


class SessionHistory:
    def __init__(self, source_db: str, *, limit_per_action: int = 100) -> None:
        self.source_db = source_db
        self.limit_per_action = limit_per_action
        self._entries: list[HistoryEntry] = []
        self._sequence = 0

    @property
    def entries(self) -> tuple[HistoryEntry, ...]:
        return tuple(self._entries)

    def append(self, action: str, settings: Mapping[str, Any]) -> HistoryEntry:
        if action not in VALID_ACTIONS:
            raise ValueError(f"unsupported history action: {action}")
        self._sequence += 1
        entry = HistoryEntry(
            id=str(uuid.uuid4()),
            sequence=self._sequence,
            action=action,
            created_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            # Settings collectors hand off a newly constructed JSON document.
            # History takes ownership and avoids a second full-tree copy on
            # the latency-sensitive APPLY completion path.
            settings=dict(settings),
        )
        self._entries.append(entry)
        same_action = [item for item in self._entries if item.action == action]
        if len(same_action) > self.limit_per_action:
            expired = same_action[: len(same_action) - self.limit_per_action]
            expired_ids = {item.id for item in expired}
            self._entries = [item for item in self._entries if item.id not in expired_ids]
        return entry

    def for_action(self, action: str) -> tuple[HistoryEntry, ...]:
        return tuple(item for item in self._entries if item.action == action)

    def remove(self, entry_id: str) -> bool:
        """Remove one retained event by stable identity."""
        retained = [item for item in self._entries if item.id != entry_id]
        if len(retained) == len(self._entries):
            return False
        self._entries = retained
        return True

    def document(self) -> dict[str, Any]:
        return {
            "_meta": {
                "format": HISTORY_FORMAT,
                "version": HISTORY_VERSION,
                "source_db": self.source_db,
                "saved_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            },
            "entries": [asdict(item) for item in self._entries],
        }
