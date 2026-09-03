"""Per-session history of successful user apply actions."""

from __future__ import annotations

import datetime as dt
import json
import uuid
from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping

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


@dataclass(frozen=True)
class HistoryDescriptionLine:
    """One stable, human-readable element in a history snapshot or diff."""

    key: str
    text: str
    removed: bool = False


_SETTING_LABELS = {
    "view_mode": "View mode",
    "mag_or_contig": "MAG/contig view",
    "one_or_all_samples": "Sample mode",
    "selection": "Selection",
    "sample": "Sample",
    "contig": "Contig",
    "mag": "MAG",
    "position_range": "Position range",
    "from": "From",
    "to": "To",
    "filtering": "Filters",
    "plotting_params": "Plotting parameters",
    "feature_widgets": "Genome features",
    "coloring": "Coloring",
    "custom_color_rows": "Annotation coloring rules",
    "mag_track_color_rows": "MAG coloring rules",
}


def _display_value(value: Any) -> str:
    if isinstance(value, bool):
        return "Enabled" if value else "Disabled"
    if value is None:
        return "None"
    if isinstance(value, str):
        return value or "None"
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(", ", ": "))


def _path_label(path: tuple[str, ...]) -> str:
    parts = []
    for token in path:
        if token.startswith("#"):
            parts.append(token)
        else:
            parts.append(_SETTING_LABELS.get(token, token.replace("_", " ").capitalize()))
    return " · ".join(parts)


def _flatten_settings(value: Any, path: tuple[str, ...] = ()) -> list[HistoryDescriptionLine]:
    lines: list[HistoryDescriptionLine] = []
    if isinstance(value, Mapping):
        for key, child in value.items():
            if key == "_meta":
                continue
            lines.extend(_flatten_settings(child, (*path, str(key))))
        return lines
    if isinstance(value, list):
        if not value:
            return [HistoryDescriptionLine(".".join(path), f"{_path_label(path)}: None")]
        if all(not isinstance(item, (Mapping, list)) for item in value):
            shown = ", ".join(_display_value(item) for item in value)
            return [HistoryDescriptionLine(".".join(path), f"{_path_label(path)}: {shown}")]
        for index, child in enumerate(value, start=1):
            lines.extend(_flatten_settings(child, (*path, f"#{index}")))
        return lines
    return [
        HistoryDescriptionLine(
            ".".join(path),
            f"{_path_label(path)}: {_display_value(value)}",
        )
    ]


def _filter_description_lines(settings: Mapping[str, Any]) -> list[HistoryDescriptionLine]:
    lines: list[HistoryDescriptionLine] = []
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
            prefix = ""
            if row_index > 1:
                prefix = f"{row.get('row_and_or') or 'AND'} · "
            condition = (
                f"{prefix}{row.get('category', '?')} · {row.get('column', '?')} "
                f"{row.get('operator', '=')} {_display_value(value)}"
            )
            lines.append(
                HistoryDescriptionLine(
                    f"filtering.section.{section_index}.row.{row_index}",
                    f"Filter #{section_index}.{row_index}: {condition}",
                )
            )
    return lines or [HistoryDescriptionLine("filtering", "Filters: None")]


def history_description_lines(entry: HistoryEntry) -> tuple[HistoryDescriptionLine, ...]:
    """Describe a complete entry from the canonical saved-settings document."""
    source: Mapping[str, Any]
    if entry.action == "apply_filters":
        return tuple(_filter_description_lines(entry.settings))
    else:
        source = entry.settings
    return tuple(_flatten_settings(source))


def history_diff_lines(
    previous: HistoryEntry | None, current: HistoryEntry
) -> tuple[HistoryDescriptionLine, ...]:
    """Describe additions, replacements, and removals from the prior same-action entry."""
    before = {line.key: line.text for line in history_description_lines(previous)} if previous else {}
    after_lines = history_description_lines(current)
    after = {line.key: line.text for line in after_lines}
    result: list[HistoryDescriptionLine] = []
    for line in after_lines:
        old = before.get(line.key)
        if old == line.text:
            continue
        if old is None:
            result.append(line)
            continue
        label, _, old_value = old.partition(": ")
        _new_label, _, new_value = line.text.partition(": ")
        result.append(HistoryDescriptionLine(line.key, f"{label}: {old_value} → {new_value}"))
    for key, text in before.items():
        if key not in after:
            result.append(HistoryDescriptionLine(key, text, removed=True))
    if not result:
        result.append(HistoryDescriptionLine("unchanged", "No changes"))
    return tuple(result)


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

    def restore(self, entries: Iterable[HistoryEntry]) -> None:
        """Replace retained entries with validated, previously saved history."""
        restored = list(entries)
        if len({entry.id for entry in restored}) != len(restored):
            raise ValueError("session history contains duplicate entry IDs")
        sequences = [entry.sequence for entry in restored]
        if any(sequence <= 0 for sequence in sequences) or sequences != sorted(set(sequences)):
            raise ValueError("session history entry sequences must be unique, positive, and increasing")
        if any(entry.action not in VALID_ACTIONS for entry in restored):
            raise ValueError("session history contains an unsupported action")
        self._entries = restored
        self._sequence = max(sequences, default=0)

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


def entries_from_session_document(document: Mapping[str, Any]) -> tuple[HistoryEntry, ...]:
    """Validate and deserialize a SAVE SESSION document."""
    meta = document.get("_meta")
    if not isinstance(meta, Mapping) or meta.get("format") != HISTORY_FORMAT:
        raise ValueError("session JSON has an unsupported format")
    if meta.get("version") != HISTORY_VERSION:
        raise ValueError(f"unsupported session history version: {meta.get('version')!r}")
    raw_entries = document.get("entries")
    if not isinstance(raw_entries, list):
        raise ValueError("session JSON 'entries' must be a list")

    entries = []
    for index, raw in enumerate(raw_entries):
        if not isinstance(raw, Mapping):
            raise ValueError(f"session history entry {index} must be an object")
        entry_id = raw.get("id")
        sequence = raw.get("sequence")
        action = raw.get("action")
        created_at = raw.get("created_at")
        settings = raw.get("settings")
        if not isinstance(entry_id, str) or not entry_id:
            raise ValueError(f"session history entry {index} has an invalid ID")
        if not isinstance(sequence, int) or isinstance(sequence, bool):
            raise ValueError(f"session history entry {index} has an invalid sequence")
        if action not in VALID_ACTIONS:
            raise ValueError(f"session history entry {index} has an unsupported action")
        if not isinstance(created_at, str) or not created_at:
            raise ValueError(f"session history entry {index} has an invalid timestamp")
        if not isinstance(settings, dict):
            raise ValueError(f"session history entry {index} has invalid settings")
        entries.append(HistoryEntry(entry_id, sequence, action, created_at, dict(settings)))

    # Reuse the model's cross-entry validation so imported and live histories
    # follow the same identity and sequencing invariants.
    validation = SessionHistory(str(meta.get("source_db", "")))
    validation.restore(entries)
    return validation.entries
