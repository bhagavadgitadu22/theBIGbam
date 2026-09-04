"""Per-session history of successful user apply actions."""

from __future__ import annotations

import datetime as dt
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
