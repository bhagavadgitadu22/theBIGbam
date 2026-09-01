"""Crash-safe recording of interactive plotting scenarios."""

from __future__ import annotations

import copy
import json
import os
import queue
import threading
import uuid
from pathlib import Path
from typing import Any, Mapping


class ScenarioFormatError(ValueError):
    """Raised when a scenario document cannot be described safely."""


def _json_value(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _changed_values(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    """Flatten settings-shaped changes into deterministic ontology paths."""
    if isinstance(value, dict):
        flattened: list[tuple[str, Any]] = []
        for key in sorted(value):
            path = f"{prefix}.{key}" if prefix else str(key)
            flattened.extend(_changed_values(value[key], path))
        return flattened
    return [(prefix, value)]


def _action_label(step: Mapping[str, Any]) -> str:
    action = step["action"]
    if action == "filter_lookup":
        details = step.get("details", {})
        category = details.get("category")
        column = details.get("column")
        target = ".".join(str(value) for value in (category, column) if value not in (None, ""))
        return f"Look up filter values for {target}" if target else "Look up filter values"
    if action == "apply_filters":
        return "Apply filters"
    if action == "apply_plot":
        return "Apply plot"
    if action == "restore_history":
        details = step.get("details", {})
        history_action = str(details.get("history_action", "history")).replace("_", " ")
        sequence = details.get("history_sequence")
        return f"Restore {history_action} history entry {sequence}"
    if action == "remove_history":
        details = step.get("details", {})
        history_action = str(details.get("history_action", "history")).replace("_", " ")
        sequence = details.get("history_sequence")
        return f"Remove {history_action} history entry {sequence}"
    labels = {
        "show_summary": "Show summary",
        "download_contig_metrics": "Download contig metrics",
        "download_mag_metrics": "Download MAG metrics",
        "show_download_command": "Show download command",
        "save_settings": "Save settings",
        "save_session": "Save session",
        "reset_position": "Reset genomic position",
        "filter_distribution_scale": "Change filter distribution scale",
    }
    if action in labels:
        return labels[action]
    if action == "state_change":
        return "Change settings"
    return str(action).replace("_", " ").capitalize()


def describe_scenario_document(document: Any) -> tuple[str, ...]:
    """Return one deterministic human-readable line per recorded scenario step."""
    if not isinstance(document, dict):
        raise ScenarioFormatError("scenario root must be a JSON object")
    metadata = document.get("_meta")
    if not isinstance(metadata, dict) or metadata.get("format") != "thebigbam-scenario":
        raise ScenarioFormatError("not a thebigbam-scenario document")
    steps = document.get("steps")
    if not isinstance(steps, list):
        raise ScenarioFormatError("'steps' must be a JSON array")

    lines = []
    seen_sequences = set()
    for index, step in enumerate(steps, start=1):
        if not isinstance(step, dict):
            raise ScenarioFormatError(f"step {index} must be a JSON object")
        sequence = step.get("sequence")
        if not isinstance(sequence, int) or isinstance(sequence, bool) or sequence < 1:
            raise ScenarioFormatError(f"step {index} has an invalid sequence")
        if sequence in seen_sequences:
            raise ScenarioFormatError(f"duplicate step sequence {sequence}")
        seen_sequences.add(sequence)
        action = step.get("action")
        if not isinstance(action, str) or not action:
            raise ScenarioFormatError(f"step {sequence} has an invalid action")

        description = _action_label(step)
        changes = step.get("changes")
        if changes is not None:
            if not isinstance(changes, dict):
                raise ScenarioFormatError(f"step {sequence} has invalid changes")
            rendered = [f"{path}={_json_value(value)}" for path, value in _changed_values(changes)]
            if rendered:
                description += ": " + "; ".join(rendered)

        annotations = []
        status = step.get("status")
        if status is not None:
            if not isinstance(status, str) or not status:
                raise ScenarioFormatError(f"step {sequence} has an invalid status")
            annotations.append(status)
        duration = step.get("duration_seconds")
        if duration is not None:
            if not isinstance(duration, (int, float)) or isinstance(duration, bool) or duration < 0:
                raise ScenarioFormatError(f"step {sequence} has an invalid duration_seconds")
            annotations.append(f"{duration:g} s")
        memory = step.get("memory")
        if memory is not None:
            if not isinstance(memory, dict):
                raise ScenarioFormatError(f"step {sequence} has invalid memory")
            for key, label in (
                ("server_rss_mb", "server"),
                ("browser_heap_mb", "browser"),
            ):
                value = memory.get(key)
                if value is None:
                    continue
                if not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0:
                    raise ScenarioFormatError(f"step {sequence} has invalid memory.{key}")
                annotations.append(f"{label} {value:g} MB")
        if annotations:
            description += f" [{', '.join(annotations)}]"
        lines.append(f"{sequence}. {description}")
    return tuple(lines)


def describe_scenario_file(path: str | Path) -> tuple[str, ...]:
    """Load and describe a scenario, translating input errors for concise CLI output."""
    scenario_path = Path(path)
    try:
        with scenario_path.open(encoding="utf-8") as stream:
            document = json.load(stream)
    except OSError as error:
        raise ScenarioFormatError(str(error)) from error
    except json.JSONDecodeError as error:
        raise ScenarioFormatError(f"invalid JSON at line {error.lineno}, column {error.colno}") from error
    return describe_scenario_document(document)


def _settings_state(settings: Mapping[str, Any]) -> dict[str, Any]:
    """Copy settings while discarding SAVE SETTINGS-specific volatile metadata."""
    state = copy.deepcopy(dict(settings))
    state.pop("_meta", None)
    return state


def _diff(previous: Any, current: Any) -> Any:
    """Return a JSON merge-like recursive diff, or ``None`` when unchanged."""
    if isinstance(previous, dict) and isinstance(current, dict):
        changed = {}
        for key in previous.keys() | current.keys():
            if key not in current:
                changed[key] = None
            elif key not in previous:
                changed[key] = copy.deepcopy(current[key])
            else:
                value = _diff(previous[key], current[key])
                if value is not None:
                    changed[key] = value
        return changed or None
    if previous == current:
        return None
    return copy.deepcopy(current)


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(document, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


class ScenarioRecorder:
    """Record settings-shaped state transitions without blocking Bokeh callbacks."""

    def __init__(self, path: str | Path, db_path: str | Path, initial_settings: Mapping[str, Any]) -> None:
        self.path = Path(path)
        initial_state = _settings_state(initial_settings)
        self._document: dict[str, Any] = {
            "_meta": {
                "format": "thebigbam-scenario",
                "version": 1,
                "source_db": Path(db_path).name,
            },
            "initial_state": initial_state,
            "steps": [],
            "final_state": copy.deepcopy(initial_state),
        }
        self._state = initial_state
        self._lock = threading.Lock()
        self._writes: queue.Queue[dict[str, Any] | None] = queue.Queue()
        self._closed = False
        self._error: Exception | None = None
        self._worker = threading.Thread(target=self._write_loop, name="thebigbam-scenario-writer", daemon=True)
        self._worker.start()
        self._publish()

    @property
    def error(self) -> Exception | None:
        return self._error

    def _write_loop(self) -> None:
        while True:
            document = self._writes.get()
            try:
                if document is None:
                    return
                if self._error is None:
                    _atomic_write_json(self.path, document)
            except Exception as error:  # recording must never break the UI
                self._error = error
                print(f"[scenario] Recording disabled after write failure: {error}", flush=True)
            finally:
                self._writes.task_done()

    def _publish(self) -> None:
        self._writes.put(copy.deepcopy(self._document))

    def record_state(self, settings: Mapping[str, Any], *, action: str = "state_change") -> bool:
        current = _settings_state(settings)
        with self._lock:
            if self._closed:
                return False
            changes = _diff(self._state, current)
            if changes is None:
                return False
            self._document["steps"].append(
                {
                    "sequence": len(self._document["steps"]) + 1,
                    "action": action,
                    "changes": changes,
                }
            )
            self._state = current
            self._document["final_state"] = copy.deepcopy(current)
            self._publish()
        return True

    def record_action(
        self,
        action: str,
        settings: Mapping[str, Any] | None = None,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> bool:
        """Record a semantic user action, including any state not observed yet."""
        current = self._state if settings is None else _settings_state(settings)
        with self._lock:
            if self._closed:
                return False
            step = {
                "sequence": len(self._document["steps"]) + 1,
                "action": action,
            }
            if details is not None:
                step["details"] = json.loads(json.dumps(details))
            changes = _diff(self._state, current)
            if changes is not None:
                step["changes"] = changes
                self._state = copy.deepcopy(current)
                self._document["final_state"] = copy.deepcopy(current)
            self._document["steps"].append(step)
            self._publish()
        return True

    def flush(self) -> None:
        self._writes.join()

    def close(self, settings: Mapping[str, Any] | None = None) -> None:
        if settings is not None:
            self.record_state(settings)
        with self._lock:
            if self._closed:
                return
            self._closed = True
            self._publish()
        self.flush()
        self._writes.put(None)
        self._writes.join()
        self._worker.join(timeout=1)


class ScenarioPathAllocator:
    """Give concurrent Panel sessions distinct output files."""

    def __init__(self, requested_path: str | Path) -> None:
        self.path = Path(requested_path)
        self._count = 0
        self._lock = threading.Lock()

    def next_path(self) -> Path:
        with self._lock:
            self._count += 1
            if self._count == 1:
                return self.path
            return self.path.with_name(f"{self.path.stem}.session-{self._count}{self.path.suffix or '.json'}")
