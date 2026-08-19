"""Lightweight structured diagnostics for plotting sessions."""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from typing import TextIO


@dataclass(frozen=True)
class DiagnosticEvent:
    operation: str
    elapsed_seconds: float
    generation: int = 0
    models: int = 0
    data_sources: int = 0
    payload_bytes: int = 0
    query_count: int = 0
    callbacks: int = 0
    property_updates: int = 0


class PlotDiagnostics:
    def __init__(self, enabled: bool = False, stream: TextIO | None = None, max_events: int = 1_000) -> None:
        self.enabled = enabled
        self.stream = stream
        self.max_events = max(1, max_events)
        self.events: list[DiagnosticEvent] = []
        self.generation = 0

    def next_generation(self) -> int:
        self.generation += 1
        return self.generation

    def record(self, operation: str, started_at: float, **metrics: int) -> DiagnosticEvent:
        event = DiagnosticEvent(
            operation=operation,
            elapsed_seconds=time.perf_counter() - started_at,
            generation=self.generation,
            **metrics,
        )
        if self.enabled:
            self.events.append(event)
            if len(self.events) > self.max_events:
                del self.events[: len(self.events) - self.max_events]
            line = json.dumps(asdict(event), sort_keys=True)
            print(f"[diagnostic] {line}", file=self.stream, flush=True)
        return event

    def snapshot(self) -> tuple[DiagnosticEvent, ...]:
        return tuple(self.events)
