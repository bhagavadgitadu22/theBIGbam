"""Typed mutable state shared by one plotting session."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..repositories.preload import PreloadedPlotData
from ..shared.diagnostics import PlotDiagnostics


@dataclass
class CurrentPlotState:
    contig: str | None = None
    sample: str | None = None
    is_all: bool | None = None
    shared_xrange: Any = None
    data_xstart: int | None = None
    data_xend: int | None = None
    range_callbacks: tuple[tuple[Any, str, Any], ...] = ()

    def __getitem__(self, key: str):
        return getattr(self, key)

    def __setitem__(self, key: str, value) -> None:
        setattr(self, key, value)

    def get(self, key: str, default=None):
        return getattr(self, key, default)


@dataclass
class PlotSessionContext:
    connection: Any
    preloaded: PreloadedPlotData
    widgets: dict[str, Any]
    diagnostics: PlotDiagnostics
    plot_state: CurrentPlotState = field(default_factory=CurrentPlotState)
