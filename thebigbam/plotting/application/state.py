"""Typed, database-free state transitions for the plotting application."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import IntEnum
from typing import Callable, Iterable


class SampleScope(IntEnum):
    ONE = 0
    ALL = 1


class SubjectScope(IntEnum):
    MAG = 0
    CONTIG = 1


@dataclass(frozen=True)
class ModeSelections:
    """Variable labels selected independently in each sample scope."""

    one: tuple[str, ...] = ()
    all: tuple[str, ...] = ()


@dataclass(frozen=True)
class PlotUiState:
    sample_scope: SampleScope = SampleScope.ONE
    subject_scope: SubjectScope = SubjectScope.CONTIG
    sample: str = ""
    contig: str = ""
    mag: str = ""
    filter_revision: int = 0
    variables: ModeSelections = field(default_factory=ModeSelections)


@dataclass(frozen=True)
class AvailabilitySnapshot:
    samples: tuple[str, ...] = ()
    contigs: tuple[str, ...] = ()
    mags: tuple[str, ...] = ()

    def reconcile(self, state: PlotUiState) -> PlotUiState:
        """Clear selections which are no longer present in this snapshot."""
        return replace(
            state,
            sample=state.sample if state.sample in self.samples else "",
            contig=state.contig if state.contig in self.contigs else "",
            mag=state.mag if state.mag in self.mags else "",
        )


@dataclass(frozen=True)
class PlotRequest:
    """Validated identity of a requested plot.

    Display-specific options remain in the renderer during milestone one; this
    type establishes the stable state boundary used by controllers and tests.
    """

    sample_scope: SampleScope
    subject_scope: SubjectScope
    sample: str | None
    contig: str | None
    mag: str | None
    variables: tuple[str, ...]
    x_start: int
    x_end: int | None


def switch_sample_scope(state: PlotUiState, scope: SampleScope) -> PlotUiState:
    """Return a new state without performing I/O or changing selections."""
    return replace(state, sample_scope=SampleScope(scope))


def enforce_all_sample_variable(selections: Iterable[str], newly_selected: str | None = None) -> tuple[str, ...]:
    """Keep the most recently selected variable in ALL mode."""
    values = tuple(dict.fromkeys(selections))
    if newly_selected in values:
        return (newly_selected,)
    return values[-1:] if values else ()


class PlotController:
    """Own sample-scope state and cached availability without database access."""

    def __init__(
        self,
        state: PlotUiState,
        snapshots: dict[SampleScope, AvailabilitySnapshot],
        project: Callable[[PlotUiState, AvailabilitySnapshot], None] | None = None,
    ) -> None:
        self.state = state
        self._snapshots = dict(snapshots)
        self._project = project
        self.transition_count = 0

    def set_snapshot(self, scope: SampleScope, snapshot: AvailabilitySnapshot) -> None:
        self._snapshots[SampleScope(scope)] = snapshot

    def snapshot(self, scope: SampleScope | None = None) -> AvailabilitySnapshot:
        key = self.state.sample_scope if scope is None else SampleScope(scope)
        return self._snapshots[key]

    def switch_sample_scope(self, scope: SampleScope) -> PlotUiState:
        scope = SampleScope(scope)
        self.state = self.snapshot(scope).reconcile(switch_sample_scope(self.state, scope))
        self.transition_count += 1
        if self._project is not None:
            self._project(self.state, self.snapshot(scope))
        return self.state
