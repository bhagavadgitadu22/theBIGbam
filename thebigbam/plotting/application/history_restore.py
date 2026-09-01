"""Atomic restoration of a previously successful history snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping


@dataclass(frozen=True)
class HistoryRestoreBindings:
    begin: Callable[[], bool]
    end: Callable[[], None]
    set_loading: Callable[[bool], None]
    schedule: Callable[[Callable[[], None]], None]
    snapshot: Callable[[], Mapping[str, Any]]
    install: Callable[[Mapping[str, Any]], None]
    apply_plot: Callable[[], bool]


class HistoryRestoreCoordinator:
    """Run one previously validated action as a single interaction."""

    def __init__(self, bindings: HistoryRestoreBindings) -> None:
        self.bindings = bindings

    def restore(self, entry, on_complete=None, on_error=None) -> bool:
        bindings = self.bindings
        if not bindings.begin():
            return False
        bindings.set_loading(True)

        def run() -> None:
            error: Exception | None = None
            previous: Mapping[str, Any] | None = None
            try:
                previous = bindings.snapshot()
                bindings.install(entry.settings)
                if entry.action == "apply_plot" and bindings.apply_plot() is False:
                    raise RuntimeError("previously successful plot snapshot could not be regenerated")
            except Exception as caught:
                error = caught
                if previous is not None:
                    try:
                        bindings.install(previous)
                        if entry.action == "apply_plot" and bindings.apply_plot() is False:
                            raise RuntimeError("rollback plot could not be regenerated")
                    except Exception as rollback_error:
                        error = RuntimeError(f"{caught}; rollback failed: {rollback_error}")
            finally:
                bindings.set_loading(False)
                bindings.end()
            if error is not None:
                if on_error is not None:
                    on_error(error)
                else:
                    print(f"[history] Restore invariant failed: {error}", flush=True)
                return
            if on_complete is not None:
                on_complete()

        try:
            bindings.schedule(run)
        except Exception:
            bindings.set_loading(False)
            bindings.end()
            raise
        return True
