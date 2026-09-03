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
    install_filters: Callable[[Mapping[str, Any], Mapping[str, Any]], None] | None = None


class HistoryRestoreCoordinator:
    """Restore controls, settle dependencies, and replay a validated action."""

    def __init__(self, bindings: HistoryRestoreBindings) -> None:
        self.bindings = bindings

    def restore(self, entry, on_complete=None, on_error=None) -> bool:
        bindings = self.bindings
        if not bindings.begin():
            return False
        regenerates_plot = entry.action == "apply_plot"
        if regenerates_plot:
            bindings.set_loading(True)

        finished = False
        owns_lock = True

        def release_lock() -> None:
            nonlocal owns_lock
            if owns_lock:
                bindings.end()
                owns_lock = False

        def acquire_lock() -> bool:
            nonlocal owns_lock
            if owns_lock:
                return True
            owns_lock = bindings.begin()
            return owns_lock

        def finish(error: Exception | None = None) -> None:
            nonlocal finished
            if finished:
                return
            finished = True
            if regenerates_plot:
                bindings.set_loading(False)
            release_lock()
            if error is not None:
                if on_error is not None:
                    on_error(error)
                else:
                    print(f"[history] Restore invariant failed: {error}", flush=True)
            elif on_complete is not None:
                on_complete()

        def schedule(callback: Callable[[], None]) -> bool:
            try:
                bindings.schedule(callback)
                return True
            except Exception as error:
                finish(error)
                return False

        previous: Mapping[str, Any] | None = None

        def rollback(original_error: Exception) -> None:
            try:
                if previous is None:
                    raise RuntimeError("no previous settings snapshot is available")
                bindings.install(previous)
            except Exception as rollback_error:
                finish(RuntimeError(f"{original_error}; rollback failed: {rollback_error}"))
                return

            # Installation changes dependent Bokeh controls. Give their
            # callbacks one document tick to settle before redrawing too.
            def redraw_previous() -> None:
                if not acquire_lock():
                    finish(RuntimeError(f"{original_error}; rollback failed: interaction lock unavailable"))
                    return
                try:
                    if bindings.apply_plot() is False:
                        raise RuntimeError("rollback plot could not be regenerated")
                except Exception as rollback_error:
                    finish(RuntimeError(f"{original_error}; rollback failed: {rollback_error}"))
                    return
                finish(original_error)

            release_lock()
            schedule(redraw_previous)

        def regenerate_plot() -> None:
            if not acquire_lock():
                finish(RuntimeError("restored plot could not acquire the interaction lock"))
                return
            try:
                if bindings.apply_plot() is False:
                    raise RuntimeError("previously successful plot snapshot could not be regenerated")
            except Exception as error:
                rollback(error)
                return
            finish()

        def install_settings() -> None:
            nonlocal previous
            try:
                previous = bindings.snapshot()
                if not regenerates_plot and bindings.install_filters is not None:
                    bindings.install_filters(entry.settings, previous)
                else:
                    bindings.install(entry.settings)
            except Exception as error:
                if previous is not None:
                    try:
                        bindings.install(previous)
                    except Exception as rollback_error:
                        error = RuntimeError(f"{error}; rollback failed: {rollback_error}")
                finish(error)
                return
            if regenerates_plot:
                # Release the controls transaction before yielding. Scope and
                # subject callbacks reject work while it is busy, so keeping
                # this lock across the settling tick makes the first render
                # observe an intermediate state. The spinner remains active.
                release_lock()
                schedule(regenerate_plot)
            else:
                finish()

        return schedule(install_settings)
