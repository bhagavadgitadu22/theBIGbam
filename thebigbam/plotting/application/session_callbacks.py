"""Stable session callbacks for controllers which are wired later in construction."""

from __future__ import annotations

from typing import Any, Callable


class SessionCallbacks:
    def __init__(self, schedule: Callable[[Callable[[], None]], None], interactions: Any | None = None) -> None:
        self.main_placeholder: Any | None = None
        self.schedule = schedule
        self._apply: Any | None = None
        self._summary: Any | None = None
        self._subject: Any | None = None
        self._scope: Any | None = None
        self.interactions = interactions
        self._scenario_recorder: Any | None = None
        self._collect_settings: Callable[[], Any] | None = None
        self._collect_applied_settings: Callable[[], Any] | None = None
        self._plot_succeeded: Callable[[Any, int | None], None] | None = None
        self._pending_history_settings: Any | None = None
        self._pending_history_apply_step: int | None = None

    def attach_scenario(self, recorder: Any, collect_settings: Callable[[], Any]) -> None:
        self._scenario_recorder = recorder
        self._collect_settings = collect_settings

    def attach_apply(self, controller: Any) -> None:
        self._apply = controller

    def attach_history(
        self,
        collect_applied_settings: Callable[[], Any],
        plot_succeeded: Callable[[Any, int | None], None],
    ) -> None:
        self._collect_applied_settings = collect_applied_settings
        self._plot_succeeded = plot_succeeded

    def attach_summary(self, controller: Any) -> None:
        self._summary = controller

    def attach_subject(self, controller: Any) -> None:
        self._subject = controller

    def attach_scope(self, controller: Any) -> None:
        self._scope = controller

    def attach_placeholder(self, placeholder: Any) -> None:
        self.main_placeholder = placeholder

    @staticmethod
    def _required(controller: Any | None, name: str) -> Any:
        if controller is None:
            raise RuntimeError(f"{name} controller has not been attached")
        return controller

    def apply_clicked(
        self,
        _event: Any = None,
    ) -> bool:
        if self.interactions is not None and not self.interactions.begin("plot"):
            return False
        placeholder = self._required(self.main_placeholder, "Main placeholder")
        placeholder.loading = True
        try:
            if self._scenario_recorder is not None:
                self._pending_history_apply_step = self._scenario_recorder.record_action(
                    "apply_plot", self._collect_settings()
                )
            self._pending_history_settings = (
                self._collect_applied_settings()
                if self._collect_applied_settings is not None
                else None
            )
            self.schedule(self.do_apply)
        except Exception:
            self._pending_history_settings = None
            self._pending_history_apply_step = None
            placeholder.loading = False
            if self.interactions is not None:
                self.interactions.end()
            raise
        return True

    def do_apply(self) -> None:
        try:
            succeeded = self._required(self._apply, "Apply").apply()
            if succeeded is not False and self._pending_history_settings is not None and self._plot_succeeded is not None:
                self._plot_succeeded(
                    self._pending_history_settings, self._pending_history_apply_step
                )
        finally:
            self._pending_history_settings = None
            self._pending_history_apply_step = None
            self._required(self.main_placeholder, "Main placeholder").loading = False
            if self.interactions is not None:
                self.interactions.end()

    def apply_restored_now(self) -> bool:
        """Apply an already validated restored state inside its owning transaction."""
        return self._required(self._apply, "Apply").apply()

    def set_plot_loading(self, loading: bool) -> None:
        self._required(self.main_placeholder, "Main placeholder").loading = loading

    def show_summary(self, _event: Any = None) -> None:
        self._required(self._summary, "Summary").show()

    def sync_selected_contig_position(self) -> None:
        self._required(self._subject, "Subject").sync_selected_contig_position()

    def variable_callback(self, group: Any, _unused: Any = None):
        return self._required(self._scope, "Scope transition").variable_callback(group)

    def view_changed(self, attr: str, old: Any, new: Any) -> None:
        self._required(self._scope, "Scope transition").view_changed(attr, old, new)
