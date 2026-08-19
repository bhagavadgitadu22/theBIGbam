"""Stable session callbacks for controllers which are wired later in construction."""

from __future__ import annotations

from typing import Any, Callable


class SessionCallbacks:
    def __init__(self, schedule: Callable[[Callable[[], None]], None]) -> None:
        self.main_placeholder: Any | None = None
        self.schedule = schedule
        self._apply: Any | None = None
        self._summary: Any | None = None
        self._subject: Any | None = None
        self._scope: Any | None = None

    def attach_apply(self, controller: Any) -> None:
        self._apply = controller

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

    def apply_clicked(self, _event: Any = None) -> None:
        self._required(self.main_placeholder, "Main placeholder").loading = True
        self.schedule(self.do_apply)

    def do_apply(self) -> None:
        self._required(self._apply, "Apply").apply()

    def show_summary(self, _event: Any = None) -> None:
        self._required(self._summary, "Summary").show()

    def sync_selected_contig_position(self) -> None:
        self._required(self._subject, "Subject").sync_selected_contig_position()

    def variable_callback(self, group: Any, _unused: Any = None):
        return self._required(self._scope, "Scope transition").variable_callback(group)

    def view_changed(self, attr: str, old: Any, new: Any) -> None:
        self._required(self._scope, "Scope transition").view_changed(attr, old, new)
