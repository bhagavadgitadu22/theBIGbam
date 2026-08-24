"""Deferred availability-controller access and autocomplete widget projection."""

from __future__ import annotations

from typing import Any


class AvailabilityFacade:
    """Provide stable callbacks before the availability controller is wired."""

    def __init__(self) -> None:
        self._controller: Any | None = None

    def attach(self, controller: Any) -> None:
        self._controller = controller

    @property
    def controller(self) -> Any:
        if self._controller is None:
            raise RuntimeError("Availability controller has not been attached")
        return self._controller

    @staticmethod
    def update_widget(widget: Any, completions: list[str]) -> None:
        if widget.options == completions:
            return
        current_value = widget.value
        widget.options = completions
        if current_value and current_value in completions and not widget.value:
            widget.value = current_value
        elif widget.value and widget.value not in completions:
            widget.value = ""

    @classmethod
    def push_search(cls, widget: Any, completions: list[str]) -> None:
        cls.update_widget(widget, completions)
        widget.param.trigger("options")

    def compute_contigs(self, term: str = ""):
        return self.controller.compute_contigs(term)

    def refresh_contigs(self, term: str = "") -> None:
        self.controller.refresh_contigs(term)

    def compute_samples(self, term: str = ""):
        return self.controller.compute_samples(term)

    def sort_samples_for_mag(self, completions, selected_mag=None):
        return self.controller.sort_samples_for_mag(completions, selected_mag)

    def refresh_samples(self, term: str = "") -> None:
        self.controller.refresh_samples(term)

    def compute_mags(self, term: str = ""):
        return self.controller.compute_mags(term)

    def refresh_mags(self, term: str = "") -> None:
        self.controller.refresh_mags(term)

    def update_titles(self) -> None:
        self.controller.update_titles()

    def invalidate_titles(self) -> None:
        self.controller.invalidate_titles()
