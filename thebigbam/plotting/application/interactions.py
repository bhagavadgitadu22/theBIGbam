"""Scoped busy-state coordination for sidebar transitions."""

from __future__ import annotations

from typing import Any


class InteractionCoordinator:
    """Block a sidebar with one root-level update while a transition runs."""

    def __init__(self) -> None:
        self.root: Any | None = None
        self.scope: str | None = None
        self._root_classes_before: list[str] = []

    @property
    def busy(self) -> bool:
        return self.scope is not None

    def attach(self, root: Any) -> None:
        self.root = root

    def begin(self, scope: str) -> bool:
        if self.busy or self.root is None:
            return False
        self.scope = scope
        self._root_classes_before = list(getattr(self.root, "css_classes", ()))
        if hasattr(self.root, "css_classes"):
            self.root.css_classes = [*self._root_classes_before, "sidebar-busy"]
        return True

    def end(self) -> None:
        if self.root is not None and hasattr(self.root, "css_classes"):
            self.root.css_classes = list(self._root_classes_before)
        self._root_classes_before = []
        self.scope = None
