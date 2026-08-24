"""Scoped busy-state coordination for sidebar transitions."""

from __future__ import annotations

from typing import Any


class InteractionCoordinator:
    """Disable a registered control tree while one serialized transition runs."""

    def __init__(self) -> None:
        self.root: Any | None = None
        self.scope: str | None = None
        self._disabled_before: list[tuple[Any, bool]] = []
        self._root_classes_before: list[str] = []

    @property
    def busy(self) -> bool:
        return self.scope is not None

    def attach(self, root: Any) -> None:
        self.root = root

    def _walk(self, node: Any, seen: set[int]):
        if node is None or id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attribute in ("objects", "children"):
            for child in getattr(node, attribute, ()) or ():
                yield from self._walk(child, seen)
        wrapped = getattr(node, "object", None)
        if wrapped is not None and wrapped is not node:
            yield from self._walk(wrapped, seen)

    def begin(self, scope: str) -> bool:
        if self.busy or self.root is None:
            return False
        self.scope = scope
        self._disabled_before = []
        self._root_classes_before = list(getattr(self.root, "css_classes", ()))
        if hasattr(self.root, "css_classes"):
            self.root.css_classes = [*self._root_classes_before, "sidebar-busy"]
        for control in self._walk(self.root, set()):
            if not hasattr(control, "disabled"):
                continue
            previous = bool(control.disabled)
            self._disabled_before.append((control, previous))
            control.disabled = True
        return True

    def end(self) -> None:
        for control, previous in reversed(self._disabled_before):
            control.disabled = previous
        self._disabled_before = []
        if self.root is not None and hasattr(self.root, "css_classes"):
            self.root.css_classes = list(self._root_classes_before)
        self._root_classes_before = []
        self.scope = None
