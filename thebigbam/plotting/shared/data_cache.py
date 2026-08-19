"""Bounded session-local cache for plain plotting data."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Hashable


@dataclass(frozen=True)
class CacheStats:
    hits: int
    misses: int
    entries: int
    invalidations: int
    last_invalidation: str | None


class SessionDataCache:
    """Store only plain decoded/service data; never Bokeh models."""

    def __init__(self, max_entries: int = 256) -> None:
        self.max_entries = max_entries
        self._values: OrderedDict[Hashable, Any] = OrderedDict()
        self.hits = 0
        self.misses = 0
        self.invalidations = 0
        self.last_invalidation: str | None = None

    def get(self, key: Hashable) -> tuple[bool, Any]:
        if key not in self._values:
            self.misses += 1
            return False, None
        self.hits += 1
        self._values.move_to_end(key)
        return True, self._values[key]

    def put(self, key: Hashable, value: Any) -> None:
        self._values[key] = value
        self._values.move_to_end(key)
        while len(self._values) > self.max_entries:
            self._values.popitem(last=False)

    def invalidate(self, reason: str) -> None:
        self._values.clear()
        self.invalidations += 1
        self.last_invalidation = reason

    def stats(self) -> CacheStats:
        return CacheStats(self.hits, self.misses, len(self._values), self.invalidations, self.last_invalidation)
