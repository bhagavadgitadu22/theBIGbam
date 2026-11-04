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

    def __init__(self, max_entries: int = 256, *, max_retained_bytes: int = 64 * 1024 * 1024) -> None:
        self.max_entries = max_entries
        self._values: OrderedDict[Hashable, Any] = OrderedDict()
        self.max_retained_bytes = max_retained_bytes
        self.retained_bytes = 0
        self._weights: dict[Hashable, int] = {}
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

    def _remove(self, key: Hashable) -> None:
        self._values.pop(key)
        self.retained_bytes -= self._weights.pop(key, 0)

    def discard(self, key: Hashable) -> None:
        if key in self._values:
            self._remove(key)

    def put(self, key: Hashable, value: Any, *, retained_bytes: int = 0) -> bool:
        if retained_bytes > self.max_retained_bytes:
            return False
        if key in self._values:
            self._remove(key)
        while self.retained_bytes + retained_bytes > self.max_retained_bytes:
            oldest = next(key for key in self._values if self._weights.get(key, 0))
            self._remove(oldest)
        self._values[key] = value
        if retained_bytes:
            self._weights[key] = retained_bytes
            self.retained_bytes += retained_bytes
        self._values.move_to_end(key)
        while len(self._values) > self.max_entries:
            self._remove(next(iter(self._values)))
        return True

    def invalidate(self, reason: str, preserve_prefixes: tuple[Hashable, ...] = ()) -> None:
        if preserve_prefixes:
            self._values = OrderedDict(
                (key, value)
                for key, value in self._values.items()
                if isinstance(key, tuple) and key and key[0] in preserve_prefixes
            )
        else:
            self._values.clear()
        self._weights = {key: weight for key, weight in self._weights.items() if key in self._values}
        self.retained_bytes = sum(self._weights.values())
        self.invalidations += 1
        self.last_invalidation = reason

    def stats(self) -> CacheStats:
        return CacheStats(self.hits, self.misses, len(self._values), self.invalidations, self.last_invalidation)
