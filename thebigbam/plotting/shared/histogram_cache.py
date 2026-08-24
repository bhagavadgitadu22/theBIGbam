"""Process-wide, memory-only cache for immutable filtering histograms."""

from concurrent.futures import Future
from threading import Lock


class SharedHistogramCache:
    """Deduplicate completed and concurrent histogram calculations."""

    def __init__(self):
        self._lock = Lock()
        self._results = {}
        self._inflight = {}

    def get_or_load(self, key, loader):
        with self._lock:
            if key in self._results:
                return self._results[key], "cache_hit"
            future = self._inflight.get(key)
            owner = future is None
            if owner:
                future = Future()
                self._inflight[key] = future
        if not owner:
            return future.result(), "shared_inflight"
        try:
            result = loader()
        except BaseException as exc:
            with self._lock:
                self._inflight.pop(key, None)
                future.set_exception(exc)
            raise
        with self._lock:
            self._inflight.pop(key, None)
            # None represents a failed/empty calculation and is intentionally
            # retryable. Successful immutable results are retained.
            if result is not None:
                self._results[key] = result
            future.set_result(result)
        return result, "calculated"


SERVER_HISTOGRAM_CACHE = SharedHistogramCache()
