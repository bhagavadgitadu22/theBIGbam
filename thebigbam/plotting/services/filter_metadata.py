"""Filtering metadata values prepared for UI presentation."""

from __future__ import annotations


class FilterMetadataService:
    """Cache plain distribution results returned by a metadata repository."""

    def __init__(self, repository, histogram_cache=None):
        self.repository = repository
        self._cache = {}
        self.histogram_cache = histogram_cache

    def _cached(self, key, load):
        if key not in self._cache:
            self._cache[key] = load()
        return self._cache[key]

    def distinct_values(self, category, column):
        return self._cached(
            ("distinct", category, column),
            lambda: self.repository.distinct_values(category, column),
        )

    def search_distinct_values(self, category, column, search_term="", limit=100):
        return self._cached(
            ("distinct_search", category, column, search_term, limit),
            lambda: self.repository.search_distinct_values(category, column, search_term, limit),
        )

    def histogram_bins(self, category, column, *, n_bins=50, log_mode=False, scale=None):
        import os
        db_path = os.path.realpath(self.repository.db_path)
        stat = os.stat(db_path)
        db_identity = (db_path, stat.st_size, stat.st_mtime_ns)
        key = (db_identity, category, column, log_mode, scale, n_bins)
        def load():
            return self.repository.histogram_bins(
                category, column, n_bins=n_bins, log_mode=log_mode, scale=scale
            )
        if self.histogram_cache is None:
            return self._cached(("histogram", *key), load)
        result, outcome = self.histogram_cache.get_or_load(key, load)
        if self.repository.enable_timing and outcome != "calculated":
            mode = "log" if log_mode else "linear"
            bins = result.bin_count if result else 0
            sampled = result.sampled_rows if result else 0
            print(
                f"[timing] resolve_histogram_bins({category}.{column}, {mode}): "
                f"{outcome} ({bins} bins, {sampled} rows)",
                flush=True,
            )
        return result

    def value_counts(self, category, column):
        return self._cached(
            ("counts", category, column),
            lambda: self.repository.value_counts(category, column),
        )

    def column_null_stats(self, category, column):
        return self._cached(
            ("nulls", category, column),
            lambda: self.repository.column_null_stats(category, column),
        )
