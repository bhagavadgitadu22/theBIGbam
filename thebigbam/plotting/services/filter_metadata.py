"""Filtering metadata values prepared for UI presentation."""

from __future__ import annotations


class FilterMetadataService:
    """Cache plain distribution results returned by a metadata repository."""

    def __init__(self, repository):
        self.repository = repository
        self._cache = {}

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

    def histogram_bins(self, category, column, *, n_bins=50, log_mode=False):
        return self._cached(
            ("histogram", category, column, n_bins, log_mode),
            lambda: self.repository.histogram_bins(category, column, n_bins=n_bins, log_mode=log_mode),
        )

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
