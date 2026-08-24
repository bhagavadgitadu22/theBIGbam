"""Database-backed metadata and distribution queries for filtering controls."""

from __future__ import annotations

from ...database.database_getters import (
    resolve_column_null_stats,
    resolve_distinct_values,
    resolve_histogram_bins,
    resolve_value_counts,
    search_distinct_values,
)


class FilterMetadataRepository:
    """Expose filtering data as plain Python values behind one repository API."""

    def __init__(self, db_path, metadata, *, enable_timing=False):
        self.db_path = db_path
        self.metadata = metadata
        self.enable_timing = enable_timing

    def distinct_values(self, category, column):
        return resolve_distinct_values(
            self.db_path,
            self.metadata,
            category,
            column,
            enable_timing=self.enable_timing,
        )

    def search_distinct_values(self, category, column, search_term="", limit=100):
        return search_distinct_values(
            self.db_path,
            self.metadata,
            category,
            column,
            search_term,
            limit,
        )

    def histogram_bins(self, category, column, *, n_bins=50, log_mode=False, scale=None):
        return resolve_histogram_bins(
            self.db_path,
            self.metadata,
            category,
            column,
            n_bins=n_bins,
            log_mode=log_mode,
            scale=scale,
            enable_timing=self.enable_timing,
        )

    def value_counts(self, category, column):
        return resolve_value_counts(
            self.db_path,
            self.metadata,
            category,
            column,
            enable_timing=self.enable_timing,
        )

    def column_null_stats(self, category, column):
        return resolve_column_null_stats(self.db_path, self.metadata, category, column)
