"""Revision-aware cache for filtering and autocomplete availability."""

from __future__ import annotations

from typing import Callable, Hashable

from ..models.filters import FilterExpression, FilterResult
from ..repositories.filtering import FilteringRepository
from .filter_query import FilterQueryBuilder


class FilteringAvailabilityService:
    def __init__(self, repository: FilteringRepository) -> None:
        self.repository = repository
        self.revision = 0
        self._cache: dict[tuple[Hashable, ...], object] = {}

    def invalidate(self) -> None:
        self.revision += 1
        self._cache.clear()

    def _cached(self, key: tuple[Hashable, ...], load: Callable[[], object]):
        revision_key = (self.revision, *key)
        if revision_key not in self._cache:
            self._cache[revision_key] = load()
        return self._cache[revision_key]

    def contigs_for_sample(self, sample_id: int, search_term: str = "", preserve: str = "") -> tuple[str, ...]:
        return self._cached(
            ("contigs", sample_id, search_term, preserve),
            lambda: self.repository.contigs_for_sample(sample_id, search_term, preserve),
        )

    def samples_for_contig(self, contig_id: int, search_term: str = "") -> tuple[str, ...]:
        return self._cached(
            ("samples", contig_id, search_term),
            lambda: self.repository.samples_for_contig(contig_id, search_term),
        )

    def mags_for_sample(self, sample_id: int, search_term: str = "") -> tuple[str, ...]:
        return self._cached(
            ("mags", sample_id, search_term),
            lambda: self.repository.mags_for_sample(sample_id, search_term),
        )

    def count_contigs_for_sample(self, sample_id: int) -> int:
        return self._cached(
            ("contig_count", sample_id),
            lambda: self.repository.count_contigs_for_sample(sample_id),
        )

    def count_samples_for_contig(self, contig_id: int) -> int:
        return self._cached(
            ("sample_count", contig_id),
            lambda: self.repository.count_samples_for_contig(contig_id),
        )

    def filtered_contigs(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_contigs", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_contigs(sql, parameters, **scope))

    def filtered_samples(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_samples", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_samples(sql, parameters, **scope))

    def filtered_mags(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_mags", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_mags(sql, parameters, **scope))

    def filtered_mag_counts(self, sql: str, parameters: list, mag_name: str) -> tuple[int, int]:
        return self._cached(
            ("filtered_mag_counts", sql, tuple(parameters), mag_name),
            lambda: self.repository.filtered_mag_counts(sql, parameters, mag_name),
        )


class FilterExpressionService:
    """Compile and evaluate an immutable expression once per revision."""

    def __init__(self, repository: FilteringRepository, builder: FilterQueryBuilder, *, has_mags: bool) -> None:
        self.repository = repository
        self.builder = builder
        self.has_mags = has_mags
        self._expression: FilterExpression | None = None
        self._result: FilterResult | None = None

    def invalidate(self) -> None:
        self._expression = None
        self._result = None

    def evaluate(self, expression: FilterExpression) -> FilterResult | None:
        if expression == self._expression:
            return self._result
        compiled = self.builder.compile(expression)
        self._expression = expression
        self._result = self.repository.evaluate(compiled, has_mags=self.has_mags) if compiled else None
        return self._result
