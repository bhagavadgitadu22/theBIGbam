"""Revision-aware cache for filtering and autocomplete availability."""

from __future__ import annotations

from collections import OrderedDict
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

    def contigs(self, search_term: str = "", preserve: str = "") -> tuple[str, ...]:
        return self._cached(
            ("all_contigs", search_term, preserve), lambda: self.repository.contigs(search_term, preserve)
        )

    def samples(self, search_term: str = "", preserve: str = "") -> tuple[str, ...]:
        return self._cached(
            ("all_samples", search_term, preserve), lambda: self.repository.samples(search_term, preserve)
        )

    def mags(self, search_term: str = "", preserve: str = "") -> tuple[str, ...]:
        return self._cached(
            ("all_mags", search_term, preserve), lambda: self.repository.mags(search_term, preserve)
        )

    def _cached(self, key: tuple[Hashable, ...], load: Callable[[], object]):
        revision_key = (self.revision, *key)
        if revision_key not in self._cache:
            self._cache[revision_key] = load()
        return self._cache[revision_key]

    def contigs_for_sample(
        self, sample_id: int, search_term: str = "", preserve: str = "", mag_name: str | None = None
    ) -> tuple[str, ...]:
        return self._cached(
            ("contigs", sample_id, search_term, preserve, mag_name),
            lambda: self.repository.contigs_for_sample(sample_id, search_term, preserve, mag_name),
        )

    def samples_for_contig(
        self, contig_id: int, search_term: str = "", preserve: str = ""
    ) -> tuple[str, ...]:
        return self._cached(
            ("samples", contig_id, search_term, preserve),
            lambda: self.repository.samples_for_contig(contig_id, search_term, preserve),
        )

    def mags_for_sample(
        self, sample_id: int, search_term: str = "", preserve: str = ""
    ) -> tuple[str, ...]:
        return self._cached(
            ("mags", sample_id, search_term, preserve),
            lambda: self.repository.mags_for_sample(sample_id, search_term, preserve),
        )

    def count_contigs_for_sample(self, sample_id: int, mag_name: str | None = None) -> int:
        return self._cached(
            ("contig_count", sample_id, mag_name),
            lambda: self.repository.count_contigs_for_sample(sample_id, mag_name),
        )

    def count_samples_for_contig(self, contig_id: int) -> int:
        return self._cached(
            ("sample_count", contig_id),
            lambda: self.repository.count_samples_for_contig(contig_id),
        )

    def filtered_contigs(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_contigs", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_contigs(sql, parameters, **scope))

    def count_filtered_contigs(self, sql, parameters, **scope) -> int:
        key = ("count_filtered_contigs", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.count_filtered_contigs(sql, parameters, **scope))

    def filtered_samples(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_samples", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_samples(sql, parameters, **scope))

    def count_filtered_samples(self, sql, parameters, **scope) -> int:
        key = ("count_filtered_samples", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.count_filtered_samples(sql, parameters, **scope))

    def filtered_mags(self, sql, parameters, **scope) -> tuple[str, ...]:
        key = ("filtered_mags", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.filtered_mags(sql, parameters, **scope))

    def count_filtered_mags(self, sql, parameters, **scope) -> int:
        key = ("count_filtered_mags", sql, tuple(parameters), *sorted(scope.items()))
        return self._cached(key, lambda: self.repository.count_filtered_mags(sql, parameters, **scope))

class FilterExpressionService:
    """Compile and evaluate an immutable expression once per revision."""

    def __init__(self, repository: FilteringRepository, builder: FilterQueryBuilder, *, has_mags: bool) -> None:
        self.repository = repository
        self.builder = builder
        self.has_mags = has_mags
        self._expression: FilterExpression | None = None
        self._result: FilterResult | None = None
        self._cache: OrderedDict[FilterExpression, FilterResult | None] = OrderedDict()

    def invalidate(self) -> None:
        self._expression = None
        self._result = None

    def evaluate(self, expression: FilterExpression) -> FilterResult | None:
        if expression == self._expression:
            return self._result
        if expression in self._cache:
            cached = self._cache.pop(expression)
            if self.repository.reuse_filter_result(cached):
                self._expression = expression
                self._result = cached
                self._cache[expression] = cached
                return self._result
        compiled = self.builder.compile(expression)
        self._expression = expression
        self._result = self.repository.evaluate(compiled, has_mags=self.has_mags) if compiled else None
        self._cache[expression] = self._result
        while len(self._cache) > self.repository.result_cache_size:
            self._cache.popitem(last=False)
        return self._result
