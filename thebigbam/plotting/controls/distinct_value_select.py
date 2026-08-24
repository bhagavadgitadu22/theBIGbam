"""Shared server-backed autocomplete for categorical database values."""

from __future__ import annotations

from typing import Any, Callable

from .searchable_select import SearchableSelect


def build_distinct_value_select(
    metadata_service: Any,
    category: str,
    column: str,
    *,
    name: str,
    on_value: Callable[[], None] | None = None,
    **kwargs: Any,
) -> SearchableSelect:
    """Build the canonical two-character, live database autocomplete."""
    widget = SearchableSelect(
        name=name,
        value="",
        options=[],
        placeholder="Search...",
        server_search=True,
        min_search_chars=2,
        **kwargs,
    )

    def search(_event: Any) -> None:
        request_nonce = widget.search_nonce
        query = widget.search_query
        values: list[str] = []
        if len(query) >= widget.min_search_chars:
            values = list(metadata_service.search_distinct_values(category, column, query, limit=100))
            if widget.value and widget.value not in values:
                values.insert(0, widget.value)
        widget.options = values
        widget.search_result_query = query
        widget.search_result_nonce = request_nonce

    widget.param.watch(search, "search_nonce")
    if on_value is not None:
        widget.param.watch(lambda _event: on_value(), "value")
    return widget
