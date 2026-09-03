"""Shared server-backed autocomplete for categorical database values."""

from __future__ import annotations

from typing import Any, Callable

from ..shared.defaults import AUTOCOMPLETE_LIMIT
from .searchable_select import SearchableSelect, decode_search_request


def build_distinct_value_select(
    metadata_service: Any,
    category: str,
    column: str,
    *,
    name: str,
    on_value: Callable[[], None] | None = None,
    **kwargs: Any,
) -> SearchableSelect:
    """Build the canonical bounded, live database autocomplete."""
    search_values = getattr(metadata_service, "search_distinct_values", None)
    initial_values = (
        list(search_values(category, column, "", limit=AUTOCOMPLETE_LIMIT))
        if search_values is not None
        else []
    )
    widget = SearchableSelect(
        name=name,
        value="",
        options=initial_values,
        placeholder="Search...",
        server_search=True,
        min_search_chars=0,
        **kwargs,
    )

    def search(_event: Any) -> None:
        request_nonce, query = decode_search_request(widget.search_request)
        values: list[str] = []
        if len(query) >= widget.min_search_chars:
            values = (
                list(search_values(category, column, query, limit=AUTOCOMPLETE_LIMIT))
                if search_values is not None
                else []
            )
        widget.options = values
        widget.search_result_query = query
        widget.search_result_nonce = request_nonce

    widget.param.watch(search, "search_request")
    if on_value is not None:
        widget.param.watch(lambda _event: on_value(), "value")
    return widget
