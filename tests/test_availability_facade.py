from types import SimpleNamespace

import pytest

from thebigbam.plotting.application.availability_facade import AvailabilityFacade


def test_widget_projection_preserves_valid_value_and_clears_invalid_value():
    widget = SimpleNamespace(options=["a"], value="a")
    AvailabilityFacade.update_widget(widget, ["a", "b"])
    assert widget.value == "a"

    AvailabilityFacade.update_widget(widget, ["b"])
    assert widget.value == ""


def test_search_projection_forces_options_event_even_for_same_values():
    events = []
    widget = SimpleNamespace(
        options=["a"],
        value="a",
        param=SimpleNamespace(trigger=events.append),
    )
    AvailabilityFacade.push_search(widget, ["a"])
    assert events == ["options"]


def test_scope_projection_replaces_equal_frontend_pool_and_clears_search():
    widget = SimpleNamespace(
        options=["only"],
        value="only",
        search_query="old query",
        scope_nonce=4,
    )

    AvailabilityFacade.replace_widget(widget, ["only"])

    assert widget.options == ["only"]
    assert widget.value == "only"
    assert widget.search_query == ""
    assert widget.scope_nonce == 5


def test_facade_requires_attachment_and_delegates_after_attachment():
    facade = AvailabilityFacade()
    with pytest.raises(RuntimeError, match="not been attached"):
        facade.compute_contigs()

    controller = SimpleNamespace(compute_contigs=lambda term: [term])
    facade.attach(controller)
    assert facade.compute_contigs("abc") == ["abc"]
