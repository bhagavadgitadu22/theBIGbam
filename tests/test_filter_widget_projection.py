from types import SimpleNamespace

from thebigbam.plotting.controls.filter_projection import FilterWidgetProjection, expression_from_widgets
from thebigbam.plotting.models.filters import CompiledFilter, FilterResult


def _widget(value):
    return SimpleNamespace(value=value)


def _row(value, connector=None):
    row = {
        "input_ref": {"widget": _widget(value), "is_panel": False},
        "category_select": _widget("Contig"),
        "subcategory_select": _widget("Length"),
        "comparison_select": _widget(">"),
    }
    if connector is not None:
        row["and_div"] = _widget(connector)
    return row


def test_expression_from_widgets_ignores_empty_rows_and_preserves_connectors():
    expression = expression_from_widgets(
        [{"rows": [_row(10), _row(20, "OR")]}, {"rows": [_row(""), _row(30)]}],
        [_widget("AND")],
    )

    assert len(expression.sections) == 2
    assert expression.sections[0].connectors == ("OR",)
    assert expression.connectors == ("AND",)


def test_projection_caches_result_until_invalidated():
    class Service:
        def __init__(self):
            self.calls = 0
            self.invalidations = 0

        def evaluate(self, _expression):
            self.calls += 1
            return FilterResult(CompiledFilter("SELECT 1", (7,)), 4, 3, 2, 1)

        def invalidate(self):
            self.invalidations += 1

    service = Service()
    projection = FilterWidgetProjection(service, [{"rows": [_row(10)]}], [])

    assert projection.evaluate()["count_pairs"] == 4
    assert projection.evaluate()["params"] == [7]
    assert service.calls == 1
    projection.invalidate()
    projection.evaluate()
    assert service.calls == 2
    assert service.invalidations == 1
