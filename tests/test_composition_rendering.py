from types import SimpleNamespace

import pytest
from bokeh.models import Range1d
from bokeh.plotting import figure

from thebigbam.plotting.composers.layout import assemble_grid
from thebigbam.plotting.renderers.composition import apply_per_feature_y_ranges, apply_primary_relative_y_range


def test_grid_assembly_filters_missing_optional_tracks():
    plot = figure()
    grid = assemble_grid([None, plot])
    assert plot in grid.references()
    with pytest.raises(ValueError, match="nothing"):
        assemble_grid([None], empty_message="nothing")


def test_subplot_titles_setting_hides_feature_legends():
    visible_plot = figure()
    visible_plot.line([1, 2], [1, 2], legend_label="Coverage")
    assemble_grid([visible_plot])
    assert visible_plot.legend[0].visible

    hidden_plot = figure()
    hidden_plot.line([1, 2], [1, 2], legend_label="Coverage")
    assemble_grid([hidden_plot], show_subplot_titles=False)
    assert not hidden_plot.legend[0].visible
    assert hidden_plot.renderers


def test_y_range_policies_are_sql_free_and_deterministic():
    first = SimpleNamespace(y_range=None)
    second = SimpleNamespace(y_range=None)
    apply_per_feature_y_ranges([("Coverage", first, 4), ("Coverage", second, 9)])
    assert isinstance(first.y_range, Range1d)
    assert first.y_range.end == second.y_range.end == 9

    relative = SimpleNamespace(y_range=None)
    unrelated = SimpleNamespace(y_range=None)
    apply_primary_relative_y_range([("Clippings", relative, 2), ("GC content", unrelated, 5)], 12)
    assert relative.y_range.end == 12
    assert unrelated.y_range is None
