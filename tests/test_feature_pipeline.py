import duckdb
from bokeh.models import HoverTool, Range1d

from thebigbam.plotting.renderers.features import FeatureRenderer
from thebigbam.plotting.repositories.features import FeatureRepository, split_features_by_storage
from thebigbam.plotting.services.features import add_valid_termini, apply_minimum_frequency


def test_valid_termini_query_keeps_both_passing_dtr_copies_and_respects_status():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Phage_mechanisms (Packaging_id INTEGER, Contig_id INTEGER, Sample_id INTEGER)")
    connection.execute(
        "CREATE TABLE Phage_termini (Packaging_id INTEGER, Status VARCHAR, Center INTEGER, "
        "Passed_PoissonTest BOOLEAN, Passed_ClippingTest BOOLEAN)"
    )
    connection.execute("INSERT INTO Phage_mechanisms VALUES (1, 7, 3), (2, 7, 4)")
    connection.execute(
        "INSERT INTO Phage_termini VALUES "
        "(1, 'start', 10, true, true), (1, 'start', 90, true, true), "
        "(1, 'start', 20, true, false), (1, 'end', 30, true, true), "
        "(2, 'start', 40, true, true)"
    )
    repository = FeatureRepository(connection)
    passing = repository.valid_termini((7,), (3,))
    assert passing[(7, 3, "start")] == {10, 90}
    assert passing[(7, 3, "end")] == {30}
    before = repository.query_count
    assert repository.valid_termini((7,), (3,)) is passing
    assert repository.query_count == before


def test_valid_termini_marks_exact_bars_and_aggregated_bins():
    raw = add_valid_termini({"x": [10, 20, 90]}, {10, 90})
    assert raw["terminus_valid"] == ["Yes", "No", "Yes"]
    zoom = add_valid_termini(
        {"x": [15, 25], "first_pos": [11, 21], "last_pos": [20, 30], "width": [10, 10]},
        {20},
    )
    assert zoom["terminus_valid"] == ["Yes", "No"]


def test_metadata_batch_uses_one_query_and_preserves_storage_split():
    connection = duckdb.connect(":memory:")
    connection.execute(
        'CREATE TABLE Variable (Subplot VARCHAR, "Type" VARCHAR, Color VARCHAR, Alpha REAL, '
        'Fill_alpha REAL, "Size" REAL, Title VARCHAR, Feature_table_name VARCHAR, Module_order INTEGER)'
    )
    connection.executemany(
        "INSERT INTO Variable VALUES (?, 'curve', '#000', 1, 0.2, 1, ?, ?, ?)",
        [("GC", "GC", "Contig_blob", 1), ("Coverage", "Coverage", "Feature_blob", 2)],
    )
    repository = FeatureRepository(connection)
    metadata = repository.metadata_batch(["GC", "Coverage"])
    assert repository.query_count == 1
    assert split_features_by_storage(metadata, ["Coverage", "GC", "Missing"]) == (
        ["GC"],
        ["Coverage", "Missing"],
    )


def test_minimum_frequency_is_absolute_pure_and_handles_zero_series():
    original = [0.2, 0.5, 0.8]
    assert apply_minimum_frequency(original, 0.5) == [0, 0.5, 0.8]
    assert original == [0.2, 0.5, 0.8]
    assert apply_minimum_frequency([0, 0], 0.5) == [0, 0]


def test_feature_renderer_handles_curves_bars_and_relative_scale():
    common = {"color": "#123456", "alpha": 1.0, "fill_alpha": 0.2, "size": 1, "title": "Track"}
    curve = {**common, "type": "curve", "x": [1, 2], "y": [0.25, 0.75], "is_relative_scaled": True}
    bars = {
        **common,
        "type": "bars",
        "x": [1, 3],
        "y": [0.5, 0.75],
        "width": [1, 2],
        "first_pos": [1, 2],
        "last_pos": [1, 4],
        "is_relative_scaled": True,
    }
    plot = FeatureRenderer().render([curve, bars], 100, Range1d(1, 4), sample_title="sample")
    assert plot is not None
    assert plot.y_range.end == 1
    assert len(plot.renderers) >= 3


def test_termini_hover_only_appears_on_annotated_bars():
    common = {"type": "bars", "color": "#123456", "alpha": 1, "size": 1, "title": "Read starts"}
    bars = {**common, "x": [10, 20], "y": [5, 3], "terminus_valid": ["Yes", "No"]}
    plot = FeatureRenderer().render([bars], 100, Range1d(1, 30))
    hover = next(tool for tool in plot.tools if isinstance(tool, HoverTool))
    assert ("Valid terminus", "@terminus_valid") in hover.tooltips
    assert plot.renderers[0].data_source.data["terminus_valid"] == ["Yes", "No"]

    zoom = {**bars, "width": [10, 10], "first_pos": [5, 15], "last_pos": [14, 24]}
    plot = FeatureRenderer().render([zoom], 100, Range1d(1, 30))
    hover = next(tool for tool in plot.tools if isinstance(tool, HoverTool))
    assert ("Contains valid terminus", "@terminus_valid") in hover.tooltips
