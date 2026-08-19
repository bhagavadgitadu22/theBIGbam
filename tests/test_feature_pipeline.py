import inspect

import duckdb
from bokeh.models import Range1d

from thebigbam.plotting.renderers.features import FeatureRenderer
from thebigbam.plotting.repositories.features import FeatureRepository, split_features_by_storage
from thebigbam.plotting.services.features import apply_relative_threshold


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


def test_relative_threshold_is_pure_and_handles_zero_series():
    original = [1.0, 5.0, 10.0]
    assert apply_relative_threshold(original, 0.5) == [0, 5.0, 10.0]
    assert original == [1.0, 5.0, 10.0]
    assert apply_relative_threshold([0, 0], 0.5) == [0, 0]


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


def test_renderer_and_service_boundaries_have_no_database_or_ui_coupling():
    import thebigbam.plotting.renderers.all_samples as all_renderer
    import thebigbam.plotting.renderers.features as renderer
    import thebigbam.plotting.services.features as service

    renderer_source = inspect.getsource(renderer).lower()
    service_source = inspect.getsource(service).lower()
    assert "execute(" not in renderer_source
    assert "duckdb" not in renderer_source
    assert "bokeh" not in service_source
    assert "panel" not in service_source
    assert "plotting_data_per_sample" not in inspect.getsource(all_renderer)
