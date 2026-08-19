from dataclasses import FrozenInstanceError
from pathlib import Path

import duckdb
import pytest

from thebigbam.plotting.models.composition import (
    CompositionDisplay,
    MagCompositionRequest,
    SingleSampleCompositionRequest,
    TrackSelection,
)
from thebigbam.plotting.repositories.mag_overview import MagOverviewRepository
from thebigbam.plotting.services.mag_overview import MagOverviewService


def test_composition_requests_are_immutable_and_group_related_options():
    request = SingleSampleCompositionRequest(
        contig_name="c1",
        sample_name="s1",
        xstart=1,
        xend=100,
        tracks=TrackSelection(features=("Coverage",), plot_sequence=True),
        display=CompositionDisplay(subplot_height=120, max_base_resolution=10_000),
    )

    assert request.tracks.features == ("Coverage",)
    assert request.display.subplot_height == 120
    with pytest.raises(FrozenInstanceError):
        request.contig_name = "c2"


def test_mag_request_has_independent_track_and_display_contracts():
    request = MagCompositionRequest(
        mag_name="m1",
        sample_name=None,
        is_all=True,
        allowed_samples=frozenset({"s1", "s2"}),
        tracks=TrackSelection(features=("Primary alignments",)),
    )

    assert request.is_all is True
    assert request.allowed_samples == frozenset({"s1", "s2"})


def test_mag_overview_repository_retrieves_and_limits_annotation_dots(tmp_path):
    connection = duckdb.connect(str(tmp_path / "overview.db"))
    connection.execute(
        "CREATE TABLE Contig_annotation_core "
        '(Annotation_id INTEGER, Contig_id INTEGER, "Start" INTEGER, "End" INTEGER, Type VARCHAR)'
    )
    connection.execute("CREATE TABLE Annotation_qualifier (Annotation_id INTEGER, Key VARCHAR, Value VARCHAR)")
    connection.execute("INSERT INTO Contig_annotation_core VALUES (1, 10, 5, 9, 'CDS'), (2, 11, 2, 4, 'CDS')")
    connection.execute("INSERT INTO Annotation_qualifier VALUES (1, 'product', 'capsid'), (2, 'product', 'capsid')")

    repository = MagOverviewRepository(connection)
    xs, colors, total = MagOverviewService(repository).annotation_dots(
        ((10, 100, 0), (11, 50, 100)),
        ({"qualifier_key": "product", "match_mode": "exact", "value": "capsid", "color": "red"},),
        max_dots=1,
    )

    assert xs == [7]
    assert colors == ["red"]
    assert total == 2
    assert repository.query_count == 2
    connection.close()


def test_legacy_plot_function_adapters_are_removed():
    source = Path("thebigbam/plotting/application/sample_mag_pipeline.py").read_text(encoding="utf-8")
    assert "def generate_bokeh_plot_per_sample" not in source
    assert "def generate_bokeh_plot_mag_view" not in source
    assert "def compute_mag_track_dots" not in source
    assert "def make_bokeh_mag_track" not in source


def test_server_uses_typed_composers_and_mag_repository_is_rendering_free():
    controller = Path("thebigbam/plotting/application/apply_controller.py").read_text(encoding="utf-8")
    handlers = Path("thebigbam/plotting/application/apply_handlers/rendering.py").read_text(encoding="utf-8")
    repository = Path("thebigbam/plotting/repositories/mag_overview.py").read_text(encoding="utf-8")

    assert "build_single_sample_plot(conn, request)" in handlers
    assert "build_mag_plot(conn, request)" in handlers
    assert "compose_" not in controller
    assert "generate_bokeh_plot_per_sample(" not in handlers
    assert "generate_bokeh_plot_mag_view(" not in handlers
    assert "genome_tracks.mag_members(" not in handlers
    assert "from bokeh" not in repository.lower()
    assert "import bokeh" not in repository.lower()
    assert "import panel" not in repository.lower()


def test_composers_are_sql_free_and_services_are_bokeh_free():
    composers = Path("thebigbam/plotting/composers/layout.py").read_text(encoding="utf-8")
    repository = Path("thebigbam/plotting/repositories/composition.py").read_text(encoding="utf-8")
    service = Path("thebigbam/plotting/services/composition.py").read_text(encoding="utf-8")

    assert ".execute(" not in composers
    assert "SELECT " not in composers
    assert "from bokeh" not in repository.lower()
    assert "import bokeh" not in repository.lower()
    assert "from bokeh" not in service.lower()
    assert "import bokeh" not in service.lower()
