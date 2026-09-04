from dataclasses import FrozenInstanceError

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
from thebigbam.plotting.shared.data_cache import SessionDataCache


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


def test_mag_overview_cache_reuses_queries_when_only_color_changes(tmp_path):
    connection = duckdb.connect(str(tmp_path / "overview-cache.db"))
    connection.execute(
        "CREATE TABLE Contig_annotation_core "
        '(Annotation_id INTEGER, Contig_id INTEGER, "Start" INTEGER, "End" INTEGER, Type VARCHAR)'
    )
    connection.execute("CREATE TABLE Annotation_qualifier (Annotation_id INTEGER, Key VARCHAR, Value VARCHAR)")
    connection.execute("INSERT INTO Contig_annotation_core VALUES (1, 10, 5, 9, 'CDS')")
    connection.execute("INSERT INTO Annotation_qualifier VALUES (1, 'product', 'capsid')")
    repository = MagOverviewRepository(connection)
    service = MagOverviewService(repository, cache=SessionDataCache())
    members = ((10, 100, 0),)
    rule = {"qualifier_key": "product", "match_mode": "exact", "value": "capsid", "color": "red"}

    assert service.annotation_dots(members, (rule,))[1] == ["red"]
    assert service.annotation_dots(members, ({**rule, "color": "blue"},))[1] == ["blue"]
    assert repository.query_count == 2
    connection.close()
