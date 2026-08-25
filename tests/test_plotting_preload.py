from pathlib import Path

import duckdb

from thebigbam.plotting.controls.base import build_controls
from thebigbam.plotting.models.preload import PreloadedPlotData
from thebigbam.plotting.models.session import CurrentPlotState, PlotSessionContext
from thebigbam.plotting.repositories.filtering import FilteringRepository
from thebigbam.plotting.repositories.preload import PreloadRepository
from thebigbam.plotting.repositories.summary import SummaryRepository
from thebigbam.plotting.services.filtering import FilteringAvailabilityService
from thebigbam.plotting.services.summary import column_scales


def _create_preload_database(path):
    connection = duckdb.connect(str(path))
    connection.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR, Contig_length INTEGER)")
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    connection.execute("CREATE TABLE Coverage (Contig_id INTEGER, Sample_id INTEGER)")
    connection.execute(
        "CREATE TABLE Variable (Feature_table_name VARCHAR, Module VARCHAR, Subplot VARCHAR, "
        "Module_order INTEGER, Title VARCHAR, Help VARCHAR, Variable_name VARCHAR, Encoding VARCHAR)"
    )
    connection.execute("INSERT INTO Contig VALUES (1, 'c1', 100), (2, 'c2', 200)")
    connection.execute("INSERT INTO Sample VALUES (1, 's1'), (2, 's2')")
    connection.execute("INSERT INTO Coverage VALUES (1, 1), (2, 1), (2, 2)")
    connection.execute(
        "INSERT INTO Variable VALUES "
        "('Feature_blob', 'Coverage', 'Primary alignments', 1, 'Coverage', 'help', 'primary', 'dense')"
    )
    connection.close()


def test_preload_repository_returns_typed_data_and_closes_connection(tmp_path):
    path = tmp_path / "preload.db"
    _create_preload_database(path)
    repository = PreloadRepository(
        str(path),
        scale_initializer=lambda _connection: None,
        filter_encoder=lambda _connection: {"dense": 1},
        filtering_metadata_loader=lambda _path: {"Contig": {}},
        mag_mode_loader=lambda _connection: False,
        mag_map_loader=lambda _connection: ({}, {}),
    )

    data = repository.load()

    assert isinstance(data, PreloadedPlotData)
    assert data.contigs == ("c1", "c2")
    assert data.samples == ("s1", "s2")
    assert data["module_names"] == ("Coverage",)
    assert data.total_coverage_count == 3
    assert data.filter_encode == {"dense": 1}
    assert repository.query_count > 0
    writable = duckdb.connect(str(path))
    writable.close()


def test_availability_service_caches_until_filter_revision_changes(tmp_path):
    path = tmp_path / "availability.db"
    _create_preload_database(path)
    connection = duckdb.connect(str(path), read_only=True)
    repository = FilteringRepository(connection)
    service = FilteringAvailabilityService(repository)

    assert service.contigs_for_sample(1) == ("c1", "c2")
    assert service.contigs_for_sample(1) == ("c1", "c2")
    assert repository.query_count == 1
    assert service.samples_for_contig(2) == ("s1", "s2")
    assert service.count_contigs_for_sample(1) == 2
    assert service.count_samples_for_contig(2) == 2
    service.invalidate()
    assert service.contigs_for_sample(1) == ("c1", "c2")
    assert repository.query_count == 5
    connection.close()


def test_filtered_availability_queries_are_repository_backed_and_cached(tmp_path):
    path = tmp_path / "filtered-availability.db"
    _create_preload_database(path)
    connection = duckdb.connect(str(path))
    connection.execute("CREATE TABLE MAG (MAG_id INTEGER, MAG_name VARCHAR)")
    connection.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    connection.execute("INSERT INTO MAG VALUES (1, 'm1'), (2, 'm2')")
    connection.execute("INSERT INTO MAG_contigs_association VALUES (1, 1), (1, 2), (2, 2)")
    repository = FilteringRepository(connection)
    service = FilteringAvailabilityService(repository)
    filter_sql = "SELECT Contig_id, Sample_id FROM Coverage WHERE Sample_id = ?"
    parameters = [1]

    assert service.filtered_contigs(filter_sql, parameters, mag_name="m1") == ("c1", "c2")
    assert service.filtered_samples(filter_sql, parameters, mag_name="m2") == ("s1",)
    assert service.filtered_mags(filter_sql, parameters, sample_id=1) == ("m1", "m2")
    assert service.count_filtered_contigs(filter_sql, parameters, mag_name="m1") == 2
    assert service.count_filtered_samples(filter_sql, parameters, mag_name="m1") == 1
    query_count = repository.query_count
    assert service.count_filtered_contigs(filter_sql, parameters, mag_name="m1") == 2
    assert service.count_filtered_samples(filter_sql, parameters, mag_name="m1") == 1
    assert repository.query_count == query_count
    connection.close()


def test_session_context_groups_mutable_session_dependencies():
    state = CurrentPlotState(contig="c1", data_xstart=1, data_xend=100)
    context = PlotSessionContext(None, None, {}, None, state)  # type: ignore[arg-type]
    assert context.plot_state.contig == "c1"


def test_server_delegates_preload_and_unfiltered_availability_queries():
    source = Path("thebigbam/plotting/application/composition_root.py").read_text(encoding="utf-8")
    source += Path("thebigbam/plotting/application/availability.py").read_text(encoding="utf-8")
    public_preload = source[source.index("def preload_db_data") : source.index("def create_layout")]

    assert "PreloadRepository(db_path)" in public_preload
    assert "repository.load()" in public_preload
    assert "availability_service.contigs_for_sample" in source
    assert "availability_service.samples_for_contig" in source
    assert "availability_service.mags_for_sample" in source
    assert ".execute(" not in Path("thebigbam/plotting/application/availability.py").read_text(encoding="utf-8")
    assert "def _legacy_preload_db_data" not in source
    assert "source_table_map" not in source
    assert "get_pairs_for_condition" not in source
    assert "from ..controls.base import build_controls" in source


def test_server_entry_point_contains_no_sql_or_widget_construction():
    source = Path("thebigbam/plotting/start_bokeh_server.py").read_text(encoding="utf-8")
    assert len(source.splitlines()) <= 100
    assert ".execute(" not in source
    assert "pn." not in source
    assert "Checkbox" not in source


def test_core_control_factory_uses_only_preloaded_data(tmp_path):
    path = tmp_path / "controls.db"
    _create_preload_database(path)
    data = PreloadRepository(
        str(path),
        scale_initializer=lambda _connection: None,
        filter_encoder=lambda _connection: {},
        filtering_metadata_loader=lambda _path: {},
        mag_mode_loader=lambda _connection: False,
        mag_map_loader=lambda _connection: ({}, {}),
    ).load()

    widgets = build_controls(data)

    assert widgets["contig_select"].options == ["c1", "c2"]
    assert widgets["sample_select"].options == ["s1", "s2"]
    assert widgets["view_radio"].labels == ["MAG VIEW", "CONTIG VIEW"]
    assert len(widgets["variables_widgets_one"]) == 1
    assert ".execute(" not in Path("thebigbam/plotting/controls/base.py").read_text(encoding="utf-8")


def test_column_scales_are_scoped_to_each_database_connection():
    first = duckdb.connect(":memory:")
    second = duckdb.connect(":memory:")
    for connection, scale in ((first, 10), (second, 100)):
        connection.execute("CREATE TABLE Column_scales (Column_name VARCHAR, Scale DOUBLE, Feature_name VARCHAR)")
        connection.execute("INSERT INTO Column_scales VALUES ('Depth', ?, 'Coverage')", [scale])

    assert column_scales(SummaryRepository(first)) == {"Depth": 10.0}
    assert column_scales(SummaryRepository(second)) == {"Depth": 100.0}
    first.close()
    second.close()
