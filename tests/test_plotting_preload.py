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


def test_initial_browser_selector_payload_is_capped_at_one_hundred(tmp_path):
    path = tmp_path / "large-controls.db"
    _create_preload_database(path)
    connection = duckdb.connect(str(path))
    connection.executemany(
        "INSERT INTO Contig VALUES (?, ?, 100)",
        [(index, f"contig-{index:03d}") for index in range(3, 153)],
    )
    connection.executemany(
        "INSERT INTO Sample VALUES (?, ?)",
        [(index, f"sample-{index:03d}") for index in range(3, 153)],
    )
    connection.close()
    data = PreloadRepository(
        str(path),
        scale_initializer=lambda _connection: None,
        filter_encoder=lambda _connection: {},
        filtering_metadata_loader=lambda _path: {},
        mag_mode_loader=lambda _connection: False,
        mag_map_loader=lambda _connection: ({}, {}),
    ).load()

    widgets = build_controls(data)

    assert len(data.contigs) == 152
    assert len(data.samples) == 152
    assert len(widgets["contig_select"].options) == 100
    assert len(widgets["sample_select"].options) == 100


def test_unscoped_name_searches_are_bounded_and_prioritize_matching_modes(tmp_path):
    path = tmp_path / "search.db"
    _create_preload_database(path)
    connection = duckdb.connect(str(path))
    connection.execute("INSERT INTO Contig VALUES (3, 'target', 100), (4, 'target-prefix', 100), (5, 'x-target', 100)")
    repository = FilteringRepository(connection)

    assert repository.contigs("target") == ("target", "target-prefix", "x-target")
    assert len(repository.contigs()) <= 100
    connection.close()


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
