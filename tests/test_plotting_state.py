import io
import time
from pathlib import Path

import duckdb
import pytest

from thebigbam.plotting.application.state import (
    AvailabilitySnapshot,
    PlotController,
    PlotUiState,
    SampleScope,
    enforce_all_sample_variable,
)
from thebigbam.plotting.models.sample_scope import (
    apply_changed_attributes,
    visibility_for_sample_scope,
)
from thebigbam.plotting.shared.diagnostics import PlotDiagnostics


@pytest.fixture
def plotting_db(tmp_path):
    """Small deterministic contig/MAG database for plotting-controller tests."""
    path = tmp_path / "plotting.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name VARCHAR, Contig_length INTEGER)")
    conn.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    conn.execute("CREATE TABLE Coverage (Contig_id INTEGER, Sample_id INTEGER)")
    conn.execute("CREATE TABLE MAG (MAG_id INTEGER, MAG_name VARCHAR)")
    conn.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    conn.execute("INSERT INTO Contig VALUES (1, 'c1', 1000), (2, 'c2', 2000), (3, 'orphan', 50)")
    conn.execute("INSERT INTO Sample VALUES (1, 's1'), (2, 's2')")
    conn.execute("INSERT INTO Coverage VALUES (1, 1), (1, 2), (2, 2)")
    conn.execute("INSERT INTO MAG VALUES (1, 'm1')")
    conn.execute("INSERT INTO MAG_contigs_association VALUES (1, 1), (1, 2)")
    conn.close()
    return path


def test_synthetic_plotting_database_models_contig_and_mag_views(plotting_db):
    conn = duckdb.connect(str(plotting_db), read_only=True)
    pairs = conn.execute(
        "SELECT c.Contig_name, s.Sample_name FROM Coverage p "
        "JOIN Contig c USING (Contig_id) JOIN Sample s USING (Sample_id) ORDER BY 1, 2"
    ).fetchall()
    mag_members = conn.execute(
        "SELECT mg.MAG_name, c.Contig_name FROM MAG_contigs_association mca "
        "JOIN MAG mg USING (MAG_id) JOIN Contig c USING (Contig_id) ORDER BY 2"
    ).fetchall()
    conn.close()
    assert pairs == [("c1", "s1"), ("c1", "s2"), ("c2", "s2")]
    assert mag_members == [("m1", "c1"), ("m1", "c2")]


def test_scope_switch_is_io_free_and_preserves_independent_state():
    projections = []
    one = AvailabilitySnapshot(samples=("s1",), contigs=("c1",), mags=("m1",))
    all_samples = AvailabilitySnapshot(samples=("s1", "s2"), contigs=("c1", "c2"), mags=("m1",))
    state = PlotUiState(sample="s1", contig="c1", mag="m1")
    controller = PlotController(
        state,
        {SampleScope.ONE: one, SampleScope.ALL: all_samples},
        project=lambda projected, available: projections.append((projected, available)),
    )

    controller.switch_sample_scope(SampleScope.ALL)
    controller.switch_sample_scope(SampleScope.ONE)

    assert controller.state.sample_scope is SampleScope.ONE
    assert controller.state.sample == "s1"
    assert controller.transition_count == 2
    assert [state.sample_scope for state, _ in projections] == [SampleScope.ALL, SampleScope.ONE]


def test_scope_switch_reconciles_invalid_values():
    controller = PlotController(
        PlotUiState(sample="s1", contig="c1", mag="m1"),
        {
            SampleScope.ONE: AvailabilitySnapshot(samples=("s1",), contigs=("c1",), mags=("m1",)),
            SampleScope.ALL: AvailabilitySnapshot(samples=("s2",), contigs=("c2",), mags=()),
        },
    )
    state = controller.switch_sample_scope(SampleScope.ALL)
    assert (state.sample, state.contig, state.mag) == ("", "", "")


def test_all_sample_selection_keeps_latest_variable():
    assert enforce_all_sample_variable(("Coverage", "Mismatch"), "Mismatch") == ("Mismatch",)
    assert enforce_all_sample_variable(("Coverage", "Mismatch")) == ("Mismatch",)
    assert enforce_all_sample_variable(()) == ()


def test_sample_scope_visibility_is_a_pure_complete_decision():
    one = visibility_for_sample_scope(SampleScope.ONE, mag_sort_category="Coverage")
    all_samples = visibility_for_sample_scope(SampleScope.ALL, mag_sort_category="Coverage")
    all_contig_sort = visibility_for_sample_scope(SampleScope.ALL, mag_sort_category="Contig")

    assert (one.sample_section, one.variables_one, one.variables_all) == (True, True, False)
    assert (all_samples.sample_section, all_samples.variables_one, all_samples.variables_all) == (False, False, True)
    assert all_samples.sample_params_header is True
    assert all_samples.sample_params_content is False
    assert all_samples.mag_sort_sample_row is True
    assert all_contig_sort.mag_sort_sample_row is False


def test_scope_projection_skips_unchanged_properties():
    class Target:
        visible = False

    first = Target()
    second = Target()
    changed = apply_changed_attributes(((first, "visible", False), (second, "visible", True)))

    assert changed == 1
    assert first.visible is False
    assert second.visible is True


def test_diagnostics_are_structured_and_generation_aware():
    output = io.StringIO()
    diagnostics = PlotDiagnostics(enabled=True, stream=output)
    diagnostics.next_generation()
    event = diagnostics.record(
        "view_change",
        time.perf_counter(),
        query_count=0,
        models=12,
        property_updates=4,
    )
    assert event.generation == 1
    assert event.query_count == 0
    assert event.models == 12
    assert event.property_updates == 4
    assert '"operation": "view_change"' in output.getvalue()


def test_diagnostics_are_disabled_without_retention_and_bounded_when_enabled():
    disabled = PlotDiagnostics()
    disabled.record("view_change", time.perf_counter())
    assert disabled.snapshot() == ()

    bounded = PlotDiagnostics(enabled=True, stream=io.StringIO(), max_events=2)
    for _ in range(3):
        bounded.record("view_change", time.perf_counter())
    assert len(bounded.snapshot()) == 2


def test_view_switch_does_not_replace_plot_and_timing_ignores_stale_acknowledgements():
    source = Path("thebigbam/plotting/application/composition_root.py").read_text(encoding="utf-8")
    scope_source = Path("thebigbam/plotting/application/scope_transition.py").read_text(encoding="utf-8")
    timing_source = Path("thebigbam/plotting/shared/timing.py").read_text(encoding="utf-8")
    assert 'views.on_change("active", scope_transition.view_changed)' in source
    assert "main_placeholder.objects" not in scope_source
    assert 'send_timing_ping("view_change")' in scope_source or "send_timing_ping('view_change')" in scope_source
    assert 'token != self.state.get("token")' in timing_source
    assert "Browser paint after patch" in timing_source
    assert "__thebigbam_view_change_started" in timing_source
