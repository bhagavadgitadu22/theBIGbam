from types import SimpleNamespace

from thebigbam.plotting.controls import genome as genome_controls
from thebigbam.plotting.repositories.genome_controls import GenomeControlCapabilities, GenomeControlRepository


class CapabilityConnection:
    def execute(self, sql):
        self.sql = sql
        return self

    def fetchone(self):
        if "Constants_for_plotting" in self.sql:
            return (True,)
        if "information_schema.tables" in self.sql:
            return (1,)
        return None


def test_genome_control_repository_reports_database_capabilities():
    capabilities = GenomeControlRepository(CapabilityConnection()).capabilities()
    assert capabilities == GenomeControlCapabilities(True, True, False)


def test_genome_control_factory_builds_position_controls_without_optional_tracks(monkeypatch):
    monkeypatch.setattr(
        genome_controls.ColorTemplateRepository,
        "load",
        lambda self: {},
    )
    monkeypatch.setattr(
        genome_controls.GenomeControlRepository,
        "capabilities",
        lambda self: GenomeControlCapabilities(False, False, False),
    )
    widgets = {
        "annotation_types": [],
        "custom_contig_subplots": [],
        "has_mags": False,
        "view_radio": SimpleNamespace(active=1),
        "contig_select": SimpleNamespace(value="c1"),
        "contig_lengths": {"c1": 42},
        "mag_select": SimpleNamespace(value=""),
        "mag_to_contigs": {},
        "helps_widgets": [],
    }

    controls = genome_controls.build_genome_controls(
        connection=object(),
        db_path="example.db",
        widgets=widgets,
        filtering_metadata={},
        genome_checkbox=None,
        genome_index=None,
        stylesheet="",
        toggle_stylesheet="",
        make_toggle_callback=lambda button, content: lambda event: None,
        enable_timing=False,
        interaction_lock={"locked": False},
    )

    assert controls.from_position.value == "1"
    assert controls.to_position.value == "42"
    assert controls.sequence is None
    assert controls.translated_sequence is None


def test_genome_control_projects_immutable_annotation_types_to_bokeh_list(monkeypatch):
    monkeypatch.setattr(genome_controls.ColorTemplateRepository, "load", lambda self: {})
    monkeypatch.setattr(
        genome_controls.GenomeControlRepository,
        "capabilities",
        lambda self: GenomeControlCapabilities(False, False, False),
    )
    widgets = {
        "annotation_types": ("CDS", "tRNA"),
        "custom_contig_subplots": [],
        "has_mags": False,
        "view_radio": SimpleNamespace(active=1),
        "contig_select": SimpleNamespace(value="c1"),
        "contig_lengths": {"c1": 42},
        "mag_select": SimpleNamespace(value=""),
        "mag_to_contigs": {},
        "helps_widgets": [],
    }

    controls = genome_controls.build_genome_controls(
        connection=object(),
        db_path="example.db",
        widgets=widgets,
        filtering_metadata={},
        genome_checkbox=None,
        genome_index=None,
        stylesheet="",
        toggle_stylesheet="",
        make_toggle_callback=lambda button, content: lambda event: None,
        enable_timing=False,
        interaction_lock={"locked": False},
    )

    assert controls.feature_types.options == ["CDS", "tRNA"]
    assert controls.feature_types.value == ["CDS"]
