from pathlib import Path

PLOTTING = Path("thebigbam/plotting")


def _python_sources(folder):
    return tuple((path, path.read_text(encoding="utf-8").lower()) for path in (PLOTTING / folder).rglob("*.py"))


def test_executable_sql_is_confined_to_repositories():
    for folder in ("application", "services", "renderers", "composers", "downloads"):
        for path, source in _python_sources(folder):
            assert ".execute(" not in source, f"Executable SQL leaked into {path}"


def test_repositories_and_services_do_not_import_ui_frameworks():
    for folder in ("repositories", "services"):
        for path, source in _python_sources(folder):
            assert "import bokeh" not in source and "from bokeh" not in source, f"Bokeh leaked into {path}"
            assert "import panel" not in source and "from panel" not in source, f"Panel leaked into {path}"


def test_database_access_is_confined_to_repositories_and_application_connection_setup():
    forbidden = (".execute(", ".cursor(", "database_getters", "import duckdb", "from duckdb")
    for folder in ("controls", "services", "renderers", "composers", "reports", "settings", "downloads"):
        for path, source in _python_sources(folder):
            for marker in forbidden:
                assert marker not in source, f"Database access marker {marker!r} leaked into {path}"


def test_presentation_layers_do_not_depend_on_repositories():
    for folder in ("controls", "renderers", "composers", "reports", "settings", "downloads"):
        for path, source in _python_sources(folder):
            assert "from ..repositories" not in source, f"Repository dependency leaked into {path}"
            assert "from ...repositories" not in source, f"Repository dependency leaked into {path}"


def test_composers_only_assemble_already_rendered_models():
    sources = _python_sources("composers")
    assert {path.name for path, _source in sources} == {"__init__.py", "layout.py"}
    for path, source in sources:
        for marker in ("connection", "repository", "service", "duckdb", "request"):
            assert marker not in source, f"Pipeline responsibility {marker!r} leaked into {path}"


def test_application_owns_plot_pipeline_orchestration():
    application = PLOTTING / "application"
    for filename in ("sample_mag_pipeline.py", "all_samples_pipeline.py", "genome_track_pipeline.py"):
        assert (application / filename).is_file()
    assert "build_mag_plot" in (application / "sample_mag_pipeline.py").read_text(encoding="utf-8")
    assert "assemble_grid" in (application / "all_samples_pipeline.py").read_text(encoding="utf-8")


def test_four_apply_handlers_have_independent_typed_modules():
    handler_dir = PLOTTING / "application" / "apply_handlers"
    expected = {
        "contig_one.py": "ContigOneApplyRequest",
        "contig_all.py": "ContigAllApplyRequest",
        "mag_one.py": "MagOneApplyRequest",
        "mag_all.py": "MagAllApplyRequest",
    }
    for filename, request_type in expected.items():
        source = (handler_dir / filename).read_text(encoding="utf-8")
        assert request_type in source
        assert "is_all" not in source
        assert ".execute(" not in source

    assert len((handler_dir / "rendering.py").read_text(encoding="utf-8").splitlines()) <= 550
