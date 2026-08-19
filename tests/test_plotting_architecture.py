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
