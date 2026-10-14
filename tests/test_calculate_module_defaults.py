import argparse

import duckdb
import pytest

from thebigbam.database import calculating_data


@pytest.mark.parametrize(
    "view,explicit,expected",
    [
        ("mag", None, ["Coverage", "Misalignment"]),
        ("contig", None, ["Coverage", "Misalignment"]),
        ("contig", "coverage,rna", ["Coverage", "RNA"]),
        ("mag", "coverage,rna", ["Coverage", "RNA"]),
        ("mag", "coverage,phagetermini", ["Coverage", "Phage termini"]),
        ("mag", "phagetermini", ["Phage termini"]),
    ],
)
def test_calculate_module_selection(tmp_path, monkeypatch, view, explicit, expected):
    bam = tmp_path / "sample.bam"
    bam.touch()
    monkeypatch.setattr(calculating_data, "resolve_mag_inputs", lambda *args: ([], [], []))
    monkeypatch.setattr(calculating_data, "validate_unique_contig_names", lambda *args: None)
    calls = []
    monkeypatch.setattr(calculating_data, "calculating_all_features_parallel",
                        lambda modules, *args, **kwargs: calls.append(modules))
    parser = argparse.ArgumentParser()
    calculating_data.add_calculate_args(parser)
    argv = ["-b", str(bam), "-o", str(tmp_path / "output.db"), "--view", view]
    if explicit is not None:
        argv += ["-m", explicit]
    calculating_data.run_calculate_args(parser.parse_args(argv))
    assert calls == [expected]


@pytest.mark.parametrize("recorded,expected", [
    ("Coverage,Misalignment", ["Coverage", "Misalignment"]),
    ("coverage,phagetermini", ["Coverage", "Phage termini"]),
    (None, sorted(set(calculating_data.MODULE_ALIASES.values()))),
])
def test_extend_uses_recorded_modules_before_variable_cleanup(tmp_path, monkeypatch, recorded, expected):
    db = tmp_path / "interrupted.db"
    bam = tmp_path / "new.bam"
    bam.touch()
    with duckdb.connect(str(db)) as conn:
        conn.execute("CREATE TABLE Variable (Module VARCHAR)")
        conn.executemany("INSERT INTO Variable VALUES (?)",
                         [(module,) for module in set(calculating_data.MODULE_ALIASES.values())])
        conn.execute("CREATE TABLE Database_metadata (Key VARCHAR, Value VARCHAR)")
        if recorded is not None:
            conn.execute("INSERT INTO Database_metadata VALUES ('Modules', ?)", [recorded])
    monkeypatch.setattr(calculating_data, "resolve_mag_inputs", lambda *args: ([], [], []))
    monkeypatch.setattr(calculating_data, "validate_unique_contig_names", lambda *args: None)
    calls = []
    monkeypatch.setattr(calculating_data, "calculating_all_features_parallel",
                        lambda modules, *args, **kwargs: calls.append(modules))
    parser = argparse.ArgumentParser()
    calculating_data.add_calculate_args(parser)
    args = parser.parse_args(["--extend", "-b", str(bam), "-o", str(db), "--view", "mag"])
    calculating_data.run_calculate_args(args)
    assert calls == [expected]
