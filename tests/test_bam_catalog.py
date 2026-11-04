"""BAM discovery and warning behavior for serve read plotting."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from thebigbam.plotting.application import server
from thebigbam.plotting.services.bam_catalog import build_bam_catalog

FIXTURE = Path(__file__).parent / "HK97" / "HK97_illumina_circular.bam"
CONTIG = "NC_002167.1"
LENGTH = 39732


def _copy_bam(directory: Path, name: str, *, index: bool = True) -> Path:
    path = directory / f"{name}.bam"
    shutil.copyfile(FIXTURE, path)
    if index:
        shutil.copyfile(Path(f"{FIXTURE}.bai"), Path(f"{path}.bai"))
    return path


def test_single_bam_matches_calculate_sample_name_and_reference(tmp_path):
    path = _copy_bam(tmp_path, "sample_with_MD")

    catalog = build_bam_catalog(path, ("sample",), {CONTIG: LENGTH})

    assert catalog.warnings == ()
    assert catalog.by_sample["sample"].path == path
    assert catalog.by_sample["sample"].available
    assert catalog.warning_for_sample("sample", CONTIG, LENGTH) is None
    assert "length 39,732 bp" in catalog.warning_for_sample("sample", CONTIG, LENGTH + 1)


def test_directory_reports_missing_samples_and_unmatched_bams(tmp_path):
    _copy_bam(tmp_path, "sample")
    _copy_bam(tmp_path, "extra")

    catalog = build_bam_catalog(tmp_path, ("sample", "missing"), {CONTIG: LENGTH})

    assert sorted(catalog.by_sample) == ["extra", "sample"]
    assert any("no matching database sample 'extra'" in warning for warning in catalog.warnings)
    assert any("No BAM file matches database sample 'missing'" in warning for warning in catalog.warnings)
    assert "No BAM file matches" in catalog.warning_for_sample("missing", CONTIG, LENGTH)


def test_missing_index_and_reference_mismatch_are_sample_warnings(tmp_path):
    path = _copy_bam(tmp_path, "sample", index=False)
    catalog = build_bam_catalog(path, ("sample",), {CONTIG: LENGTH + 1})

    assert not catalog.by_sample["sample"].available
    assert any("no .bai or .csi index" in warning for warning in catalog.warnings)
    assert any("reference length mismatch" in warning for warning in catalog.warnings)


def test_duplicate_normalized_names_are_rejected(tmp_path):
    _copy_bam(tmp_path, "sample")
    _copy_bam(tmp_path, "sample_with_MD")

    with pytest.raises(ValueError, match="Multiple BAM files map to sample 'sample'"):
        build_bam_catalog(tmp_path, ("sample",), {CONTIG: LENGTH})


def test_serve_passes_catalog_into_each_session(monkeypatch, tmp_path):
    path = _copy_bam(tmp_path, "sample")
    preloaded = SimpleNamespace(samples=("sample",), contig_lengths={CONTIG: LENGTH})
    serve_args = {}
    layout_args = {}
    monkeypatch.setattr(server, "_print_database_metadata", lambda _path: None)
    monkeypatch.setattr(server.application, "preload_db_data", lambda *_args, **_kwargs: preloaded)
    monkeypatch.setattr(server, "warm_plot_pipeline_imports", lambda: None)
    monkeypatch.setattr(server.application, "create_layout", lambda *_args, **kwargs: layout_args.update(kwargs))
    monkeypatch.setattr(server.pn, "serve", lambda app, **_kwargs: serve_args.update(app=app))
    parser = argparse.ArgumentParser()
    server.add_serve_args(parser)

    assert server.run_serve(parser.parse_args(["--db", "example.db", "-b", str(path)])) == 0
    serve_args["app"]()

    assert layout_args["bam_catalog"].by_sample["sample"].path == path
