import sys
from types import SimpleNamespace

import pytest

from thebigbam.plotting.repositories.exports import (
    download_contig_metrics_csv,
    download_mag_metrics_csv,
)


class _FailingConnection:
    def __init__(self):
        self.closed = False

    def execute(self, _query):
        raise RuntimeError("query failed")

    def close(self):
        self.closed = True


@pytest.mark.parametrize(
    ("export", "subject"),
    [
        (download_contig_metrics_csv, "contig-1"),
        (download_mag_metrics_csv, "mag-1"),
    ],
)
def test_failed_export_always_closes_duckdb_connection(monkeypatch, export, subject):
    connection = _FailingConnection()
    monkeypatch.setitem(
        sys.modules,
        "duckdb",
        SimpleNamespace(connect=lambda *_args, **_kwargs: connection),
    )

    assert export("example.db", subject, ["sample-1"]) is None
    assert connection.closed is True
