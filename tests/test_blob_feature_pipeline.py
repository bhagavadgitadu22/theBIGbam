import sys
import types

import duckdb
import numpy as np

if "thebigbam_rs" not in sys.modules:
    rust_stub = types.ModuleType("thebigbam_rs")
    rust_stub.decode_dense_chunk = lambda data: []
    rust_stub.decode_sparse_chunk = lambda data: ([], [])
    sys.modules["thebigbam_rs"] = rust_stub

from thebigbam.plotting.models.blobs import (
    FeatureLoadRequest,
    FeatureRegion,
    MagFeatureLoadRequest,
    MagMember,
)
from thebigbam.plotting.repositories.features import FeatureRepository
from thebigbam.plotting.services.blob_features import FeatureDataService, _slice_with_anchors


def test_feature_region_allows_display_window_before_first_base():
    region = FeatureRegion(-100, 100)

    assert region.start == -100
    assert region.length == 201


class FakeBlobRepository:
    def __init__(self):
        self.chunk_calls = 0
        self.last_contig_ids = ()

    def metadata(self, feature):
        return [("curve", "#000", 1.0, 0.2, 1, feature, "Feature_blob")]

    def variable_name(self, title, table):
        return title.lower()

    def feature_id(self, variable):
        return 1

    def storage_settings(self, variable):
        return {"scale": 1, "chunk_size": 10, "zoom_bins": (100, 1000), "window_size": 1}

    def zoom_blobs(self, contig_ids, sample_ids, feature_id, contig_level=False):
        return {}

    def chunks(self, contig_ids, sample_ids, feature_id, first, last, contig_level=False):
        self.chunk_calls += 1
        self.last_contig_ids = tuple(contig_ids)
        sample = None if contig_level else sample_ids[0]
        return {(contig_id, sample): [(0, bytes([contig_id]))] for contig_id in contig_ids}

    def contig_names(self, ids):
        return {identifier: f"contig-{identifier}" for identifier in ids}


def test_sparse_clipping_keeps_boundary_context_and_metadata_alignment():
    data = {
        "x": np.asarray([1, 5, 9]),
        "y": np.asarray([2.0, 3.0, 4.0]),
        "sparse": True,
        "sequence": ["A", "C", "G"],
    }
    result = _slice_with_anchors(data, 5, 7, "curve")
    assert result["x"] == [1, 2, 3, 5, 6, 7, 9, 10, 11]
    assert len(result["sequence"]) == len(result["x"])


def test_contig_service_decodes_and_thresholds_without_sql(monkeypatch):
    import thebigbam.plotting.services.blob_features as module

    monkeypatch.setattr(
        module,
        "decode_raw_chunks",
        lambda rows, scale, chunk_size: {
            "x": np.asarray([0, 1, 2]),
            "y": np.asarray([0.2, 0.5, 0.8]),
            "sparse": False,
        },
    )
    repository = FakeBlobRepository()
    service = FeatureDataService(repository)
    result = service.load_contig(FeatureLoadRequest("Coverage", 1, 7, FeatureRegion(1, 3), minimum_relative_value=0.5))
    assert result[0]["x"] == [1, 2, 3]
    assert result[0]["y"] == [0, 0.5, 0.8]
    assert repository.chunk_calls == 1
    assert service.decoded_points == 3


def test_mag_service_fetches_all_member_chunks_once(monkeypatch):
    import thebigbam.plotting.services.blob_features as module

    monkeypatch.setattr(
        module,
        "decode_raw_chunks",
        lambda rows, scale, chunk_size: {
            "x": np.asarray([0, 1]),
            "y": np.asarray([float(rows[0][1][0])] * 2),
            "sparse": False,
        },
    )
    repository = FakeBlobRepository()
    service = FeatureDataService(repository)
    request = MagFeatureLoadRequest(
        "Coverage",
        7,
        FeatureRegion(1, 6),
        (MagMember(1, 3, 0), MagMember(2, 3, 3)),
    )
    result = service.load_mag(request)
    assert repository.chunk_calls == 1
    assert result[0]["x"] == [1, 2, 4, 5]


def test_mag_raw_service_fetches_only_members_overlapping_region(monkeypatch):
    import thebigbam.plotting.services.blob_features as module

    monkeypatch.setattr(
        module,
        "decode_raw_chunks",
        lambda rows, scale, chunk_size: {
            "x": np.asarray([0, 1]),
            "y": np.asarray([float(rows[0][1][0])] * 2),
            "sparse": False,
        },
    )
    repository = FakeBlobRepository()
    service = FeatureDataService(repository)
    request = MagFeatureLoadRequest(
        "Coverage",
        7,
        FeatureRegion(4, 6),
        (MagMember(1, 3, 0), MagMember(2, 3, 3), MagMember(3, 3, 6)),
    )

    result = service.load_mag(request)

    assert repository.last_contig_ids == (2,)
    assert result[0]["x"] == [4, 5]
    assert service.mag_members_total == 3
    assert service.mag_members_loaded == 1


def test_repository_batches_chunks_for_samples_and_contigs_in_one_query():
    connection = duckdb.connect(":memory:")
    connection.execute(
        "CREATE TABLE Feature_blob_chunk (Contig_id INTEGER, Sample_id INTEGER, "
        "Feature_id INTEGER, Chunk_idx INTEGER, Data BLOB)"
    )
    connection.executemany(
        "INSERT INTO Feature_blob_chunk VALUES (?, ?, 3, 0, ?)",
        [(1, 10, b"a"), (1, 20, b"b"), (2, 10, b"c"), (2, 20, b"d")],
    )
    repository = FeatureRepository(connection)
    rows = repository.chunks((1, 2), (10, 20), 3, 0, 0)
    assert repository.query_count == 1
    assert set(rows) == {(1, 10), (1, 20), (2, 10), (2, 20)}
