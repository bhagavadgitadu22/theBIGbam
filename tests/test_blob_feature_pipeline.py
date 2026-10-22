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


def test_mag_merge_preserves_partner_metadata_when_another_member_has_no_hits():
    from thebigbam.plotting.services.blob_features import _merge

    merged = _merge(
        [
            {"x": np.array([0, 9]), "y": np.array([0, 0]), "sparse": False},
            {"x": np.array([11]), "y": np.array([98.5]), "sparse": True, "partner_contig_id": np.array([7])},
        ]
    )
    assert merged["partner_contig_id"].tolist() == [None, None, 7]


def test_missing_identity_scale_uses_variable_metadata():
    conn = duckdb.connect()
    conn.execute("CREATE TABLE Column_scales (Feature_name TEXT, Column_name TEXT, Scale INTEGER)")
    conn.execute("CREATE TABLE Variable (Variable_name TEXT, Value_scale TEXT)")
    conn.execute("INSERT INTO Variable VALUES ('hit_identity', 'times100'), ('hit_identity_within_mag', 'times100')")
    conn.execute("CREATE TABLE Database_metadata (Key TEXT, Value TEXT)")
    conn.execute("INSERT INTO Database_metadata VALUES ('Blast_hit_storage','query-v1')")
    repository = FeatureRepository(conn)
    assert repository.storage_settings("hit_identity")["scale"] == 100
    assert repository.storage_settings("hit_identity_within_mag")["scale"] == 100


def test_zoom_best_partner_respects_mag_scope_offsets_and_identity():
    from thebigbam.plotting.services.blob_features import _format_zoom

    conn = duckdb.connect()
    conn.execute("CREATE TABLE Database_metadata (Key TEXT, Value TEXT)")
    conn.execute("INSERT INTO Database_metadata VALUES ('Blast_hit_storage','query-v1')")
    conn.execute("INSERT INTO Database_metadata VALUES ('View_mode', 'mag')")
    conn.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_name TEXT)")
    conn.execute("INSERT INTO Contig VALUES (1,'query'),(2,'within'),(3,'outside'),(4,'lower')")
    conn.execute("CREATE TABLE MAG_contigs_association (MAG_id INTEGER, Contig_id INTEGER)")
    conn.execute("INSERT INTO MAG_contigs_association VALUES (1,1),(1,2),(2,3),(2,4)")
    conn.execute(
        "CREATE TABLE Contig_blast_hits (Contig_id INTEGER, Contig_id_other INTEGER, Position1 INTEGER, Position2 INTEGER, Pident INTEGER)"
    )
    conn.execute(
        "INSERT INTO Contig_blast_hits VALUES (1,2,1,10,9900),(1,3,10,1,10000),(1,4,1,10,9700),(1,1,1,10,10000)"
    )
    conn.execute("CREATE VIEW Contig_blast_hits_query AS SELECT * FROM Contig_blast_hits")
    service = FeatureDataService(FeatureRepository(conn))

    def values(identity):
        return _format_zoom(
            {"bin_start": np.array([100]), "bin_end": np.array([109]), "mean": np.array([identity]), "sparse": True},
            "curve",
            101,
            111,
        )

    outside = values(100.0)
    service._resolve_partners(outside, "hit_identity", ((1, 100, 10),))
    assert outside["partner_contig"] == [None, "outside", None]
    assert outside["y"] == [0.0, 100.0, 0.0]
    within = values(99.0)
    service._resolve_partners(within, "hit_identity_within_mag", ((1, 100, 10),))
    assert within["partner_contig"] == [None, "within", None]
    conn.execute("UPDATE Database_metadata SET Value='contig'")
    assert service.repository.best_hit_names([(0, 1, 1, 10)], False) == {0: "outside"}


def test_identity_zoom_decodes_stored_10000_as_100_percent():
    import struct

    import zstandard

    # One sparse zoom bin: delta index 0 and zigzag identity 10000.
    identity = bytes([0xA0, 0x9C, 0x01])
    level = struct.pack("<II", 1, 1) + b"\x00" + struct.pack("<I", len(identity)) + identity
    compressed = zstandard.ZstdCompressor().compress(level)
    blob = b"TBZ\x01" + bytes([1, 1, 1, 0]) + struct.pack("<II", 10000, 0)
    blob += struct.pack("<I", len(compressed)) + compressed
    values = FeatureDataService._decode_zoom(blob, {"scale": 100, "zoom_bins": (10000,)}, "curve", 1, 10001)
    assert max(values["y"]) == 100.0


def test_repeat_intervals_preserve_second_plateau_next_to_empty_mag_member(monkeypatch):
    from thebigbam.plotting.services.blob_features import _merge

    monkeypatch.setattr(
        "thebigbam.plotting.services.blob_features.decode_raw_sparse_chunks",
        lambda *_: {
            "x": np.array([0, 38719, 38720, 82513, 82514, 119994]),
            "y": np.array([88.49, 88.49, 0, 0, 88.49, 88.49]),
            "sparse": True,
            "partner_contig_id": np.array([82514, 119994, -1, -1, 0, 38719]),
        },
    )
    decoded = FeatureDataService._raw_chunks(
        [(0, b"chunk")], {"scale": 100, "chunk_size": 65536, "intervals": True}, True
    )
    assert decoded["repeat_positions"][-2:] == ["1-38720", "1-38720"]
    empty = {"x": np.array([119995, 125422]), "y": np.array([0.0, 0.0]), "sparse": False}
    values = _slice_with_anchors(_merge([decoded, empty]), 1, 125423, "curve")
    second = [(x, y) for x, y in zip(values["x"], values["y"]) if 82515 <= x <= 119995]
    assert second == [(82515, 88.49), (119995, 88.49)]
    # A window entirely inside the second repeat retains both plateau endpoints.
    inside = _slice_with_anchors(decoded, 100000, 105000, "curve")
    assert inside["x"] == [100000, 105000]
    assert inside["y"] == [88.49, 88.49]


def test_repeat_window_clipping_preserves_plateau_and_metadata_at_both_edges():
    data = {
        "x": np.array([9, 10, 99, 100]),
        "y": np.array([0.0, 88.49, 88.49, 0.0]),
        "sparse": False,
        "intervals": True,
        "repeat_positions": [None, "201-290", "201-290", None],
    }
    for start, end, expected_x, expected_y in [
        (50, 110, [50, 100, 101], [88.49, 88.49, 0.0]),
        (1, 50, [10, 11, 50], [0.0, 88.49, 88.49]),
        (50, 75, [50, 75], [88.49, 88.49]),
    ]:
        values = _slice_with_anchors(data, start, end, "curve")
        assert values["x"] == expected_x
        assert values["y"] == expected_y
        for value, label in zip(values["y"], values["repeat_positions"]):
            assert label == ("201-290" if value else None)


def test_repeat_zoom_clips_plateau_to_window_without_end_peak():
    from thebigbam.plotting.services.blob_features import _format_zoom

    data = {
        "bin_start": np.array([0, 100, 200]),
        "bin_end": np.array([99, 199, 299]),
        "mean": np.array([88.49, 88.49, 88.49]),
        "sparse": True,
    }
    values = _format_zoom(data, "curve", 150, 250, intervals=True)
    assert values["x"][0] == 150 and values["x"][-1] == 250
    assert all(value == 88.49 for value in values["y"])


def test_repeat_chunk_and_mag_merging_keep_interval_clipping(monkeypatch):
    from thebigbam.plotting.services.blob_features import _merge

    # A fetched interior chunk carries only its two plateau endpoints.
    monkeypatch.setattr(
        "thebigbam.plotting.services.blob_features.decode_raw_sparse_chunks",
        lambda *_: {"x": np.array([65536, 131071]), "y": np.array([1.0, 1.0]), "sparse": True},
    )
    decoded = FeatureDataService._raw_chunks(
        [(1, b"chunk")], {"scale": 1, "chunk_size": 65536, "intervals": True}, True
    )
    decoded["x"] += 1000
    empty = {"x": np.array([201000, 209999]), "y": np.array([0.0, 0.0]), "sparse": False}
    values = _slice_with_anchors(_merge([decoded, empty]), 100000, 105000, "curve")
    assert values["x"] == [100000, 105000]
    assert values["y"] == [1.0, 1.0]
