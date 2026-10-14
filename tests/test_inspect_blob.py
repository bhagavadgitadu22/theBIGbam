import duckdb
import numpy as np
import pytest

from thebigbam.database import inspect_blob


@pytest.mark.parametrize("feature", [
    "hit_count_within_mag", "hit_count",
    "hit_identity_within_mag", "hit_identity", "gc_content",
])
def test_contig_inspection_selects_blob_encoding_and_preserves_mag_coordinates(monkeypatch, feature):
    conn = duckdb.connect()
    conn.execute("CREATE TABLE Database_metadata (Key TEXT, Value TEXT)")
    conn.execute("INSERT INTO Database_metadata VALUES ('Blast_hit_storage','query-v1')")
    sparse = feature != "gc_content"
    scale = 100 if "identity" in feature else 1
    zoom = b"TBZ\x01" + bytes([1 if scale == 100 else 0, int(sparse), 3, 0]) + bytes(8)
    conn.execute("CREATE TABLE Contig_blob (Contig_id INTEGER, Feature_id INTEGER, Zoom_data BLOB)")
    conn.execute("CREATE TABLE Contig_blob_chunk (Contig_id INTEGER, Feature_id INTEGER, Chunk_idx INTEGER, Data BLOB)")
    conn.execute("INSERT INTO Contig_blob VALUES (1,1,?)", [zoom])
    conn.execute("INSERT INTO Contig_blob_chunk VALUES (1,1,0,?)", [b"chunk"])
    monkeypatch.setattr(inspect_blob, "get_blob_scale", lambda *_: scale)
    monkeypatch.setattr(inspect_blob, "get_chunk_size", lambda *_: 65536)
    monkeypatch.setattr(inspect_blob, "get_gc_window_size", lambda *_: 1 if sparse else 500)
    def wrong_decoder(*_):
        pytest.fail("Selected the wrong decoder for the blob encoding")
    def correct_decoder(rows, divisor, chunk_size):
        assert rows == [(0, b"chunk")]
        assert divisor == scale
        return {"x": np.array([0, 1]), "y": np.array([10000, 9900]) / divisor}
    monkeypatch.setattr(inspect_blob, "decode_raw_sparse_chunks", correct_decoder if sparse else wrong_decoder)
    monkeypatch.setattr(inspect_blob, "decode_raw_chunks", wrong_decoder if sparse else correct_decoder)
    captured = []
    monkeypatch.setattr(inspect_blob, "_print_rows", lambda *args, **kwargs: captured.append((args, kwargs)))
    inspect_blob._inspect_contig_feature(conn, 1, "contig", feature, 1, None, None, None, mag_name="mag", mag_offset=200)
    args, kwargs = captured[0]
    assert args[3].tolist() == ([0, 1] if sparse else [0, 500])
    assert args[4].tolist() == ([100.0, 99.0] if scale == 100 else [10000.0, 9900.0])
    assert kwargs["mag_offset"] == 200
    assert kwargs["mag_name"] == "mag"
    conn.close()


@pytest.mark.parametrize("feature", ["direct_repeat_count", "inverted_repeat_count", "direct_repeat_identity", "inverted_repeat_identity"])
def test_repeat_inspection_exports_full_intervals_and_clips_regions(monkeypatch, feature):
    conn = duckdb.connect()
    table = "Contig_directRepeats" if feature.startswith("direct_") else "Contig_invertedRepeats"
    conn.execute(f"CREATE TABLE {table} (Contig_id INTEGER, Position1 INTEGER, Position2 INTEGER, Pident INTEGER)")
    conn.execute(f"INSERT INTO {table} VALUES (1,1,38720,884900),(1,119994,82515,884900)")
    conn.execute("CREATE TABLE Column_scales (Feature_name TEXT, Column_name TEXT, Scale INTEGER)")
    conn.execute("INSERT INTO Column_scales VALUES (?, 'Pident', 10000)", [table])
    captured = []
    monkeypatch.setattr(inspect_blob, "_write_rle", lambda *args, **kwargs: captured.append((args, kwargs)))
    inspect_blob._inspect_contig_feature(conn, 1, "contig", feature, 1, 100, 100000, None, mag_name="mag", mag_offset=200)
    args, kwargs = captured[0]
    assert args[3].tolist() == [300, 82715]
    assert args[4].tolist() == [38920, 100200]
    assert args[5].tolist() == ([88.49, 88.49] if feature.endswith("identity") else [1, 1])
    assert kwargs["mag_name"] == "mag"
    conn.close()
