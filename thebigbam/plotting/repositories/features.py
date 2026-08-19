"""DuckDB access for feature metadata used by plot composers."""

from __future__ import annotations


class FeatureRepository:
    def __init__(self, connection_or_cursor) -> None:
        self.database = connection_or_cursor
        self.query_count = 0

    def _execute(self, sql, params=()):
        self.query_count += 1
        return self.database.execute(sql, params)

    def metadata(self, subplot):
        return self._execute(
            'SELECT "Type", Color, Alpha, Fill_alpha, "Size", Title, Feature_table_name FROM Variable WHERE Subplot=?',
            (subplot,),
        ).fetchall()

    def metadata_batch(self, subplots):
        if not subplots:
            return {}
        placeholders = ", ".join("?" for _ in subplots)
        rows = self._execute(
            f'SELECT Subplot, "Type", Color, Alpha, Fill_alpha, "Size", Title, Feature_table_name '
            f"FROM Variable WHERE Subplot IN ({placeholders}) ORDER BY Module_order",
            tuple(subplots),
        ).fetchall()
        result = {}
        for row in rows:
            result.setdefault(row[0], []).append(row[1:])
        return result

    def variable_name(self, title, table):
        row = self._execute(
            "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name=?",
            (title, table),
        ).fetchone()
        return str(row[0]) if row else None

    def feature_id(self, variable_name):
        row = self._execute("SELECT Variable_id FROM Variable WHERE Variable_name=?", (variable_name,)).fetchone()
        return int(row[0]) if row else None

    def storage_settings(self, variable_name):
        scale = self._execute(
            "SELECT Scale FROM Column_scales WHERE Feature_name=? AND Column_name='Value'",
            (variable_name,),
        ).fetchone()
        rows = self._execute(
            "SELECT Key, Value FROM Database_metadata WHERE Key IN "
            "('Chunk_size','Zoom_bin_sizes','GC_content_window_size','GC_skew_window_size')"
        ).fetchall()
        values = {str(key): value for key, value in rows}
        gc_key = "GC_content_window_size" if variable_name == "gc_content" else "GC_skew_window_size"
        gc_window = int(values.get(gc_key, 1)) if variable_name in {"gc_content", "gc_skew"} else 1
        return {
            "scale": int(scale[0]) if scale else 1,
            "chunk_size": int(values.get("Chunk_size", 65536)),
            "zoom_bins": tuple(int(value) for value in str(values.get("Zoom_bin_sizes", "100,1000,10000")).split(",")),
            "window_size": gc_window,
        }

    def zoom_blobs(self, contig_ids, sample_ids, feature_id, contig_level=False):
        if not contig_ids:
            return {}
        contig_marks = ",".join("?" for _ in contig_ids)
        if contig_level:
            rows = self._execute(
                f"SELECT Contig_id, NULL, Zoom_data FROM Contig_blob WHERE Feature_id=? "
                f"AND Contig_id IN ({contig_marks})",
                [feature_id, *contig_ids],
            ).fetchall()
        else:
            if not sample_ids:
                return {}
            sample_marks = ",".join("?" for _ in sample_ids)
            rows = self._execute(
                f"SELECT Contig_id, Sample_id, Zoom_data FROM Feature_blob WHERE Feature_id=? "
                f"AND Contig_id IN ({contig_marks}) AND Sample_id IN ({sample_marks})",
                [feature_id, *contig_ids, *sample_ids],
            ).fetchall()
        return {(int(cid), int(sid) if sid is not None else None): blob for cid, sid, blob in rows}

    def chunks(self, contig_ids, sample_ids, feature_id, first_chunk, last_chunk, contig_level=False):
        if not contig_ids:
            return {}
        contig_marks = ",".join("?" for _ in contig_ids)
        if contig_level:
            rows = self._execute(
                f"SELECT Contig_id, NULL, Chunk_idx, Data FROM Contig_blob_chunk WHERE Feature_id=? "
                f"AND Contig_id IN ({contig_marks}) AND Chunk_idx BETWEEN ? AND ? ORDER BY Contig_id, Chunk_idx",
                [feature_id, *contig_ids, first_chunk, last_chunk],
            ).fetchall()
        else:
            if not sample_ids:
                return {}
            sample_marks = ",".join("?" for _ in sample_ids)
            rows = self._execute(
                f"SELECT Contig_id, Sample_id, Chunk_idx, Data FROM Feature_blob_chunk WHERE Feature_id=? "
                f"AND Contig_id IN ({contig_marks}) AND Sample_id IN ({sample_marks}) "
                "AND Chunk_idx BETWEEN ? AND ? ORDER BY Contig_id, Sample_id, Chunk_idx",
                [feature_id, *contig_ids, *sample_ids, first_chunk, last_chunk],
            ).fetchall()
        result = {}
        for cid, sid, chunk_idx, data in rows:
            result.setdefault((int(cid), int(sid) if sid is not None else None), []).append((int(chunk_idx), data))
        return result

    def contig_names(self, contig_ids):
        if not contig_ids:
            return {}
        marks = ",".join("?" for _ in contig_ids)
        rows = self._execute(
            f"SELECT Contig_id, Contig_name FROM Contig WHERE Contig_id IN ({marks})", list(contig_ids)
        ).fetchall()
        return {int(cid): str(name) for cid, name in rows}

    def contig_length(self, contig_id):
        row = self._execute("SELECT Contig_length FROM Contig WHERE Contig_id=?", (contig_id,)).fetchone()
        return int(row[0]) if row else 1


def split_features_by_storage(metadata, requested_features):
    """Return contig-only and sample-dependent subplot names in request order."""
    contig_features, sample_features = [], []
    for feature in requested_features:
        rows = metadata.get(feature)
        target = contig_features if rows and all(row[6] == "Contig_blob" for row in rows) else sample_features
        target.append(feature)
    return contig_features, sample_features
