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
        if scale is None:
            declared = self._execute("SELECT Value_scale FROM Variable WHERE Variable_name=?", (variable_name,)).fetchone()
            scale = ({"raw": 1, "times100": 100, "times1000": 1000, "times10": 10}.get(str(declared[0]), 1),) if declared else (1,)
        rows = self._execute(
            "SELECT Key, Value FROM Database_metadata WHERE Key IN "
            "('Chunk_size','Zoom_bin_sizes','GC_content_window_size','GC_skew_window_size','Repeat_blob_encoding')"
        ).fetchall()
        values = {str(key): value for key, value in rows}
        if variable_name in {"hit_count", "hit_count_within_mag", "hit_identity", "hit_identity_within_mag"}:
            marker = self._execute("SELECT Value FROM Database_metadata WHERE Key='Blast_hit_storage'").fetchone()
            if not marker or marker[0] != "query-v1":
                raise ValueError("Legacy BLAST hit tracks are unsupported. Rebuild the database for query-only hits.")
        if variable_name in {"direct_repeat_count", "inverted_repeat_count", "direct_repeat_identity", "inverted_repeat_identity"} and values.get("Repeat_blob_encoding") not in {"interval-v1", "interval-v2"}:
            raise ValueError("Repeat blobs use the old endpoint encoding. Regenerate the database to plot complete repeat intervals.")
        gc_key = "GC_content_window_size" if variable_name == "gc_content" else "GC_skew_window_size"
        gc_window = int(values.get(gc_key, 1)) if variable_name in {"gc_content", "gc_skew"} else 1
        return {
            "scale": int(scale[0]) if scale else 1,
            "chunk_size": int(values.get("Chunk_size", 65536)),
            "zoom_bins": tuple(int(value) for value in str(values.get("Zoom_bin_sizes", "100,1000,10000")).split(",")),
            "window_size": gc_window,
            "intervals": variable_name in {"direct_repeat_count", "inverted_repeat_count", "direct_repeat_identity", "inverted_repeat_identity"} and values.get("Repeat_blob_encoding") in {"interval-v1", "interval-v2"},
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

    def best_hit_names(self, targets, within_mag):
        """Resolve the maximum-identity partner per visible bin, without loading DNA."""
        if not targets:
            return {}
        exists = self._execute("SELECT EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name='Contig_blast_hits')").fetchone()[0]
        if not exists:
            return {}
        mode = self._execute("SELECT Value FROM Database_metadata WHERE Key='View_mode'").fetchone()
        predicate = "TRUE"
        if mode and mode[0] == "mag":
            same_mag = "EXISTS (SELECT 1 FROM MAG_contigs_association a JOIN MAG_contigs_association b ON a.MAG_id=b.MAG_id WHERE a.Contig_id=t.cid AND b.Contig_id=h.Contig_id_other)"
            predicate = same_mag if within_mag else f"NOT ({same_mag})"
        columns = [list(column) for column in zip(*targets)]
        rows = self._execute(f"""
            WITH targets AS (SELECT unnest(?::BIGINT[]) AS idx, unnest(?::BIGINT[]) AS cid,
                unnest(?::BIGINT[]) AS lo, unnest(?::BIGINT[]) AS hi)
            SELECT t.idx, best.Contig_name FROM targets t
            LEFT JOIN LATERAL (
                SELECT c.Contig_name FROM Contig_blast_hits_query h
                JOIN Contig c ON c.Contig_id=h.Contig_id_other
                WHERE h.Contig_id=t.cid AND h.Contig_id_other<>t.cid
                  AND LEAST(h.Position1,h.Position2)<=t.hi
                  AND GREATEST(h.Position1,h.Position2)>=t.lo AND {predicate}
                ORDER BY h.Pident DESC, h.Contig_id_other, h.Position1, h.Position2 LIMIT 1
            ) best ON TRUE
        """, columns).fetchall()
        return {int(index): name for index, name in rows}

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
