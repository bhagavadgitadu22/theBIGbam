"""DuckDB access for the ALL SAMPLES contig plot pipeline."""

from __future__ import annotations

import re
from typing import Iterable

from ..models.plots import FeatureStyle
from .features import FeatureRepository

_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class AllSamplesRepository(FeatureRepository):
    def __init__(self, conn) -> None:
        super().__init__(conn)
        self.conn = conn
        self._source_columns: dict[str, frozenset[str]] = {}

    def _execute(self, sql, params=()):
        self.query_count += 1
        return self.conn.execute(sql, params)

    def source_columns(self, source: str) -> frozenset[str]:
        if source not in self._source_columns:
            self._source_columns[source] = frozenset(
                row[0]
                for row in self._execute(
                    "SELECT column_name FROM information_schema.columns WHERE table_name = ?",
                    [source],
                ).fetchall()
            )
        return self._source_columns[source]

    def resolve_contig(self, name: str) -> tuple[int, str, int]:
        row = self._execute(
            "SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name = ?", [name]
        ).fetchone()
        if row is None:
            raise ValueError(f"Unknown contig: {name}")
        return int(row[0]), str(row[1]), int(row[2])

    def available_samples(self, contig_id: int) -> list[tuple[int, str]]:
        rows = self._execute(
            "SELECT DISTINCT s.Sample_id, s.Sample_name FROM Coverage c "
            "JOIN Sample s ON s.Sample_id = c.Sample_id WHERE c.Contig_id = ? ORDER BY s.Sample_name",
            [contig_id],
        ).fetchall()
        return [(int(sid), str(name)) for sid, name in rows]

    def order_samples(self, contig_id, samples, ordering):
        source = ordering.source or "Sample"
        column = ordering.column
        if not _IDENTIFIER.fullmatch(source) or not _IDENTIFIER.fullmatch(column):
            raise ValueError("Unsafe sample ordering identifier")
        valid_column = self._execute(
            "SELECT 1 FROM information_schema.columns WHERE table_name = ? AND column_name = ?",
            [source, column],
        ).fetchone()
        if valid_column is None:
            raise ValueError(f"Unknown ordering column: {source}.{column}")
        ids = [sid for sid, _ in samples]
        placeholders = ",".join("?" for _ in ids)
        direction = "ASC" if ordering.ascending else "DESC"
        if source == "Sample":
            ordered = self._execute(
                f"SELECT Sample_id FROM Sample WHERE Sample_id IN ({placeholders}) "
                f'ORDER BY "{column}" {direction} NULLS LAST, Sample_name ASC',
                ids,
            ).fetchall()
            ordered_ids = [int(row[0]) for row in ordered]
        else:
            source_columns = self.source_columns(source)
            if "Contig_id" in source_columns:
                ordered = self._execute(
                    f'SELECT Sample_id FROM "{source}" WHERE Contig_id = ? '
                    f'AND Sample_id IN ({placeholders}) ORDER BY "{column}" {direction} NULLS LAST',
                    [contig_id, *ids],
                ).fetchall()
            elif "MAG_id" in source_columns:
                ordered = self._execute(
                    f'SELECT DISTINCT v.Sample_id, v."{column}" FROM "{source}" v '
                    "JOIN MAG_contigs_association mca ON mca.MAG_id = v.MAG_id "
                    f"WHERE mca.Contig_id = ? AND v.Sample_id IN ({placeholders}) "
                    f'ORDER BY v."{column}" {direction} NULLS LAST',
                    [contig_id, *ids],
                ).fetchall()
            else:
                raise ValueError(f"Ordering source has no Contig_id or MAG_id: {source}")
            ordered_ids = [int(row[0]) for row in ordered]
            ordered_ids.extend(sid for sid in ids if sid not in set(ordered_ids))
        by_id = {sid: (sid, name) for sid, name in samples}
        return [by_id[sid] for sid in ordered_ids]

    def feature_styles(self, subplot: str, feature_table: str | None = None) -> list[FeatureStyle]:
        params: list[object] = [subplot]
        where = "Subplot = ?"
        if feature_table:
            where += " AND Feature_table_name = ?"
            params.append(feature_table)
        rows = self._execute(
            'SELECT Variable_name, "Type", Color, Alpha, Fill_alpha, "Size", Title, Feature_table_name, Encoding '
            f"FROM Variable WHERE {where} ORDER BY Variable_id",
            params,
        ).fetchall()
        return [
            FeatureStyle(str(t), str(c), float(a), float(fa), float(sz), str(title), str(name), str(table), enc)
            for name, t, c, a, fa, sz, title, table, enc in rows
        ]

    def feature_id(self, variable_name: str) -> int | None:
        row = self._execute("SELECT Variable_id FROM Variable WHERE Variable_name = ?", [variable_name]).fetchone()
        return int(row[0]) if row else None

    def blob_settings(self, variable_name: str) -> tuple[float, int, list[int], int]:
        settings = self.storage_settings(variable_name)
        return settings["scale"], settings["chunk_size"], list(settings["zoom_bins"]), settings["window_size"]

    def sample_zoom_blobs(self, contig_id: int, sample_ids: Iterable[int], feature_id: int):
        rows = self.zoom_blobs((contig_id,), tuple(sample_ids), feature_id)
        return {sample_id: blob for (_contig_id, sample_id), blob in rows.items()}

    def sample_chunks(self, contig_id, sample_ids, feature_id, first_chunk, last_chunk):
        rows = self.chunks((contig_id,), tuple(sample_ids), feature_id, first_chunk, last_chunk)
        return {sample_id: chunks for (_contig_id, sample_id), chunks in rows.items()}

    def contig_zoom_blob(self, contig_id, feature_id):
        return self.zoom_blobs((contig_id,), (), feature_id, contig_level=True).get((contig_id, None))

    def contig_chunks(self, contig_id, feature_id, first_chunk, last_chunk):
        return self.chunks((contig_id,), (), feature_id, first_chunk, last_chunk, contig_level=True).get(
            (contig_id, None), []
        )
