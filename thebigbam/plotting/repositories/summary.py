"""DuckDB retrieval for summary reports."""

from __future__ import annotations

from typing import Any, Iterable


class SummaryRepository:
    def __init__(self, connection: Any) -> None:
        self.connection = connection
        self.query_count = 0

    def _execute(self, sql: str, parameters: Iterable[Any] = ()):
        self.query_count += 1
        return self.connection.execute(sql, list(parameters))

    def column_scales(self):
        try:
            return self._execute(
                "SELECT Column_name, Scale FROM Column_scales "
                "WHERE Feature_name IN ('Contig', 'Coverage', 'Side_misassembly')"
            ).fetchall()
        except Exception:
            return []

    def named_row(self, table: str, name_column: str, name: str):
        cursor = self._execute(f'SELECT * FROM "{table}" WHERE "{name_column}" = ?', [name])
        return cursor.fetchone(), tuple(item[0] for item in cursor.description)

    def named_rows(self, table: str, name_column: str, names):
        names = tuple(names)
        if not names:
            return [], ()
        placeholders = ",".join("?" for _ in names)
        cursor = self._execute(
            f'SELECT * FROM "{table}" WHERE "{name_column}" IN ({placeholders}) ORDER BY "{name_column}"',
            names,
        )
        return cursor.fetchall(), tuple(item[0] for item in cursor.description)

    def identifier(self, table: str, id_column: str, name_column: str, name: str):
        row = self._execute(f'SELECT "{id_column}" FROM "{table}" WHERE "{name_column}" = ?', [name]).fetchone()
        return row[0] if row else None

    def sample_ids(self, names) -> dict[str, int]:
        names = tuple(names)
        if not names:
            return {}
        placeholders = ",".join("?" for _ in names)
        return dict(
            self._execute(
                f"SELECT Sample_name, Sample_id FROM Sample WHERE Sample_name IN ({placeholders})", names
            ).fetchall()
        )

    def view_exists(self, view: str) -> bool:
        try:
            self._execute(f'SELECT 1 FROM "{view}" LIMIT 0')
            return True
        except Exception:
            return False

    def metric_rows(self, view: str, subject_column: str, subject_id: int, sample_ids, columns):
        sample_ids = tuple(sample_ids)
        if not sample_ids or not columns:
            return []
        placeholders = ",".join("?" for _ in sample_ids)
        return self._execute(
            f'SELECT Sample_name, {", ".join(columns)} FROM "{view}" '
            f'WHERE "{subject_column}" = ? AND Sample_id IN ({placeholders})',
            (subject_id, *sample_ids),
        ).fetchall()

    def has_paired_samples(self, names) -> bool:
        names = tuple(names)
        if not names:
            return False
        placeholders = ",".join("?" for _ in names)
        row = self._execute(
            f"SELECT COUNT(*) FROM Sample WHERE Sample_name IN ({placeholders}) AND Sequencing_type = 'paired-short'",
            names,
        ).fetchone()
        return bool(row and row[0])

    def mag_for_contig(self, contig_name: str):
        try:
            cursor = self._execute(
                "SELECT mg.* FROM MAG mg JOIN MAG_contigs_association mca ON mca.MAG_id=mg.MAG_id "
                "JOIN Contig c ON c.Contig_id=mca.Contig_id WHERE c.Contig_name=?",
                [contig_name],
            )
            return cursor.fetchone(), tuple(item[0] for item in cursor.description)
        except Exception:
            return None, ()

    def contigs_for_mag(self, mag_name: str):
        try:
            cursor = self._execute(
                "SELECT c.* FROM Contig c JOIN MAG_contigs_association mca ON mca.Contig_id=c.Contig_id "
                "JOIN MAG mg ON mg.MAG_id=mca.MAG_id WHERE mg.MAG_name=? ORDER BY c.Contig_name",
                [mag_name],
            )
            return cursor.fetchall(), tuple(item[0] for item in cursor.description)
        except Exception:
            return [], ()
