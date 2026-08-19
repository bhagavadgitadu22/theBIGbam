"""Read-only availability queries used by plotting controls."""

from __future__ import annotations

from typing import Any

import duckdb

from ..models.filters import CompiledFilter, FilterResult


class FilteringRepository:
    def __init__(self, connection) -> None:
        self.connection = connection
        self.query_count = 0

    def _values(self, sql, parameters=()):
        self.query_count += 1
        return tuple(row[0] for row in self.connection.cursor().execute(sql, parameters).fetchall())

    @staticmethod
    def _search_clause(alias: str, column: str, search_term: str) -> tuple[str, list[str]]:
        if not search_term:
            return "", []
        return f" AND {alias}.{column} ILIKE '%' || ? || '%'", [search_term]

    def contigs_for_sample(self, sample_id: int, search_term: str = "", preserve: str = "") -> tuple[str, ...]:
        search, parameters = self._search_clause("c", "Contig_name", search_term)
        return self._values(
            "SELECT DISTINCT c.Contig_name FROM Coverage p "
            "JOIN Contig c ON c.Contig_id = p.Contig_id "
            f"WHERE p.Sample_id = ?{search} "
            "ORDER BY (c.Contig_name = ?) DESC, c.Contig_name LIMIT 100",
            [sample_id, *parameters, preserve],
        )

    def samples_for_contig(self, contig_id: int, search_term: str = "") -> tuple[str, ...]:
        search, parameters = self._search_clause("s", "Sample_name", search_term)
        return self._values(
            "SELECT DISTINCT s.Sample_name FROM Coverage p "
            "JOIN Sample s ON s.Sample_id = p.Sample_id "
            f"WHERE p.Contig_id = ?{search} ORDER BY s.Sample_name LIMIT 100",
            [contig_id, *parameters],
        )

    def samples_for_contig_unbounded(self, contig_id: int) -> tuple[str, ...]:
        """Return the semantic sample set for plotting, without the UI payload cap."""
        return self._values(
            "SELECT DISTINCT s.Sample_name FROM Coverage p "
            "JOIN Sample s ON s.Sample_id = p.Sample_id "
            "WHERE p.Contig_id = ? ORDER BY s.Sample_name",
            [contig_id],
        )

    def sample_names_for_filter(self, sql: str, parameters: list[Any]) -> tuple[str, ...]:
        """Return sample names present in a compiled filter result."""
        return self._values(
            f"SELECT DISTINCT s.Sample_name FROM ({sql}) f "
            "JOIN Sample s ON s.Sample_id = f.Sample_id ORDER BY s.Sample_name",
            parameters,
        )

    def filtered_contigs(
        self,
        sql: str,
        parameters: list[Any],
        *,
        sample_id: int | None = None,
        mag_name: str | None = None,
        search_term: str = "",
        preserve: str = "",
    ) -> tuple[str, ...]:
        joins = ""
        predicates: list[str] = []
        query_parameters = list(parameters)
        if mag_name is not None:
            joins = (
                " JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
            )
            predicates.append("mg.MAG_name = ?")
            query_parameters.append(mag_name)
        if sample_id is not None:
            predicates.append("(f.Sample_id = ? OR f.Sample_id IS NULL)")
            query_parameters.append(sample_id)
        if search_term:
            predicates.append("c.Contig_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        query_parameters.append(preserve)
        return self._values(
            f"SELECT DISTINCT c.Contig_name FROM ({sql}) f"
            f" JOIN Contig c ON c.Contig_id = f.Contig_id{joins}{where}"
            " ORDER BY (c.Contig_name = ?) DESC, c.Contig_name LIMIT 100",
            query_parameters,
        )

    def filtered_samples(
        self,
        sql: str,
        parameters: list[Any],
        *,
        contig_id: int | None = None,
        mag_name: str | None = None,
        search_term: str = "",
        preserve: tuple[str, str] = ("", ""),
    ) -> tuple[str, ...]:
        joins = ""
        predicates: list[str] = []
        query_parameters = list(parameters)
        if mag_name is not None:
            joins = (
                " JOIN MAG_contigs_association mca ON mca.Contig_id = f.Contig_id JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
            )
            predicates.append("mg.MAG_name = ?")
            query_parameters.append(mag_name)
        if contig_id is not None:
            predicates.append("f.Contig_id = ?")
            query_parameters.append(contig_id)
        if search_term:
            predicates.append("s.Sample_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        query_parameters.extend(preserve)
        return self._values(
            f"SELECT DISTINCT s.Sample_name FROM ({sql}) f"
            " JOIN Sample s ON s.Sample_id = f.Sample_id"
            f"{joins}{where}"
            " ORDER BY (s.Sample_name = ? OR s.Sample_name = ?) DESC, s.Sample_name LIMIT 100",
            query_parameters,
        )

    def filtered_mags(
        self,
        sql: str,
        parameters: list[Any],
        *,
        sample_id: int | None = None,
        search_term: str = "",
    ) -> tuple[str, ...]:
        predicates: list[str] = []
        query_parameters = list(parameters)
        if sample_id is not None:
            predicates.append("(f.Sample_id = ? OR f.Sample_id IS NULL)")
            query_parameters.append(sample_id)
        if search_term:
            predicates.append("mg.MAG_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        return self._values(
            f"SELECT DISTINCT mg.MAG_name FROM ({sql}) f"
            " JOIN MAG_contigs_association mca ON mca.Contig_id = f.Contig_id"
            f" JOIN MAG mg ON mg.MAG_id = mca.MAG_id{where} ORDER BY mg.MAG_name",
            query_parameters,
        )

    def filtered_mag_counts(self, sql: str, parameters: list[Any], mag_name: str) -> tuple[int, int]:
        self.query_count += 1
        row = (
            self.connection.cursor()
            .execute(
                f"WITH _f AS ({sql}) SELECT COUNT(DISTINCT _f.Contig_id),"
                " COUNT(DISTINCT _f.Sample_id) FROM _f"
                " JOIN MAG_contigs_association mca ON mca.Contig_id = _f.Contig_id"
                " JOIN MAG mg ON mg.MAG_id = mca.MAG_id WHERE mg.MAG_name = ?",
                [*parameters, mag_name],
            )
            .fetchone()
        )
        return int(row[0]), int(row[1])

    def mags_for_sample(self, sample_id: int, search_term: str = "") -> tuple[str, ...]:
        search, parameters = self._search_clause("mg", "MAG_name", search_term)
        return self._values(
            "SELECT DISTINCT mg.MAG_name FROM MAG mg "
            "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
            "JOIN Coverage p ON p.Contig_id = mca.Contig_id "
            f"WHERE p.Sample_id = ?{search} ORDER BY mg.MAG_name LIMIT 100",
            [sample_id, *parameters],
        )

    def count_contigs_for_sample(self, sample_id: int) -> int:
        self.query_count += 1
        return int(
            self.connection.cursor()
            .execute("SELECT COUNT(DISTINCT Contig_id) FROM Coverage WHERE Sample_id = ?", [sample_id])
            .fetchone()[0]
        )

    def count_samples_for_contig(self, contig_id: int) -> int:
        self.query_count += 1
        return int(
            self.connection.cursor()
            .execute("SELECT COUNT(DISTINCT Sample_id) FROM Coverage WHERE Contig_id = ?", [contig_id])
            .fetchone()[0]
        )

    def evaluate(self, compiled: CompiledFilter, *, has_mags: bool) -> FilterResult:
        parameters = list(compiled.parameters)
        try:
            self.query_count += 1
            pair_count, contig_count, sample_count = (
                self.connection.cursor()
                .execute(
                    f"WITH _f AS ({compiled.sql}) "
                    "SELECT COUNT(*), COUNT(DISTINCT Contig_id), COUNT(DISTINCT Sample_id) FROM _f",
                    parameters,
                )
                .fetchone()
            )
            mag_pair_count = None
            if has_mags:
                self.query_count += 1
                mag_pair_count = (
                    self.connection.cursor()
                    .execute(
                        f"WITH _f AS ({compiled.sql}) SELECT COUNT(*) FROM ("
                        "SELECT DISTINCT mg.MAG_name, _f.Sample_id FROM _f "
                        "JOIN MAG_contigs_association mca ON mca.Contig_id = _f.Contig_id "
                        "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
                        "WHERE _f.Sample_id IS NOT NULL) counted",
                        parameters,
                    )
                    .fetchone()[0]
                )
        except duckdb.Error:
            pair_count = contig_count = sample_count = 0
            mag_pair_count = 0 if has_mags else None
        return FilterResult(
            compiled,
            int(pair_count),
            int(contig_count),
            int(sample_count),
            int(mag_pair_count) if mag_pair_count is not None else None,
        )
