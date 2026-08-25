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
        return tuple(row[0] for row in self.connection.execute(sql, parameters).fetchall())

    @staticmethod
    def _search_clause(alias: str, column: str, search_term: str) -> tuple[str, list[str]]:
        if not search_term:
            return "", []
        return f" AND {alias}.{column} ILIKE '%' || ? || '%'", [search_term]

    @staticmethod
    def _filtered_contig_scope(parameters, *, sample_id=None, mag_name=None, search_term=""):
        joins = ""
        predicates: list[str] = []
        query_parameters = list(parameters)
        if mag_name is not None:
            joins = (
                " JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id"
                " JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
            )
            predicates.append("mg.MAG_name = ?")
            query_parameters.append(mag_name)
        if sample_id is not None:
            predicates.append("(f.Sample_id = ? OR f.Sample_id IS NULL)")
            query_parameters.append(sample_id)
        if search_term:
            predicates.append("c.Contig_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        return joins, predicates, query_parameters

    @staticmethod
    def _filtered_sample_scope(parameters, *, contig_id=None, mag_name=None, search_term=""):
        joins = ""
        predicates: list[str] = []
        query_parameters = list(parameters)
        if mag_name is not None:
            joins = (
                " JOIN MAG_contigs_association mca ON mca.Contig_id = f.Contig_id"
                " JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
            )
            predicates.append("mg.MAG_name = ?")
            query_parameters.append(mag_name)
        if contig_id is not None:
            predicates.append("f.Contig_id = ?")
            query_parameters.append(contig_id)
        if search_term:
            predicates.append("s.Sample_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        return joins, predicates, query_parameters

    @staticmethod
    def _filtered_mag_scope(parameters, *, sample_id=None, search_term=""):
        predicates: list[str] = []
        query_parameters = list(parameters)
        if sample_id is not None:
            predicates.append("(f.Sample_id = ? OR f.Sample_id IS NULL)")
            query_parameters.append(sample_id)
        if search_term:
            predicates.append("mg.MAG_name ILIKE '%' || ? || '%'")
            query_parameters.append(search_term)
        return predicates, query_parameters

    def contigs_for_sample(
        self, sample_id: int, search_term: str = "", preserve: str = "", mag_name: str | None = None
    ) -> tuple[str, ...]:
        search, parameters = self._search_clause("c", "Contig_name", search_term)
        mag_joins = ""
        mag_predicate = ""
        if mag_name is not None:
            mag_joins = (
                " JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id"
                " JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
            )
            mag_predicate = " AND mg.MAG_name = ?"
        return self._values(
            "SELECT DISTINCT c.Contig_name FROM Coverage p "
            "JOIN Contig c ON c.Contig_id = p.Contig_id "
            f"{mag_joins} WHERE p.Sample_id = ?{mag_predicate}{search} "
            "ORDER BY (c.Contig_name = ?) DESC, c.Contig_name LIMIT 100",
            [sample_id, *([mag_name] if mag_name is not None else []), *parameters, preserve],
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
        joins, predicates, query_parameters = self._filtered_contig_scope(
            parameters, sample_id=sample_id, mag_name=mag_name, search_term=search_term
        )
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        query_parameters.append(preserve)
        return self._values(
            f"SELECT DISTINCT c.Contig_name FROM ({sql}) f"
            f" JOIN Contig c ON c.Contig_id = f.Contig_id{joins}{where}"
            " ORDER BY (c.Contig_name = ?) DESC, c.Contig_name LIMIT 100",
            query_parameters,
        )

    def count_filtered_contigs(
        self,
        sql: str,
        parameters: list[Any],
        *,
        sample_id: int | None = None,
        mag_name: str | None = None,
    ) -> int:
        joins, predicates, query_parameters = self._filtered_contig_scope(
            parameters, sample_id=sample_id, mag_name=mag_name
        )
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        self.query_count += 1
        return int(
            self.connection.execute(
                f"SELECT COUNT(DISTINCT c.Contig_name) FROM ({sql}) f"
                f" JOIN Contig c ON c.Contig_id = f.Contig_id{joins}{where}",
                query_parameters,
            ).fetchone()[0]
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
        limit: int | None = 100,
    ) -> tuple[str, ...]:
        joins, predicates, query_parameters = self._filtered_sample_scope(
            parameters, contig_id=contig_id, mag_name=mag_name, search_term=search_term
        )
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        query_parameters.extend(preserve)
        limit_clause = " LIMIT 100" if limit is not None else ""
        return self._values(
            f"SELECT DISTINCT s.Sample_name FROM ({sql}) f"
            " JOIN Sample s ON s.Sample_id = f.Sample_id"
            f"{joins}{where}"
            f" ORDER BY (s.Sample_name = ? OR s.Sample_name = ?) DESC, s.Sample_name{limit_clause}",
            query_parameters,
        )

    def count_filtered_samples(
        self,
        sql: str,
        parameters: list[Any],
        *,
        contig_id: int | None = None,
        mag_name: str | None = None,
    ) -> int:
        joins, predicates, query_parameters = self._filtered_sample_scope(
            parameters, contig_id=contig_id, mag_name=mag_name
        )
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        self.query_count += 1
        return int(
            self.connection.execute(
                f"SELECT COUNT(DISTINCT s.Sample_name) FROM ({sql}) f"
                " JOIN Sample s ON s.Sample_id = f.Sample_id"
                f"{joins}{where}",
                query_parameters,
            ).fetchone()[0]
        )

    def filtered_mags(
        self,
        sql: str,
        parameters: list[Any],
        *,
        sample_id: int | None = None,
        search_term: str = "",
    ) -> tuple[str, ...]:
        predicates, query_parameters = self._filtered_mag_scope(
            parameters, sample_id=sample_id, search_term=search_term
        )
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        return self._values(
            f"SELECT DISTINCT mg.MAG_name FROM ({sql}) f"
            " JOIN MAG_contigs_association mca ON mca.Contig_id = f.Contig_id"
            f" JOIN MAG mg ON mg.MAG_id = mca.MAG_id{where} ORDER BY mg.MAG_name",
            query_parameters,
        )

    def count_filtered_mags(
        self, sql: str, parameters: list[Any], *, sample_id: int | None = None
    ) -> int:
        predicates, query_parameters = self._filtered_mag_scope(parameters, sample_id=sample_id)
        where = f" WHERE {' AND '.join(predicates)}" if predicates else ""
        self.query_count += 1
        return int(
            self.connection.execute(
                f"SELECT COUNT(DISTINCT mg.MAG_name) FROM ({sql}) f"
                " JOIN MAG_contigs_association mca ON mca.Contig_id = f.Contig_id"
                f" JOIN MAG mg ON mg.MAG_id = mca.MAG_id{where}",
                query_parameters,
            ).fetchone()[0]
        )

    def mags_for_sample(self, sample_id: int, search_term: str = "") -> tuple[str, ...]:
        search, parameters = self._search_clause("mg", "MAG_name", search_term)
        return self._values(
            "SELECT DISTINCT mg.MAG_name FROM MAG mg "
            "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
            "JOIN Coverage p ON p.Contig_id = mca.Contig_id "
            f"WHERE p.Sample_id = ?{search} ORDER BY mg.MAG_name LIMIT 100",
            [sample_id, *parameters],
        )

    def count_contigs_for_sample(self, sample_id: int, mag_name: str | None = None) -> int:
        self.query_count += 1
        if mag_name is not None:
            return int(
                self.connection.execute(
                    "SELECT COUNT(DISTINCT p.Contig_id) FROM Coverage p"
                    " JOIN MAG_contigs_association mca ON mca.Contig_id = p.Contig_id"
                    " JOIN MAG mg ON mg.MAG_id = mca.MAG_id"
                    " WHERE p.Sample_id = ? AND mg.MAG_name = ?",
                    [sample_id, mag_name],
                ).fetchone()[0]
            )
        return int(
            self.connection
            .execute("SELECT COUNT(DISTINCT Contig_id) FROM Coverage WHERE Sample_id = ?", [sample_id])
            .fetchone()[0]
        )

    def count_samples_for_contig(self, contig_id: int) -> int:
        self.query_count += 1
        return int(
            self.connection
            .execute("SELECT COUNT(DISTINCT Sample_id) FROM Coverage WHERE Contig_id = ?", [contig_id])
            .fetchone()[0]
        )

    def evaluate(self, compiled: CompiledFilter, *, has_mags: bool) -> FilterResult:
        parameters = list(compiled.parameters)
        materialized = compiled
        try:
            self.query_count += 1
            # The same applied relation feeds counts and every availability
            # dropdown. Materialize it once per filter revision instead of
            # re-running a potentially expensive annotation/MAG expression
            # for each consumer.
            self.connection.execute(
                f"CREATE OR REPLACE TEMP TABLE _thebigbam_filter_result AS {compiled.sql}",
                parameters,
            )
            materialized = CompiledFilter("SELECT Contig_id, Sample_id FROM _thebigbam_filter_result", ())
            pair_count, contig_count, sample_count = (
                self.connection
                .execute(
                    "SELECT COUNT(*), COUNT(DISTINCT Contig_id), COUNT(DISTINCT Sample_id) "
                    "FROM _thebigbam_filter_result",
                )
                .fetchone()
            )
            mag_pair_count = None
            if has_mags:
                self.query_count += 1
                mag_pair_count = (
                    self.connection
                    .execute(
                        "WITH _f AS (SELECT * FROM _thebigbam_filter_result) SELECT COUNT(*) FROM ("
                        "SELECT DISTINCT mg.MAG_name, _f.Sample_id FROM _f "
                        "JOIN MAG_contigs_association mca ON mca.Contig_id = _f.Contig_id "
                        "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
                        "WHERE _f.Sample_id IS NOT NULL) counted",
                    )
                    .fetchone()[0]
                )
        except duckdb.Error:
            pair_count = contig_count = sample_count = 0
            mag_pair_count = 0 if has_mags else None
        return FilterResult(
            materialized,
            int(pair_count),
            int(contig_count),
            int(sample_count),
            int(mag_pair_count) if mag_pair_count is not None else None,
        )
