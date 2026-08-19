"""Pure compilation of typed filter expressions into parameterized SQL."""

from __future__ import annotations

from typing import Any, Mapping

from ..models.filters import CompiledFilter, FilterExpression, FilterPredicate, FilterSection

SOURCE_TABLES = {
    "MAG": "MAG",
    "Contig": "Contig",
    "Annotations": "Contig_annotation",
    "Sample": "Sample",
    "MAG coverage": "Explicit_coverage_per_MAG",
    "MAG misassembly": "Explicit_misassembly_per_MAG",
    "MAG microdiversity": "Explicit_microdiversity_per_MAG",
    "Coverage": "Explicit_coverage",
    "Misassembly": "Explicit_misassembly",
    "Microdiversity": "Explicit_microdiversity",
    "Side misassembly": "Explicit_side_misassembly",
    "Topology": "Explicit_topology",
    "Termini": "Explicit_phage_mechanisms",
}
MAG_CATEGORIES = {"MAG", "MAG coverage", "MAG misassembly", "MAG microdiversity"}


class FilterQueryBuilder:
    def __init__(
        self,
        filtering_metadata: Mapping[str, Any],
        filter_encode: Mapping[str, float],
        *,
        has_samples: bool,
    ) -> None:
        self.metadata = filtering_metadata
        self.filter_encode = filter_encode
        self.has_samples = has_samples

    @staticmethod
    def _combine(parts, connectors):
        sql, parameters = parts[0]
        for connector, (right_sql, right_parameters) in zip(connectors, parts[1:]):
            operation = "UNION" if connector == "OR" else "INTERSECT"
            sql = f"({sql}) {operation} ({right_sql})"
            parameters = (*parameters, *right_parameters)
        return sql, parameters

    @staticmethod
    def _qualifier_clause(alias: str, operator: str) -> str:
        if operator == "=":
            return f"('^' || {alias}.\"Value\" || '^') LIKE ('%^' || ? || '^%')"
        if operator == "!=":
            return f"('^' || {alias}.\"Value\" || '^') NOT LIKE ('%^' || ? || '^%')"
        return f'{alias}."Value" {operator} ?'

    def _normalize(self, predicate: FilterPredicate):
        operator, value = predicate.operator, predicate.value
        if operator == "has":
            operator, value = "LIKE", f"%{value}%"
        elif operator == "has not":
            operator, value = "NOT LIKE", f"%{value}%"
        info = self.metadata.get(predicate.category, {}).get("columns", {}).get(predicate.column, {})
        column_type = info.get("type", "numeric")
        is_bool = info.get("is_bool", False)
        if is_bool:
            value = value.lower() == "yes" if isinstance(value, str) else bool(value)
        elif column_type == "text" and not isinstance(value, str):
            value = str(value)
        elif column_type == "numeric" and isinstance(value, str):
            try:
                value = float(value)
            except (TypeError, ValueError):
                return None
        encoding = self.filter_encode.get(predicate.column)
        if encoding and not is_bool and column_type == "numeric" and predicate.category in ("Contig", "Sample", "MAG"):
            value = round(value * encoding)
        return operator, value, info

    def predicate(self, predicate: FilterPredicate):
        source = SOURCE_TABLES.get(predicate.category)
        normalized = self._normalize(predicate)
        if source is None or normalized is None:
            return None
        operator, value, info = normalized
        column = predicate.column
        column_source = info.get("source")
        qualifier_key = info.get("qualifier_key")
        parameters = (value,)
        sample_joins = (
            "LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id LEFT JOIN Sample s ON p.Sample_id = s.Sample_id "
            if self.has_samples
            else ""
        )
        sample_id = "s.Sample_id" if self.has_samples else "NULL"
        if column_source == "Contig_annotation":
            sql = (
                f"SELECT DISTINCT c.Contig_id, {sample_id} FROM Contig_annotation ca "
                f"JOIN Contig c ON ca.Contig_id = c.Contig_id {sample_joins}"
                f'WHERE ca."{column}" {operator} ?'
            )
        elif column_source in {"Annotation_qualifier", "Contig_qualifier"}:
            alias = "aq" if column_source == "Annotation_qualifier" else "cq"
            base = (
                "FROM Annotation_qualifier aq "
                "JOIN Contig_annotation_core cac ON aq.Annotation_id = cac.Annotation_id "
                "JOIN Contig c ON cac.Contig_id = c.Contig_id "
                if alias == "aq"
                else "FROM Contig_qualifier cq JOIN Contig c ON cq.Contig_id = c.Contig_id "
            )
            sql = (
                f"SELECT DISTINCT c.Contig_id, {sample_id} {base}{sample_joins}"
                f'WHERE {alias}."Key" = ? AND {self._qualifier_clause(alias, operator)}'
            )
            parameters = (qualifier_key, value)
        elif predicate.category == "Contig":
            sql = (
                f'SELECT DISTINCT c.Contig_id, {sample_id} FROM Contig c {sample_joins}WHERE c."{column}" {operator} ?'
            )
        elif predicate.category == "Sample":
            sql = (
                "SELECT DISTINCT c.Contig_id, s.Sample_id FROM Sample s "
                "LEFT JOIN Coverage p ON s.Sample_id = p.Sample_id "
                "LEFT JOIN Contig c ON p.Contig_id = c.Contig_id "
                f'WHERE s."{column}" {operator} ?'
            )
        elif predicate.category == "MAG":
            joins = (
                "LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id LEFT JOIN Sample s ON p.Sample_id = s.Sample_id "
                if self.has_samples
                else ""
            )
            sql = (
                f"SELECT DISTINCT c.Contig_id, {sample_id} FROM MAG mg "
                "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
                f"JOIN Contig c ON c.Contig_id = mca.Contig_id {joins}"
                f'WHERE mg."{column}" {operator} ?'
            )
        elif predicate.category in MAG_CATEGORIES:
            sql = (
                "SELECT DISTINCT c.Contig_id, s.Sample_id "
                f"FROM {source} v JOIN Sample s ON s.Sample_name = v.Sample_name "
                "JOIN MAG mg ON mg.MAG_name = v.MAG_name "
                "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
                "JOIN Contig c ON c.Contig_id = mca.Contig_id "
                f'WHERE v."{column}" {operator} ?'
            )
        else:
            sql = (
                "SELECT DISTINCT c.Contig_id, s.Sample_id "
                f"FROM {source} v JOIN Contig c ON c.Contig_name = v.Contig_name "
                "JOIN Sample s ON s.Sample_name = v.Sample_name "
                f'WHERE v."{column}" {operator} ?'
            )
        return sql, parameters

    def section(self, section: FilterSection):
        parts = []
        connectors = []
        for index, predicate in enumerate(section.predicates):
            compiled = self.predicate(predicate)
            if compiled is None:
                continue
            if parts:
                connectors.append(section.connectors[index - 1])
            parts.append(compiled)
        if not parts:
            return None
        return self._combine(parts, connectors)

    def compile(self, expression: FilterExpression) -> CompiledFilter | None:
        parts = []
        connectors = []
        for index, section in enumerate(expression.sections):
            compiled = self.section(section)
            if compiled is None:
                continue
            if parts:
                connectors.append(expression.connectors[index - 1])
            parts.append(compiled)
        if not parts:
            return None
        sql, parameters = self._combine(parts, connectors)
        return CompiledFilter(sql, tuple(parameters))
