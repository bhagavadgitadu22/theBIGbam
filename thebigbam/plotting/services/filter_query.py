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
RAW_METRIC_TABLES = {
    "Coverage": "Coverage",
    "Misassembly": "Misassembly",
    "Microdiversity": "Microdiversity",
    "Side misassembly": "Side_misassembly",
}


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
        sql, parameters, *initial_scope = parts[0]
        scope = initial_scope[0] if initial_scope else ("wildcard" if "NULL AS Sample_id" in sql else "sample")
        for connector, right_part in zip(connectors, parts[1:]):
            right_sql, right_parameters, *declared_scope = right_part
            right_scope = declared_scope[0] if declared_scope else (
                "wildcard" if "NULL AS Sample_id" in right_sql else "sample"
            )
            if connector == "OR":
                sql = f"({sql}) UNION ({right_sql})"
                scope = scope if scope == right_scope else "mixed"
            else:
                # A null-safe pair join has the same set semantics as
                # INTERSECT, while allowing DuckDB to push the selective side
                # of mixed annotation/MAG filters into the other relation.
                if scope == "wildcard" and right_scope != "mixed":
                    sample_projection = "_right.Sample_id"
                    sample_match = ""
                    scope = right_scope
                elif right_scope == "wildcard" and scope != "mixed":
                    sample_projection = "_left.Sample_id"
                    sample_match = ""
                else:
                    sample_projection = "COALESCE(_left.Sample_id, _right.Sample_id)"
                    sample_match = (
                        " AND (_left.Sample_id IS NULL OR _right.Sample_id IS NULL OR "
                        "_left.Sample_id = _right.Sample_id)"
                    )
                    scope = "mixed" if "mixed" in {scope, right_scope} else "sample"
                sql = (
                    f"SELECT DISTINCT _left.Contig_id, {sample_projection} AS Sample_id "
                    f"FROM ({sql}) _left JOIN ({right_sql}) _right ON "
                    "_left.Contig_id IS NOT DISTINCT FROM _right.Contig_id"
                    f"{sample_match}"
                )
            parameters = (*parameters, *right_parameters)
        return sql, parameters, scope

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
        if encoding and not is_bool and column_type == "numeric" and (
            predicate.category in ("Contig", "Sample", "MAG")
            or predicate.category in RAW_METRIC_TABLES
        ):
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
        # Contig-scoped predicates use NULL as an internal wildcard sample.
        # AND composition can then bind it directly to a sample-scoped metric
        # without first multiplying the contig across all Coverage rows.
        sample_id = "NULL"
        if column_source == "Contig_annotation":
            sql = (
                f"SELECT DISTINCT ca.Contig_id, {sample_id} AS Sample_id FROM Contig_annotation ca "
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
                f"SELECT DISTINCT c.Contig_id, {sample_id} AS Sample_id {base}"
                f'WHERE {alias}."Key" = ? AND {self._qualifier_clause(alias, operator)}'
            )
            parameters = (qualifier_key, value)
        elif predicate.category == "Contig":
            sql = (
                f'SELECT DISTINCT c.Contig_id, {sample_id} AS Sample_id FROM Contig c '
                f'WHERE c."{column}" {operator} ?'
            )
        elif predicate.category == "Sample":
            sql = (
                "SELECT DISTINCT c.Contig_id, s.Sample_id FROM Sample s "
                "LEFT JOIN Coverage p ON s.Sample_id = p.Sample_id "
                "LEFT JOIN Contig c ON p.Contig_id = c.Contig_id "
                f'WHERE s."{column}" {operator} ?'
            )
        elif predicate.category == "MAG":
            sql = (
                f"SELECT DISTINCT c.Contig_id, {sample_id} AS Sample_id FROM MAG mg "
                "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
                "JOIN Contig c ON c.Contig_id = mca.Contig_id "
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
        elif predicate.category in RAW_METRIC_TABLES:
            raw_source = RAW_METRIC_TABLES[predicate.category]
            sql = (
                "SELECT DISTINCT v.Contig_id, v.Sample_id "
                f"FROM {raw_source} v "
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
        sql, parameters, scope = self._combine(parts, connectors)
        if self.has_samples and scope != "sample":
            # Resolve any contig-scoped wildcard only after all AND/OR
            # composition. Sample-specific predicates have already replaced
            # it through COALESCE, so mixed annotation/MAG filters avoid the
            # all-samples intermediate expansion.
            sql = (
                "SELECT DISTINCT _result.Contig_id, "
                "COALESCE(_result.Sample_id, _coverage.Sample_id) AS Sample_id "
                f"FROM ({sql}) _result LEFT JOIN Coverage _coverage ON "
                "_result.Sample_id IS NULL AND _coverage.Contig_id = _result.Contig_id"
            )
        return CompiledFilter(sql, tuple(parameters))
