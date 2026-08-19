"""DuckDB retrieval for the MAG overview track."""

from __future__ import annotations

_CORE_ANNOTATION_COLUMNS = {"Start", "End", "Strand", "Type", "Main_isoform"}
_MODE_TO_SQL_OP = {"exact": "=", "not_equal": "!=", "lt": "<", "gt": ">"}


class MagOverviewRepository:
    """Retrieve annotation dots without importing Bokeh or Panel."""

    def __init__(self, connection) -> None:
        self.connection = connection
        self.query_count = 0

    def _execute(self, sql, parameters):
        self.query_count += 1
        return self.connection.cursor().execute(sql, parameters).fetchall()

    def annotation_positions(self, contig_ids):
        if not contig_ids:
            return ()
        placeholders = ",".join("?" * len(contig_ids))
        return tuple(
            self._execute(
                f'SELECT Annotation_id, Contig_id, "Start", "End" FROM Contig_annotation_core '
                f"WHERE Contig_id IN ({placeholders})",
                contig_ids,
            )
        )

    def matching_annotation_ids(self, rule, contig_ids):
        key = rule["qualifier_key"]
        mode = rule.get("match_mode", "exact")
        value = rule.get("value")
        operator = (
            ("LIKE" if mode == "has" else "NOT LIKE") if mode in ("has", "has_not") else _MODE_TO_SQL_OP.get(mode, "=")
        )
        sql_value = f"%{value}%" if mode in ("has", "has_not") else value
        placeholders = ",".join("?" * len(contig_ids))
        token_clause = {
            "=": "('^' || aq.Value || '^') LIKE ('%^' || ? || '^%')",
            "!=": "('^' || aq.Value || '^') NOT LIKE ('%^' || ? || '^%')",
        }.get(operator, f"aq.Value {operator} ?")
        if key in _CORE_ANNOTATION_COLUMNS:
            sql = (
                f'SELECT Annotation_id FROM Contig_annotation_core WHERE "{key}" {operator} ? '
                f"AND Contig_id IN ({placeholders})"
            )
            parameters = [sql_value, *contig_ids]
        else:
            sql = (
                "SELECT a.Annotation_id FROM Contig_annotation_core a "
                "JOIN Annotation_qualifier aq ON aq.Annotation_id = a.Annotation_id "
                f"WHERE aq.Key = ? AND {token_clause} AND a.Contig_id IN ({placeholders})"
            )
            parameters = [key, sql_value, *contig_ids]
        return tuple(annotation_id for (annotation_id,) in self._execute(sql, parameters))
