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

    def annotation_dots(self, mag_members, color_rules, max_dots=1000):
        if not color_rules or not mag_members:
            return [], [], 0

        contig_ids = [cid for cid, _length, _offset in mag_members]
        offsets = {cid: offset for cid, _length, offset in mag_members}
        placeholders = ",".join("?" * len(contig_ids))
        rows = self._execute(
            f'SELECT Annotation_id, Contig_id, "Start", "End" FROM Contig_annotation_core '
            f"WHERE Contig_id IN ({placeholders})",
            contig_ids,
        )
        positions = {annotation_id: values for annotation_id, *values in rows}
        colors_by_annotation = {}

        for rule in color_rules:
            key = rule["qualifier_key"]
            mode = rule.get("match_mode", "exact")
            color = rule.get("color", "#ff0000")
            value = rule.get("value")
            if mode in ("has", "has_not"):
                operator = "LIKE" if mode == "has" else "NOT LIKE"
                sql_value = f"%{value}%"
            else:
                operator = _MODE_TO_SQL_OP.get(mode, "=")
                sql_value = value
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
            try:
                matching = self._execute(sql, parameters)
            except Exception as error:
                print(f"[mag_track] Rule query error (key={key}): {error}", flush=True)
                continue
            for (annotation_id,) in matching:
                colors_by_annotation.setdefault(annotation_id, color)

        dots = []
        for annotation_id, color in colors_by_annotation.items():
            if annotation_id in positions:
                contig_id, start, end = positions[annotation_id]
                dots.append((offsets.get(contig_id, 0) + (start + end) / 2, color))
        total = len(dots)
        if max_dots is not None:
            dots = dots[:max_dots]
        return [x for x, _ in dots], [color for _, color in dots], total
