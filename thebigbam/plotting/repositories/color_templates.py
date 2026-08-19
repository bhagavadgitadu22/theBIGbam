"""Read-only access to saved annotation color templates."""

from __future__ import annotations

from typing import Any


class ColorTemplateRepository:
    def __init__(self, connection: Any) -> None:
        self.connection = connection

    def load(self) -> dict[str, list[dict[str, Any]]]:
        """Return rules grouped by template name; absent legacy tables yield no templates."""
        try:
            rows = self.connection.execute(
                "SELECT t.Template_name, r.Qualifier_name, r.Operator, r.Qualifier_value, r.Color "
                "FROM Color_templates t JOIN Color_rules r ON t.Template_id = r.Template_id "
                "ORDER BY t.Template_name, r.Rule_id"
            ).fetchall()
        except Exception:
            return {}
        templates: dict[str, list[dict[str, Any]]] = {}
        for template, qualifier, operator, value, color in rows:
            templates.setdefault(template, []).append(
                {
                    "qualifier_name": qualifier,
                    "operator": operator,
                    "qualifier_value": value,
                    "color": color,
                }
            )
        return templates
