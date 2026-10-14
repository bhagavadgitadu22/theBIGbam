"""Read-only access to saved annotation color templates."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


EMPATHI_COLORS = {
    "adsorption-related": "#feb2ff",
    "cell_wall_depolymerase": "#f37c76",
    "super_infection": "#dff1e2",
    "crispr": "#dff1e2",
    "anti-restriction": "#dff1e2",
    "toxin": "#dff1e2",
    "sir2": "#dff1e2",
    "DNA-associated": "#f9ff96",
    "ejection": "#c4ffff",
    "lysis": "#fad1dc",
    "nucleotide_metabolism": "#e5c492",
    "packaging_assembly": "#abaafe",
    "phosphorylation": "#e384ff",
    "pvp": "#8ac7ee",
    "reductase": "#afafaf",
    "RNA-associated": "#ddff69",
    "transcriptional_regulator": "#fad082",
    "transferase": "#c59973",
    "unknown": "#d8dfdf",
}


class ColorTemplateRepository:
    def __init__(self, connection: Any) -> None:
        self.connection = connection

    def load(self, annotation_columns: Iterable[str] = ()) -> dict[str, list[dict[str, Any]]]:
        """Load saved rules and offer Empathi for databases with its annotation column."""
        try:
            rows = self.connection.execute(
                "SELECT t.Template_name, r.Qualifier_name, r.Operator, r.Qualifier_value, r.Color "
                "FROM Color_templates t JOIN Color_rules r ON t.Template_id = r.Template_id "
                "ORDER BY t.Template_name, r.Rule_id"
            ).fetchall()
        except Exception:
            rows = []
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
        if "empathi_Annotation" in annotation_columns:
            templates.setdefault("empathi", [
                {
                    "qualifier_name": "empathi_Annotation",
                    "operator": "has",
                    "qualifier_value": value,
                    "color": color,
                }
                for value, color in EMPATHI_COLORS.items()
            ])
        return templates
