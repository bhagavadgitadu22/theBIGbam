"""BAM-backed fields and helpers for the shared Sample filter category."""

from __future__ import annotations

from typing import Any

from ..models.filters import FilterExpression

BAM_SAMPLE_TABLE = "_thebigbam_bam_samples"
BAM_SAMPLE_COLUMNS = {
    "Sample_name": {"type": "text", "source": "BAMSample", "distinct_values": None},
    "Number_of_reads": {"type": "numeric", "source": "BAMSample", "distinct_values": None},
    "Number_of_mapped_reads": {"type": "numeric", "source": "BAMSample", "distinct_values": None},
    "Circular_mapping": {
        "type": "text", "is_bool": True, "distinct_values": ["yes", "no"], "source": "BAMSample",
    },
}
BAM_COUNT_COLUMNS = frozenset(("Number_of_reads", "Number_of_mapped_reads"))
LEGACY_BAM_COLUMNS = {
    "Sample name": "Sample_name",
    "Number of reads": "Number_of_reads",
    "Number of mapped reads": "Number_of_mapped_reads",
    "Circular mapping": "Circular_mapping",
}


def filtering_metadata_for_mode(metadata: dict[str, Any], *, bam_only: bool) -> dict[str, Any]:
    """Return session-local filtering metadata without changing the database data."""
    result = {category: dict(details) for category, details in metadata.items()}
    if bam_only:
        result["Sample"] = {"source": BAM_SAMPLE_TABLE, "columns": dict(BAM_SAMPLE_COLUMNS)}
    return result


def expression_requires_bam_counts(expression: FilterExpression) -> bool:
    return any(
        predicate.category == "Sample" and predicate.column in BAM_COUNT_COLUMNS
        for section in expression.sections
        for predicate in section.predicates
    )


def migrate_legacy_bam_filtering(settings: dict[str, Any]) -> dict[str, Any]:
    """Fold the former separate BAM filter blocks into the shared expression."""
    legacy = settings.get("bam_sample_filtering") or []
    if not legacy:
        return settings
    migrated = dict(settings)
    filtering = list(settings.get("filtering") or [])
    for index, block in enumerate(legacy):
        rows = []
        for old in block.get("rows", []):
            column = old.get("column")
            normalized_column = LEGACY_BAM_COLUMNS.get(column, column)
            value = old.get("value")
            if normalized_column == "Circular_mapping" and isinstance(value, str):
                value = {"TRUE": "yes", "FALSE": "no"}.get(value.upper(), value)
            rows.append({
                **old,
                "category": "Sample",
                "column": normalized_column,
                "value": value,
            })
        if rows:
            filtering.append({
                "section_and_or": (
                    "AND" if index == 0 else block.get("section_and_or") or "AND"
                ),
                "rows": rows,
            })
    migrated["filtering"] = filtering
    migrated.pop("bam_sample_filtering", None)
    return migrated
