"""Read-filter fields and operator contracts shared by controls and evaluation."""

READ_FILTER_FIELDS = (
    "Read length", "Alignment length", "Alignment percentage", "Alignment start", "Alignment end", "Alignment includes",
    "MAPQ", "Insert sizes", "Non-inward pair", "Has mate mapped", "Circular",
)
BOOLEAN_READ_FILTERS = frozenset(("Non-inward pair", "Has mate mapped", "Circular"))
READ_FILTER_LABELS = {"Alignment percentage": "Alignment (%)"}


def read_filter_operators(field):
    if field in BOOLEAN_READ_FILTERS:
        return ["="]
    return ["=", "!="] if field == "Alignment includes" else [">", "<", "=", "!="]
