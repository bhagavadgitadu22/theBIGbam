"""Read-filter fields and operator contracts shared by controls and evaluation."""

READ_FILTER_FIELDS = (
    "Alignment class", "Alignment end", "Alignment includes", "Alignment length", "Alignment start",
    "MAPQ", "Read % aligned", "Read length", "Strand",
    "Has mate mapped", "Insert sizes", "Non-inward pair", "Circular",
)
BOOLEAN_READ_FILTERS = frozenset(("Non-inward pair", "Has mate mapped", "Circular"))
READ_FILTER_CHOICES = {
    "Strand": ("Forward", "Reverse"),
    "Alignment class": ("Primary", "Secondary", "Supplementary"),
}
SINGLE_VALUE_READ_FILTERS = BOOLEAN_READ_FILTERS | frozenset(READ_FILTER_CHOICES)
READ_FILTER_LABELS = {"Alignment percentage": "Read % aligned"}


def read_filter_operators(field):
    if field in SINGLE_VALUE_READ_FILTERS:
        return ["="]
    return ["=", "!="] if field == "Alignment includes" else [">", "<", "=", "!="]
