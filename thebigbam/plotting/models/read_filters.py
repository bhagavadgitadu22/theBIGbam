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


def read_color_rule_operators(field):
    """Operators available to color rules, independently of read filters."""
    if field == "Alignment includes":
        return ["=", "!="]
    if field in READ_FILTER_CHOICES:
        return ["=", "!=", "has", "has not", "Use random colors"]
    if field in BOOLEAN_READ_FILTERS:
        return ["=", "!=", "Use random colors"]
    return [">", "<", "=", "!=", "Use random colors", "Use color gradient"]
