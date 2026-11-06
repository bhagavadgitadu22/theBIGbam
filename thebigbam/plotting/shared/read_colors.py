"""Stable colors and built-in coloring templates for BAM alignments."""

READ_CLASS_COLORS = {
    "primary_reads": "#999999",
    "secondary_reads": "#1f77b4",
    "supplementary_reads": "#B45C1F",
    "primary_reads_plus_only": "#a1665e",
    "primary_reads_minus_only": "#5E99A1",
}

READ_CLASS_TEMPLATE_COLORS = {
    "Primary": READ_CLASS_COLORS["primary_reads"],
    "Secondary": READ_CLASS_COLORS["secondary_reads"],
    "Supplementary": READ_CLASS_COLORS["supplementary_reads"],
}

READ_STRAND_VARIABLES = {
    "Forward": "primary_reads_plus_only",
    "Reverse": "primary_reads_minus_only",
}


def read_color_templates(strand_colors=None):
    """Return independent built-in rule sets using the database strand palette."""
    strand_colors = {**READ_CLASS_COLORS, **(strand_colors or {})}
    return {
        "Alignment class": [
            {"field": "Alignment class", "operator": "=", "value": alignment_class,
             "color": color}
            for alignment_class, color in READ_CLASS_TEMPLATE_COLORS.items()
        ],
        "Strand": [
            {"field": "Strand", "operator": "=", "value": strand,
             "color": strand_colors[variable]}
            for strand, variable in READ_STRAND_VARIABLES.items()
        ],
    }
