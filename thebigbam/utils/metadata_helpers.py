"""Shared helpers for CSV-to-DuckDB metadata import (sample, contig, MAG)."""


def infer_column_type(values):
    """Infer DuckDB column type from a list of values.

    Tries INTEGER first, then DOUBLE, falls back to TEXT.
    Empty/None values are ignored during inference.

    Args:
        values: List of string values from CSV

    Returns:
        DuckDB type string: 'INTEGER', 'DOUBLE', or 'TEXT'
    """
    # Filter out empty values
    non_empty = [v for v in values if v is not None and v.strip() != '']

    if not non_empty:
        return 'TEXT'  # Default for all-empty columns

    # Try INTEGER
    try:
        for v in non_empty:
            int(v)
        return 'INTEGER'
    except ValueError:
        pass

    # Try DOUBLE
    try:
        for v in non_empty:
            float(v)
        return 'DOUBLE'
    except ValueError:
        pass

    return 'TEXT'


def convert_value(value, col_type):
    """Convert a string value to the appropriate Python type.

    Args:
        value: String value from CSV
        col_type: DuckDB column type ('INTEGER', 'DOUBLE', or 'TEXT')

    Returns:
        Converted value, or None if empty
    """
    if value is None or value.strip() == '':
        return None

    if col_type == 'INTEGER':
        return int(value)
    elif col_type == 'DOUBLE':
        return float(value)
    else:
        return value
