"""Shared CIGAR tokenizer for read metadata, summaries and discrepancies."""

import re

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
