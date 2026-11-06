"""Build display alignment blocks only after candidates have been selected."""

from ..models.read_plot import ReadAlignment
from ..shared.read_cigar import CIGAR
from ..shared.read_identity import read_identity


def summarize_candidate(read, color=None):
    position = read.start
    blocks = []
    for match in CIGAR.finditer(read.cigar):
        length, operation = int(match[1]), match[2]
        if operation in "M=X":
            blocks.append((position, position + length))
            position += length
        elif operation in "DN":
            position += length
    return ReadAlignment(
        read.name, read.start, read.end, tuple(blocks), read.forward, read.alignment_class,
        read.mapq, read.cigar, read.mate_start, read.is_read1, read.is_paired,
        alignment_id=read_identity(read),
        color=color,
    )
