"""Sequence-free alignment candidates retained for repeatable subsampling."""

from dataclasses import dataclass
from functools import cached_property

from ..shared.read_cigar import CIGAR
from ..shared.read_identity import read_score


@dataclass(frozen=True)
class ReadCandidate:
    name: str
    start: int
    end: int
    cigar: str
    alignment_class: str
    flag: int
    mapq: int
    mate_start: int | None
    virtual_offset: int | None = None
    fetch_region: tuple[int, int] | None = None
    template_length: int = 0
    mate_same_reference: bool = True
    matching_bases: int | None = None
    full_read_length: int | None = None

    @property
    def is_read1(self):
        return bool(self.flag & 64)

    @property
    def is_paired(self):
        return bool(self.flag & 1)

    @property
    def forward(self):
        return not self.flag & 16

    @cached_property
    def sampling_score(self):
        return read_score(self)

    @cached_property
    def lengths(self):
        read_length = alignment_length = 0
        for match in CIGAR.finditer(self.cigar):
            length, operation = int(match[1]), match[2]
            if operation in "MISH=X":
                read_length += length
            if operation in "MI=X":
                alignment_length += length
        return read_length, alignment_length

    @property
    def alignment_percentage(self):
        length = self.full_read_length if self.full_read_length is not None else self.lengths[0]
        return 100 * self.matching_bases / length if self.matching_bases is not None and length > 0 else None
