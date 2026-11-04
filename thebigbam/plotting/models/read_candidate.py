"""Sequence-free alignment candidates retained for repeatable subsampling."""

from dataclasses import dataclass
from functools import cached_property

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
