"""Reference-coordinate discrepancy events for sampled read alignments."""

from dataclasses import dataclass

EVENT_VARIABLES = frozenset(("mismatches", "insertions", "deletions", "left_clippings", "right_clippings"))
READ_DISCREPANCY_CHOICES = {
    "Clippings": frozenset(("left_clippings", "right_clippings")),
    "Insertions": frozenset(("insertions",)),
    "Deletions": frozenset(("deletions",)),
    "Mismatches": frozenset(("mismatches",)),
}


@dataclass(frozen=True)
class ReadEvent:
    kind: str
    start: int
    end: int
    length: int
    sequence: str = ""
    reference: str = ""
    clipping: str = ""
