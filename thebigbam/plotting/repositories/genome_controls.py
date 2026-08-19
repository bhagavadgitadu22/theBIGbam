"""Read-only capability checks needed while constructing genome controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class GenomeControlCapabilities:
    has_isoforms: bool
    has_sequence: bool
    has_translation: bool


class GenomeControlRepository:
    def __init__(self, connection: Any) -> None:
        self.connection = connection

    def capabilities(self) -> GenomeControlCapabilities:
        isoform = self.connection.execute(
            "SELECT Status FROM Constants_for_plotting WHERE Constant = 'isoforms'"
        ).fetchone()
        sequence = self.connection.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = 'Contig_sequence'"
        ).fetchone()
        translation = self.connection.execute(
            "SELECT 1 FROM information_schema.columns "
            "WHERE table_name = 'Contig_annotation' AND column_name = 'Protein_sequence'"
        ).fetchone()
        return GenomeControlCapabilities(
            has_isoforms=bool(isoform and isoform[0]),
            has_sequence=sequence is not None,
            has_translation=translation is not None,
        )
