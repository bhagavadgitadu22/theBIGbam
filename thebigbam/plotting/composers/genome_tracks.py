"""Composition facade for database-backed genome tracks."""

from __future__ import annotations

from ..renderers.genome_tracks import GeneMapRenderer, NucleotideRenderer, TranslationRenderer
from ..repositories.genome_tracks import GenomeTrackRepository
from ..services.genome_tracks import GenomeTrackService


class GenomeTrackComposer:
    """Coordinate retrieval, transformation and rendering without mixing their responsibilities."""

    def __init__(self, conn) -> None:
        self.repository = GenomeTrackRepository(conn)
        self.service = GenomeTrackService(self.repository)

    def contig_members(self, contig_name):
        contig_id, length = self.repository.contig_identity(contig_name)
        return ((contig_id, contig_name, length, 0),)

    def mag_members(self, members):
        identities = self.repository.contig_identities([name for name, _length, _offset in members])
        return tuple(
            (identities[name][0], name, int(length), int(offset))
            for name, length, offset in members
            if name in identities
        )

    def gene_map(
        self,
        subject_name,
        subject_length,
        members,
        start,
        end,
        height,
        x_range,
        feature_types=(),
        plot_isoforms=True,
        label_key=None,
        color_rules=(),
        figure_width=30,
    ):
        data = self.service.gene_map(
            subject_name,
            subject_length,
            members,
            start,
            end,
            feature_types or (),
            plot_isoforms,
            label_key,
            color_rules or (),
        )
        return GeneMapRenderer().render(data, height, x_range, figure_width)

    def nucleotides(self, members, start, end, height, x_range):
        return NucleotideRenderer().render(self.service.nucleotides(members, start, end), height, x_range)

    def translation(self, members, start, end, height, x_range):
        data = self.service.translation(members, start, end)
        return TranslationRenderer().render(data, start, end, height, x_range)

    @property
    def diagnostics(self):
        return {"queries": self.repository.query_count, "phases": dict(self.service.phase_seconds)}
