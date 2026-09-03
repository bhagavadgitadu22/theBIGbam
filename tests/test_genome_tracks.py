from bokeh.models import Range1d

from thebigbam.plotting.renderers.genome_tracks import TranslationRenderer
from thebigbam.plotting.services.genome_tracks import (
    CdsData,
    GenomeTrackService,
    TranslationData,
)


class FakeGenomeRepository:
    query_count = 2

    def annotations(self, contig_ids, offsets, start, end, feature_types, plot_isoforms):
        assert contig_ids == [10, 20]
        assert offsets == {10: 0, 20: 100}
        return [(7, 105, 113, 1, "CDS", True, "enzyme", None, "gene_7")]

    def qualifiers(self, annotation_ids, keys):
        assert annotation_ids == [7]
        assert keys == ("score",)
        return {7: {"score": 5}}

    def sequence(self, name, local_start, length):
        return {"first": "AC", "second": "GT"}[name]

    def cds(self, name, start, end):
        if name != "second":
            return []
        return [(5, 13, 1, "ATGAAATTT", "MKF", [{"start": 5, "end": 13}])]

    def codons(self):
        return {"ATG": ("M", "Methionine", "Met", "#111111")}


def test_service_uses_canonical_mag_offsets_and_batches_qualifiers():
    service = GenomeTrackService(FakeGenomeRepository())
    members = ((10, "first", 100, 0), (20, "second", 100, 100))
    data = service.gene_map(
        "mag",
        200,
        members,
        101,
        120,
        label_key="score",
        color_rules=({"qualifier_key": "score", "match_mode": "gt", "value": 4, "color": "#abcdef"},),
    )
    assert data.features[0].start == 105
    assert data.features[0].qualifiers["score"] == 5
    assert data.features[0].color == "#abcdef"


def test_nucleotide_and_translation_coordinates_span_mag_members():
    service = GenomeTrackService(FakeGenomeRepository())
    members = ((10, "first", 2, 0), (20, "second", 2, 2))
    sequence = service.nucleotides(members, 1, 4)
    assert sequence.positions == (1, 2, 3, 4)
    assert sequence.nucleotides == ("A", "C", "G", "T")

    translation = service.translation(((20, "second", 100, 100),), 101, 120)
    assert translation.entries[0].segments == ((105, 113),)


def test_translation_renderer_renders_prepared_data():
    data = TranslationData(
        (CdsData(1, 3, 1, "ATG", "M", ((1, 3),)),),
        {"ATG": ("M", "Methionine", "Met", "#111111")},
    )
    plot = TranslationRenderer().render(data, 1, 3, 80, Range1d(1, 3))
    assert plot is not None
