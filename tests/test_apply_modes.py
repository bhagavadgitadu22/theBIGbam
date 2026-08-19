from types import SimpleNamespace

from thebigbam.plotting.application.apply_modes import (
    all_sample_features,
    one_sample_features,
    track_visibility,
    without_gene_map,
)


def _group(labels, active):
    return SimpleNamespace(labels=labels, active=active)


def test_feature_projection_distinguishes_one_and_all_sample_modes():
    genome = _group(["Gene map", "Repeats"], [0, 1])
    variables = [_group(["Depth", "Breadth"], [0, 1]), _group(["Identity"], [0])]

    assert one_sample_features(genome, variables) == ["Gene map", "Repeats", "Depth", "Breadth", "Identity"]
    assert all_sample_features(genome, variables) == (["Gene map", "Repeats"], "Breadth")


def test_track_visibility_and_gene_map_filtering_use_window_thresholds():
    visibility = track_visibility(
        500,
        max_genemap_window=400,
        max_sequence_window=600,
        sequence_requested=True,
        translation_requested=False,
    )

    assert not visibility.gene_map
    assert visibility.sequence
    assert not visibility.translation
    assert without_gene_map(["Gene map", "Repeats"], visibility.gene_map) == ["Repeats"]
