from thebigbam.plotting.services.composition_data import parse_requested_features


def test_requested_feature_modules_and_legacy_aliases_expand_once():
    assert parse_requested_features(
        ["Coverage", "phage_termini", "duplications", "gc", "GC skew", "Coverage"]
    ) == [
        "Primary alignments",
        "Other alignments",
        "Coverage reduced",
        "Reads termini",
        "Read termini transformation",
        "Repeat count",
        "Max repeat identity",
        "GC content",
        "GC skew",
    ]


def test_unknown_feature_names_preserve_original_spelling():
    assert parse_requested_features(["My custom subplot"]) == ["My custom subplot"]
