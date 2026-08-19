"""Shared feature selection and loading used by plot composers."""

from __future__ import annotations

from ..models.blobs import FeatureLoadRequest, FeatureRegion, MagFeatureLoadRequest, MagMember
from ..repositories.features import split_features_by_storage
from ..shared.defaults import DEFAULT_MAX_BASE_RESOLUTION


def get_variable_metadata_batch(repository, subplot_list):
    return repository.metadata_batch(subplot_list)


def split_contig_vs_sample_features(metadata_cache, requested_features):
    return split_features_by_storage(metadata_cache, requested_features)


def get_feature_data(
    repository,
    feature,
    contig_id,
    sample_id,
    xstart=None,
    xend=None,
    variable_metadata=None,
    max_base_resolution=None,
    min_relative_value=0.0,
    window_for_zoom=None,
    enable_timing=False,
    encoding_by_feature=None,
    profiler=None,
):
    """Load one contig feature through the repository/service pipeline."""
    from .blob_features import FeatureDataService

    end = int(xend if xend is not None else repository.contig_length(contig_id))
    request = FeatureLoadRequest(
        feature=feature,
        contig_id=int(contig_id),
        sample_id=int(sample_id) if sample_id is not None else None,
        region=FeatureRegion(int(xstart if xstart is not None else 1), end),
        max_base_resolution=int(max_base_resolution or DEFAULT_MAX_BASE_RESOLUTION),
        minimum_relative_value=float(min_relative_value),
        encoding_by_feature=encoding_by_feature,
    )
    service = FeatureDataService(repository)
    result = service.load_contig(request, variable_metadata)
    if profiler is not None:
        for name, seconds in service.phase_seconds.items():
            profiler.seconds[name] = profiler.seconds.get(name, 0.0) + seconds
    if enable_timing:
        print(
            f"[timing] feature pipeline '{feature}': queries={repository.query_count} "
            f"points={service.decoded_points} phases={service.phase_seconds}",
            flush=True,
        )
    return result


def get_mag_feature_data(
    repository,
    feature,
    mag_id,
    sample_id,
    mag_length,
    xstart=None,
    xend=None,
    variable_metadata=None,
    max_base_resolution=None,
    min_relative_value=0.0,
    members=None,
    encoding_by_feature=None,
    profiler=None,
):
    """Load one MAG feature through the repository/service pipeline."""
    from .blob_features import FeatureDataService

    if members is None:
        raise ValueError("MAG members must be loaded by a repository before service transformation")
    request = MagFeatureLoadRequest(
        feature=feature,
        sample_id=int(sample_id) if sample_id is not None else None,
        region=FeatureRegion(int(xstart if xstart is not None else 1), int(xend if xend is not None else mag_length)),
        members=tuple(MagMember(int(cid), int(length), int(offset)) for cid, length, offset in members),
        max_base_resolution=int(max_base_resolution or DEFAULT_MAX_BASE_RESOLUTION),
        minimum_relative_value=float(min_relative_value),
        encoding_by_feature=encoding_by_feature,
    )
    service = FeatureDataService(repository)
    result = service.load_mag(request, variable_metadata)
    if profiler is not None:
        for name, seconds in service.phase_seconds.items():
            profiler.seconds[name] = profiler.seconds.get(name, 0.0) + seconds
        profiler.metric("mag_members_total", service.mag_members_total)
        profiler.metric("mag_members_loaded", service.mag_members_loaded)
    return result


def parse_requested_features(features):
    """Expand UI module names into deduplicated subplot names."""
    module_features = {
        "genome": ("Repeat count", "Max repeat identity", "GC content", "GC skew"),
        "coverage": ("Primary alignments", "Other alignments"),
        "phage termini": (
            "Coverage reduced",
            "Reads termini",
            "Read termini transformation",
            "Repeat count",
            "Max repeat identity",
        ),
        "assembly check": ("Clippings", "Indels", "Mismatches", "Read lengths", "Insert sizes", "Bad orientations"),
    }
    aliases = {
        "phagetermini": "phage termini",
        "phage_termini": "phage termini",
        "assemblycheck": "assembly check",
        "assembly_check": "assembly check",
        "repeats": "repeats",
        "repeat": "repeats",
        "duplications": "repeats",
        "duplication": "repeats",
        "gc_content": "gc content",
        "gccontent": "gc content",
        "gc": "gc content",
        "gc_skew": "gc skew",
        "gcskew": "gc skew",
        "skew": "gc skew",
    }
    expanded = []
    for feature in features:
        normalized = aliases.get(feature.lower().strip(), feature.lower().strip())
        if normalized == "repeats":
            expanded.extend(("Repeat count", "Max repeat identity"))
        elif normalized == "gc content":
            expanded.append("GC content")
        elif normalized == "gc skew":
            expanded.append("GC skew")
        elif normalized == "repeat count":
            expanded.append("Repeat count")
        elif normalized == "max repeat identity":
            expanded.append("Max repeat identity")
        elif normalized in module_features:
            expanded.extend(module_features[normalized])
        else:
            expanded.append(feature)
    return list(dict.fromkeys(expanded))
