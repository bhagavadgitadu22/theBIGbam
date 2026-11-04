"""Discover indexed BAMs and match them to plotting samples."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import pysam


@dataclass(frozen=True)
class BamSample:
    name: str
    path: Path
    index_path: Path | None
    reference_lengths: Mapping[str, int]
    warnings: tuple[str, ...] = ()

    @property
    def available(self) -> bool:
        return self.index_path is not None and bool(self.reference_lengths)

    def warning_for_contig(self, contig: str, expected_length: int) -> str | None:
        actual_length = self.reference_lengths.get(contig)
        if actual_length is None:
            return f"BAM '{self.path.name}' has no alignment reference named '{contig}'."
        if actual_length != expected_length:
            return (
                f"BAM '{self.path.name}' has reference '{contig}' of length {actual_length:,} bp; "
                f"the database expects {expected_length:,} bp."
            )
        return None


@dataclass(frozen=True)
class BamCatalog:
    by_sample: Mapping[str, BamSample]
    warnings: tuple[str, ...]

    def warning_for_sample(self, sample_name: str, contig: str | None = None, expected_length: int = 0) -> str | None:
        sample = self.by_sample.get(sample_name)
        if sample is None:
            return f"No BAM file matches database sample '{sample_name}'."
        if sample.warnings:
            return " ".join(sample.warnings)
        if contig is not None:
            return sample.warning_for_contig(contig, expected_length)
        return None


def sample_name_from_bam(path: Path) -> str:
    """Use the same BAM-to-sample name convention as calculate."""
    return path.stem.replace("_with_MD", "")


def _bam_paths(source: Path) -> tuple[Path, ...]:
    if source.is_dir():
        return tuple(sorted(path for path in source.iterdir() if path.is_file() and path.suffix == ".bam"))
    if not source.is_file():
        raise ValueError(f"BAM path not found: {source}")
    if source.suffix != ".bam":
        raise ValueError(f"Expected a .bam file: {source}")
    return (source,)


def _index_path(path: Path) -> Path | None:
    candidates = (Path(f"{path}.bai"), path.with_suffix(".bai"), Path(f"{path}.csi"), path.with_suffix(".csi"))
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def build_bam_catalog(
    source: str | Path,
    db_samples: tuple[str, ...],
    contig_lengths: Mapping[str, int],
) -> BamCatalog:
    """Validate supplied BAMs without retaining file handles or reading alignments."""
    paths = _bam_paths(Path(source))
    warnings: list[str] = []
    if not paths:
        warnings.append(f"No .bam files found in directory '{source}'.")
    by_sample: dict[str, BamSample] = {}
    db_names = set(db_samples)
    for path in paths:
        name = sample_name_from_bam(path)
        if name in by_sample:
            raise ValueError(f"Multiple BAM files map to sample '{name}': '{by_sample[name].path}' and '{path}'.")
        issues: list[str] = []
        index_path = _index_path(path)
        if index_path is None:
            issues.append(f"BAM '{path.name}' has no .bai or .csi index.")
        reference_lengths: dict[str, int] = {}
        try:
            with pysam.AlignmentFile(str(path), "rb") as bam:
                reference_lengths = dict(zip(bam.references, bam.lengths))
                if index_path is not None and not bam.has_index():
                    issues.append(f"BAM '{path.name}' has an index file that cannot be read.")
                    index_path = None
        except (OSError, ValueError) as error:
            issues.append(f"BAM '{path.name}' cannot be opened: {error}")
            index_path = None
        if not reference_lengths:
            issues.append(f"BAM '{path.name}' has no reference sequences.")
        if name in db_names and reference_lengths:
            shared = reference_lengths.keys() & contig_lengths.keys()
            if not shared:
                issues.append(f"BAM '{path.name}' has no reference contig matching the database.")
            mismatched = sorted(
                contig for contig in shared if reference_lengths[contig] != contig_lengths[contig]
            )
            if mismatched:
                preview = ", ".join(mismatched[:3])
                suffix = "..." if len(mismatched) > 3 else ""
                issues.append(
                    f"BAM '{path.name}' has {len(mismatched)} reference length mismatch(es): {preview}{suffix}."
                )
        record = BamSample(name, path, index_path, reference_lengths, tuple(issues))
        by_sample[name] = record
        warnings.extend(issues)
        if name not in db_names:
            warnings.append(f"BAM '{path.name}' has no matching database sample '{name}'.")
    warnings.extend(
        f"No BAM file matches database sample '{name}'." for name in sorted(db_names - by_sample.keys())
    )
    return BamCatalog(by_sample, tuple(warnings))
