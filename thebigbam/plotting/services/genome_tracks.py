"""Database and biological data layer for annotation and sequence tracks.

No Bokeh or Panel types are used here.
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
from typing import Any, Mapping

from ..repositories.genome_tracks import GenomeTrackRepository


@dataclass(frozen=True)
class AnnotationFeatureData:
    annotation_id: int
    start: int
    end: int
    strand: int
    feature_type: str
    qualifiers: Mapping[str, Any] = field(default_factory=dict)
    color: str | None = None
    tooltip_key: str | None = None


@dataclass(frozen=True)
class GeneMapData:
    subject_name: str
    subject_length: int
    features: tuple[AnnotationFeatureData, ...]


@dataclass(frozen=True)
class NucleotideSequenceData:
    positions: tuple[int, ...]
    nucleotides: tuple[str, ...]


@dataclass(frozen=True)
class CdsData:
    start: int
    end: int
    strand: int
    nucleotide_sequence: str
    protein_sequence: str
    segments: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class TranslationData:
    entries: tuple[CdsData, ...]
    codons: Mapping[str, tuple[str, str, str, str]]


def _normalize_segments(segments, start, end, offset=0):
    if not segments:
        return ((int(start) + offset, int(end) + offset),)
    normalized = sorted((int(item["start"]) + offset, int(item["end"]) + offset) for item in segments)
    return tuple(normalized)


def _matches(value, rule):
    mode = rule.get("match_mode", "exact")
    expected = rule.get("value")
    if mode == "random":
        return value is not None
    if mode == "has":
        return expected is not None and str(expected).lower() in str(value or "").lower()
    if mode == "has_not":
        return expected is not None and str(expected).lower() not in str(value or "").lower()
    if mode == "not_equal":
        return value != expected
    if mode in {"lt", "gt"}:
        try:
            return float(value) < float(expected) if mode == "lt" else float(value) > float(expected)
        except (TypeError, ValueError):
            return False
    return value == expected or str(value) == str(expected)


def _rule_color(value, rule):
    if rule.get("match_mode") == "random":
        digest = hashlib.sha256(str(value).encode()).hexdigest()
        return f"#{digest[:6]}"
    return rule.get("color")


class GenomeTrackService:
    def __init__(self, repository: GenomeTrackRepository) -> None:
        self.repository = repository
        self.phase_seconds: dict[str, float] = {}

    def _record(self, name, started):
        self.phase_seconds[name] = self.phase_seconds.get(name, 0.0) + time.perf_counter() - started

    def gene_map(
        self,
        subject_name,
        subject_length,
        members,
        start,
        end,
        feature_types=(),
        plot_isoforms=True,
        label_key=None,
        color_rules=(),
    ):
        started = time.perf_counter()
        contig_ids = [int(member[0]) for member in members]
        offsets = {int(member[0]): int(member[3]) for member in members}
        rows = self.repository.annotations(contig_ids, offsets, start, end, feature_types, plot_isoforms)
        annotation_ids = [int(row[0]) for row in rows]
        rule_keys = [rule.get("qualifier_key") for rule in color_rules if rule.get("qualifier_key")]
        wanted_keys = tuple(dict.fromkeys([key for key in [label_key, *rule_keys] if key]))
        extra = self.repository.qualifiers(annotation_ids, wanted_keys)
        features = []
        for aid, feature_start, feature_end, strand, feature_type, main_isoform, product, function, locus in rows:
            qualifiers = {
                "product": product,
                "function": function,
                "locus_tag": locus,
                "Main_isoform": main_isoform,
                "Type": feature_type,
                "Strand": strand,
                **extra.get(int(aid), {}),
            }
            color = None
            for rule in color_rules:
                value = qualifiers.get(rule.get("qualifier_key"))
                if _matches(value, rule):
                    color = _rule_color(value, rule)
                    break
            features.append(
                AnnotationFeatureData(
                    int(aid),
                    int(feature_start),
                    int(feature_end),
                    int(strand or 1),
                    str(feature_type),
                    qualifiers,
                    color,
                    label_key,
                )
            )
        result = GeneMapData(subject_name, subject_length, tuple(features))
        self._record("annotations", started)
        return result

    def nucleotides(self, members, start, end):
        started = time.perf_counter()
        positions, nucleotides = [], []
        for member in members:
            _cid, name, length, offset = member
            lo, hi = max(start, offset + 1), min(end, offset + length)
            if lo > hi:
                continue
            sequence = self.repository.sequence(name, lo - offset, hi - lo + 1)
            positions.extend(range(lo, lo + len(sequence)))
            nucleotides.extend(sequence.upper())
        result = NucleotideSequenceData(tuple(positions), tuple(nucleotides))
        self._record("nucleotides", started)
        return result

    def translation(self, members, start, end):
        started = time.perf_counter()
        entries = []
        for _cid, name, length, offset in members:
            lo, hi = max(start, offset + 1), min(end, offset + length)
            if lo > hi:
                continue
            for cds_start, cds_end, strand, nucleotide, protein, segments in self.repository.cds(
                name, lo - offset, hi - offset
            ):
                entries.append(
                    CdsData(
                        int(cds_start) + offset,
                        int(cds_end) + offset,
                        int(strand or 1),
                        str(nucleotide or ""),
                        str(protein or ""),
                        _normalize_segments(segments, cds_start, cds_end, offset),
                    )
                )
        result = TranslationData(tuple(entries), self.repository.codons())
        self._record("translation", started)
        return result
