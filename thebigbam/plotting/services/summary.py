"""Plain-data construction and transformations for summary reports."""

from __future__ import annotations

from dataclasses import dataclass
from math import floor, log10
from typing import Any, Mapping

from ..repositories.summary import SummaryRepository


@dataclass(frozen=True)
class TableData:
    rows: tuple[tuple[Any, ...], ...]
    columns: tuple[str, ...]
    skip: tuple[str, ...] = ()


@dataclass(frozen=True)
class MetricSection:
    title: str
    labels: Mapping[str, str]
    data: Mapping[str, tuple[Any, ...]]
    note: str = ""


@dataclass(frozen=True)
class SummaryReport:
    subject: str
    subject_table: TableData
    samples: TableData
    metrics: tuple[MetricSection, ...]
    parent_mag: TableData | None = None
    parent_mag_name: str | None = None
    contigs: TableData | None = None


COVERAGE = {
    "Aligned_fraction_percentage": "Aligned fraction (%)",
    "Expected_aligned_fraction": "Expected aligned fraction (%)",
    "Read_count": "Read count",
    "Coverage_mean": "Coverage mean",
    "Coverage_median": "Coverage median",
    "Coverage_trimmed_mean": "Coverage trimmed mean",
    "RPKM": "RPKM",
    "TPM": "TPM",
    "Coverage_coefficient_of_variation": "Coverage coefficient of variation",
    "Relative_coverage_roughness": "Relative coverage roughness",
}
MISASSEMBLY = {
    "Mismatches_per_100kbp": "Mismatches (per 100kbp)",
    "Deletions_per_100kbp": "Deletions (per 100kbp)",
    "Insertions_per_100kbp": "Insertions (per 100kbp)",
    "Clippings_per_100kbp": "Clippings (per 100kbp)",
    "Collapse_bp": "Collapse (bp)",
    "Collapse_per_100kbp": "Collapse (per 100kbp)",
    "Expansion_bp": "Expansion (bp)",
    "Expansion_per_100kbp": "Expansion (per 100kbp)",
}
MICRO = {
    "Micro_Mismatches_per_100kbp": "Mismatches (per 100kbp)",
    "Micro_Deletions_per_100kbp": "Deletions (per 100kbp)",
    "Micro_Insertions_per_100kbp": "Insertions (per 100kbp)",
    "Micro_Clippings_per_100kbp": "Clippings (per 100kbp)",
    "Microdiverse_bp_on_reads": "Microdiverse bp (reads)",
    "Microdiverse_bp_per_100kbp_on_reads": "Microdiverse bp/100kbp (reads)",
    "Microdiverse_bp_on_reference": "Microdiverse bp (reference)",
    "Microdiverse_bp_per_100kbp_on_reference": "Microdiverse bp/100kbp (reference)",
}
SIDE = {
    "Coverage_first_position": "Coverage first position",
    "Contig_start_collapse_prevalence": "Start collapse prevalence",
    "Contig_start_collapse_bp": "Start collapse (bp)",
    "Contig_start_expansion_bp": "Start expansion (bp)",
    "Coverage_last_position": "Coverage last position",
    "Contig_end_collapse_prevalence": "End collapse prevalence",
    "Contig_end_collapse_bp": "End collapse (bp)",
    "Contig_end_expansion_bp": "End expansion (bp)",
    "Contig_end_misjoint_mates": "End misjoint mates",
    "Normalized_contig_end_misjoint_mates": "Normalized end misjoint mates",
}
PHAGE = {
    "Packaging_mechanism": "Mechanism",
    "Repeat_length": "Repeat length",
    "Terminase_distance": "Terminase distance",
    "Terminase_percentage": "Terminase percentage",
    "Left_termini": "Left termini",
    "Median_right_termini_clippings": "Median right termini clippings",
    "Right_termini": "Right termini",
    "Median_left_termini_clippings": "Median left termini clippings",
}
TOPOLOGY = {
    "Circularising_reads": "Circularising reads",
    "Circularising_reads_prevalence": "Circularising reads prevalence",
    "Circularising_inserts": "Circularising inserts",
    "Circularising_insert_size_deviation": "Circularising insert size deviation",
    "Normalized_circularising_inserts": "Normalized circularising inserts",
}
SQL_ALIASES = {key: key.removeprefix("Micro_") for key in MICRO}


def column_scales(connection):
    return {name: float(scale) for name, scale in SummaryRepository(connection).column_scales()}


def decode_map(connection):
    result = {}
    for column, scale in column_scales(connection).items():
        if scale != 1:
            digits = max(0, len(str(int(scale))) - 1)
            result[column] = lambda value, scale=scale, digits=digits: round(value / scale, digits)
    return result


def filter_encode(connection):
    return {column: scale for column, scale in column_scales(connection).items() if scale != 1}


def round_significant(value, figures=2):
    if value is None or not isinstance(value, float):
        return value
    if value == 0:
        return 0
    return round(value, -int(floor(log10(abs(value)))) + figures - 1)


class SummaryDataService:
    def __init__(self, repository: SummaryRepository) -> None:
        self.repository = repository

    def metric_sections(self, kind: str, name: str, samples) -> tuple[MetricSection, ...]:
        samples = tuple(samples)
        subject_id = self.repository.identifier(kind, f"{kind}_id", f"{kind}_name", name)
        sample_ids = self.repository.sample_ids(samples)
        if subject_id is None or not sample_ids:
            return ()
        configs = [
            ("Coverage", COVERAGE, ""),
            ("Misassembly", MISASSEMBLY, ""),
            ("Microdiversity", MICRO, "Positions with ≥10% prevalence."),
        ]
        if kind == "Contig":
            configs += [
                ("Side misassembly", SIDE.copy(), "Clipping events at contig extremities with ≥50% prevalence."),
                ("Topology", TOPOLOGY.copy(), ""),
                ("Phage mechanism", PHAGE, ""),
            ]
        if kind == "Contig" and not self.repository.has_paired_samples(samples):
            for key in ("Contig_end_misjoint_mates", "Normalized_contig_end_misjoint_mates"):
                configs[3][1].pop(key, None)
            for key in (
                "Circularising_inserts",
                "Circularising_insert_size_deviation",
                "Normalized_circularising_inserts",
            ):
                configs[4][1].pop(key, None)
        suffix = "_per_MAG" if kind == "MAG" else ""
        index = {sample: i for i, sample in enumerate(samples)}
        sections = []
        for title, labels, note in configs:
            view_base = "phage_mechanisms" if title == "Phage mechanism" else title.lower().replace(" ", "_")
            view = f"Explicit_{view_base}{suffix}"
            if not labels or not self.repository.view_exists(view):
                continue
            sql_columns = [SQL_ALIASES.get(column, column) for column in labels]
            rows = self.repository.metric_rows(view, f"{kind}_id", subject_id, sample_ids.values(), sql_columns)
            data = {"Sample": list(samples), **{column: [None] * len(samples) for column in labels}}
            for sample, *values in rows:
                if sample not in index:
                    continue
                for column, value in zip(labels, values):
                    data[column][index[sample]] = value
            sections.append(MetricSection(title, labels, {key: tuple(value) for key, value in data.items()}, note))
        return tuple(sections)

    def report(self, contig_name, samples, *, mag_name=None, is_mag_view=False) -> SummaryReport:
        samples = tuple(samples)
        sample_rows, sample_columns = self.repository.named_rows("Sample", "Sample_name", samples)
        if is_mag_view:
            subject = mag_name or contig_name
            row, columns = self.repository.named_row("MAG", "MAG_name", subject)
            contigs, contig_columns = self.repository.contigs_for_mag(subject)
            return SummaryReport(
                subject,
                TableData((row,) if row else (), columns, ("MAG_id",)),
                TableData(tuple(sample_rows), sample_columns, ("Sample_id",)),
                self.metric_sections("MAG", subject, samples),
                contigs=TableData(tuple(contigs), contig_columns, ("Contig_id",)),
            )
        row, columns = self.repository.named_row("Contig", "Contig_name", contig_name)
        mag_row, mag_columns = self.repository.mag_for_contig(contig_name)
        resolved_mag = mag_name
        if mag_row and "MAG_name" in mag_columns:
            resolved_mag = mag_row[mag_columns.index("MAG_name")]
        metrics = list(self.metric_sections("Contig", contig_name, samples))
        if resolved_mag:
            metrics = [*self.metric_sections("MAG", resolved_mag, samples), *metrics]
        return SummaryReport(
            contig_name,
            TableData((row,) if row else (), columns, ("Contig_id",)),
            TableData(tuple(sample_rows), sample_columns, ("Sample_id",)),
            tuple(metrics),
            TableData((mag_row,) if mag_row else (), mag_columns, ("MAG_id",)) if mag_row else None,
            resolved_mag,
        )
