"""Read-only startup retrieval for the plotting application."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterator, Mapping

import duckdb

from ...database.database_getters import get_filtering_metadata, get_mag_contig_map, is_mag_mode
from ..services.summary import filter_encode


@dataclass(frozen=True)
class PreloadedPlotData(Mapping[str, Any]):
    annotation_types: tuple[str, ...]
    has_mags: bool
    mag_to_contigs: Mapping[str, list[str]]
    contig_to_mag: Mapping[str, str]
    mags: tuple[str, ...]
    contigs: tuple[str, ...]
    contig_lengths: Mapping[str, int]
    contig_name_to_id: Mapping[str, int]
    contig_id_to_name: Mapping[int, str]
    sample_name_to_id: Mapping[str, int]
    sample_id_to_name: Mapping[int, str]
    samples: tuple[str, ...]
    has_samples: bool
    mag_to_sample_ids: Mapping[str, set[int]]
    module_names: tuple[str, ...]
    module_variables: tuple[tuple[str, ...], ...]
    module_helps: tuple[str, ...]
    custom_contig_subplots: tuple[str, ...]
    filtering_metadata: Mapping[str, Any]
    subplot_to_varnames: Mapping[str, list[str]]
    encoding_by_feature: Mapping[str, str]
    total_coverage_count: int
    filter_encode: Mapping[str, Any]
    mag_to_contig_offsets: dict[str, dict[str, int]] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        try:
            return getattr(self, key)
        except AttributeError as error:
            raise KeyError(key) from error

    def __iter__(self) -> Iterator[str]:
        return iter(self.__dataclass_fields__)

    def __len__(self) -> int:
        return len(self.__dataclass_fields__)


class PreloadRepository:
    """Load immutable application metadata using a short-lived read-only connection."""

    _MODULE_ORDER = ("Genome", "Coverage", "Misalignment", "RNA", "Long-reads", "Paired-reads", "Phage termini")

    def __init__(
        self,
        db_path: str,
        *,
        scale_initializer=None,
        filter_encoder=filter_encode,
        filtering_metadata_loader=get_filtering_metadata,
        mag_mode_loader=is_mag_mode,
        mag_map_loader=get_mag_contig_map,
    ) -> None:
        self.db_path = db_path
        self.query_count = 0
        self.scale_initializer = scale_initializer
        self.filter_encoder = filter_encoder
        self.filtering_metadata_loader = filtering_metadata_loader
        self.mag_mode_loader = mag_mode_loader
        self.mag_map_loader = mag_map_loader

    def _execute(self, cursor, sql, parameters=()):
        self.query_count += 1
        return cursor.execute(sql, parameters)

    def load(self) -> PreloadedPlotData:
        connection = duckdb.connect(self.db_path, read_only=True)
        try:
            initializer = self.scale_initializer
            if initializer is None:
                from ...database.blob_decoder import init_scales

                initializer = init_scales
            initializer(connection)
            cursor = connection.cursor()
            filter_encode = self.filter_encoder(connection)

            exists = self._execute(
                cursor, "SELECT 1 FROM information_schema.tables WHERE table_name = 'Annotated_types'"
            ).fetchone()
            annotation_types = (
                tuple(
                    row[0]
                    for row in self._execute(
                        cursor, "SELECT Type_name FROM Annotated_types ORDER BY Frequency DESC"
                    ).fetchall()
                )
                if exists
                else ()
            )

            has_mags = self.mag_mode_loader(connection)
            mag_to_contigs, contig_to_mag = self.mag_map_loader(connection)
            mags = tuple(sorted(mag_to_contigs))

            contig_rows = self._execute(
                cursor, "SELECT Contig_id, Contig_name, Contig_length FROM Contig ORDER BY Contig_name"
            ).fetchall()
            contigs = tuple(row[1] for row in contig_rows)
            contig_lengths = {row[1]: row[2] for row in contig_rows}
            contig_name_to_id = {row[1]: row[0] for row in contig_rows}
            contig_id_to_name = {row[0]: row[1] for row in contig_rows}

            has_sample_table = (
                self._execute(cursor, "SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone()
                is not None
            )
            sample_rows = (
                self._execute(cursor, "SELECT Sample_id, Sample_name FROM Sample ORDER BY Sample_name").fetchall()
                if has_sample_table
                else []
            )
            samples = tuple(row[1] for row in sample_rows)
            sample_name_to_id = {row[1]: row[0] for row in sample_rows}
            sample_id_to_name = {row[0]: row[1] for row in sample_rows}

            mag_to_sample_ids: dict[str, set[int]] = {}
            if has_sample_table and mag_to_contigs:
                has_mag_coverage = (
                    self._execute(
                        cursor, "SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG_coverage'"
                    ).fetchone()
                    is not None
                )
                sql = (
                    "SELECT mg.MAG_name, mc.Sample_id FROM MAG mg JOIN MAG_coverage mc ON mc.MAG_id = mg.MAG_id"
                    if has_mag_coverage
                    else "SELECT DISTINCT mg.MAG_name, p.Sample_id FROM MAG mg "
                    "JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id "
                    "JOIN Coverage p ON p.Contig_id = mca.Contig_id"
                )
                for mag_name, sample_id in self._execute(cursor, sql).fetchall():
                    mag_to_sample_ids.setdefault(mag_name, set()).add(sample_id)

            tables = tuple(
                row[0] for row in self._execute(cursor, "SELECT DISTINCT Feature_table_name FROM Variable").fetchall()
            )
            placeholders = ",".join("?" * len(tables))
            modules = []
            if tables:
                modules = [
                    row[0]
                    for row in self._execute(
                        cursor,
                        f"SELECT DISTINCT Module FROM Variable WHERE Feature_table_name IN ({placeholders})",
                        tables,
                    ).fetchall()
                ]
            modules.sort(
                key=lambda name: (
                    self._MODULE_ORDER.index(name) if name in self._MODULE_ORDER else len(self._MODULE_ORDER)
                )
            )
            custom_contig_subplots = tuple(
                row[0]
                for row in self._execute(
                    cursor, "SELECT Subplot FROM Variable WHERE Module='Custom' AND Feature_table_name LIKE 'Contig_%'"
                ).fetchall()
            )
            module_names, module_variables, module_helps = [], [], []
            for module in modules:
                variables = [
                    row[0]
                    for row in self._execute(
                        cursor,
                        f"SELECT DISTINCT Subplot FROM Variable WHERE Module=? "
                        f"AND Feature_table_name IN ({placeholders}) ORDER BY Module_order",
                        (module, *tables),
                    ).fetchall()
                ]
                if module == "Custom":
                    variables = [value for value in variables if value not in custom_contig_subplots]
                if not variables:
                    continue
                module_names.append(module)
                module_variables.append(tuple(variables))
                help_rows = self._execute(
                    cursor,
                    f"SELECT DISTINCT Subplot, Title, Help FROM Variable WHERE Module=? "
                    f"AND Feature_table_name IN ({placeholders}) ORDER BY Subplot",
                    (module, *tables),
                ).fetchall()
                module_helps.append(
                    "".join(
                        f"{title} ({subplot} subplot): {help_text}\n"
                        for subplot, title, help_text in help_rows
                        if help_text and help_text.strip()
                    )
                )

            subplot_to_varnames: dict[str, list[str]] = {}
            for variable, subplot in self._execute(
                cursor,
                "SELECT Variable_name, Subplot FROM Variable WHERE Variable_name IS NOT NULL AND Subplot IS NOT NULL",
            ).fetchall():
                subplot_to_varnames.setdefault(subplot, []).append(variable)
            encoding_by_feature = dict(
                self._execute(
                    cursor,
                    "SELECT Variable_name, Encoding FROM Variable "
                    "WHERE Variable_name IS NOT NULL AND Encoding IS NOT NULL",
                ).fetchall()
            )
            total_coverage_count = self._execute(cursor, "SELECT COUNT(*) FROM Coverage").fetchone()[0]
            return PreloadedPlotData(
                annotation_types=annotation_types,
                has_mags=has_mags,
                mag_to_contigs=mag_to_contigs,
                contig_to_mag=contig_to_mag,
                mags=mags,
                contigs=contigs,
                contig_lengths=contig_lengths,
                contig_name_to_id=contig_name_to_id,
                contig_id_to_name=contig_id_to_name,
                sample_name_to_id=sample_name_to_id,
                sample_id_to_name=sample_id_to_name,
                samples=samples,
                has_samples=bool(samples),
                mag_to_sample_ids=mag_to_sample_ids,
                module_names=tuple(module_names),
                module_variables=tuple(module_variables),
                module_helps=tuple(module_helps),
                custom_contig_subplots=custom_contig_subplots,
                filtering_metadata=self.filtering_metadata_loader(self.db_path),
                subplot_to_varnames=subplot_to_varnames,
                encoding_by_feature=encoding_by_feature,
                total_coverage_count=total_coverage_count,
                filter_encode=filter_encode,
            )
        finally:
            connection.close()
