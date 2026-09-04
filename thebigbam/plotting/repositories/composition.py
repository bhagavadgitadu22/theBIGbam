"""Entity validation and ordering queries used by typed plot composers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ...database.database_getters import (
    get_mag_contigs,
    get_mag_id,
    get_mag_members_full,
    get_mag_members_full_sorted,
)


@dataclass(frozen=True)
class ContigContext:
    contig_id: int
    name: str
    length: int


@dataclass(frozen=True)
class MagContext:
    mag_id: int
    members: tuple[tuple[str, int, int], ...]
    feature_members: tuple[tuple[int, int, int], ...]
    total_length: int
    xstart: int | None
    xend: int | None

    @property
    def genome_members(self) -> tuple[tuple[int, str, int, int], ...]:
        """Combine the two aligned repository projections without another query."""
        names_by_offset = {int(offset): str(name) for name, _length, offset in self.members}
        return tuple(
            (int(contig_id), names_by_offset[int(offset)], int(length), int(offset))
            for contig_id, length, offset in self.feature_members
            if int(offset) in names_by_offset
        )


class CompositionRepository:
    def __init__(self, connection: Any) -> None:
        self.connection = connection

    def contig(self, name: str) -> ContigContext:
        row = self.connection.execute(
            "SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name=?", [name]
        ).fetchone()
        if row is None:
            raise ValueError(f"Contig not found: {name}")
        return ContigContext(int(row[0]), row[1], int(row[2]))

    def sample(self, name: str | None) -> tuple[int | None, str | None]:
        if not name:
            return None, name
        row = self.connection.execute(
            "SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name=?", [name]
        ).fetchone()
        if row is None:
            raise ValueError(f"Sample not found: {name}")
        return int(row[0]), row[1]

    def mag_contigs(self, name: str) -> tuple[tuple[str, int, int], ...]:
        """Return ordered MAG members as plain plotting-domain values."""
        return tuple(get_mag_contigs(self.connection, name))

    def mag(self, request) -> MagContext:
        mag_id = get_mag_id(self.connection, request.mag_name)
        if mag_id is None:
            raise ValueError(f"MAG_id lookup failed for: {request.mag_name}")
        ordering = request.ordering
        if ordering.source is not None and ordering.metric:
            full = get_mag_members_full_sorted(
                self.connection,
                mag_id,
                ordering.source,
                ordering.metric,
                ordering.ascending,
                ordering.sample_name,
            )
        else:
            full = get_mag_members_full(self.connection, mag_id)
        members = [(name, length, offset) for _cid, name, length, offset in full]
        feature_members = [(cid, length, offset) for cid, _name, length, offset in full]
        if not members:
            raise ValueError(f"MAG not found or empty: {request.mag_name}")
        xstart, xend = request.xstart, request.xend
        if request.focus_contig is not None:
            focused = next(
                ((length, offset) for name, length, offset in members if name == request.focus_contig),
                None,
            )
            if focused is not None:
                length, offset = focused
                xstart, xend = offset + 1, offset + length
        return MagContext(
            int(mag_id),
            tuple(members),
            tuple(feature_members),
            sum(length for _name, length, _offset in members),
            xstart,
            xend,
        )

    def mag_samples(self, mag_id: int) -> list[tuple[int, str]]:
        rows = self.connection.execute(
            "SELECT DISTINCT s.Sample_id, s.Sample_name FROM Sample s "
            "JOIN Feature_blob fb ON fb.Sample_id = s.Sample_id "
            "JOIN MAG_contigs_association mca ON mca.Contig_id = fb.Contig_id "
            "WHERE mca.MAG_id = ?",
            [mag_id],
        ).fetchall()
        return [(int(sid), name) for sid, name in rows]

    def mag_sample_order_values(self, rows, mag_name, ordering):
        """Return immutable sample ordering values for in-memory presentation sorting."""
        if ordering.source == "Sample":
            ids = [sid for sid, _name in rows]
            if not ids:
                return {}
            placeholders = ",".join("?" for _ in ids)
            values = self.connection.execute(
                f'SELECT Sample_id, "{ordering.metric}" FROM Sample '
                f"WHERE Sample_id IN ({placeholders})",
                ids,
            ).fetchall()
            return {int(sid): value for sid, value in values}
        names = [name for _sid, name in rows]
        if not names:
            return {}
        placeholders = ",".join("?" for _ in names)
        values = self.connection.execute(
            f'SELECT Sample_name, "{ordering.metric}" FROM "{ordering.source}" '
            f"WHERE MAG_name = ? AND Sample_name IN ({placeholders})",
            [mag_name, *names],
        ).fetchall()
        id_by_name = {name: sid for sid, name in rows}
        return {int(id_by_name[name]): value for name, value in values if name in id_by_name}

    def order_mag_samples(self, rows, mag_name, ordering):
        direction = "ASC" if ordering.ascending else "DESC"
        if ordering.source == "Sample":
            ids = [sid for sid, _name in rows]
            placeholders = ",".join("?" for _ in ids)
            ordered = self.connection.execute(
                f"SELECT Sample_id FROM Sample WHERE Sample_id IN ({placeholders}) "
                f'ORDER BY "{ordering.metric}" {direction} NULLS LAST',
                ids,
            ).fetchall()
            by_id = {sid: name for sid, name in rows}
            return [(sid, by_id[sid]) for (sid,) in ordered]
        names = [name for _sid, name in rows]
        placeholders = ",".join("?" for _ in names)
        ordered_names = [
            row[0]
            for row in self.connection.execute(
                f'SELECT Sample_name FROM "{ordering.source}" WHERE MAG_name = ? '
                f"AND Sample_name IN ({placeholders}) "
                f'ORDER BY "{ordering.metric}" {direction} NULLS LAST',
                [mag_name, *names],
            ).fetchall()
        ]
        ordered_names.extend(name for name in names if name not in set(ordered_names))
        by_name = {name: (sid, name) for sid, name in rows}
        return [by_name[name] for name in ordered_names]
