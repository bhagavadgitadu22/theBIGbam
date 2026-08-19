"""DuckDB retrieval for genome annotation and sequence tracks."""

from __future__ import annotations


class GenomeTrackRepository:
    def __init__(self, conn) -> None:
        self.conn = conn
        self.query_count = 0

    def _execute(self, sql, params=()):
        self.query_count += 1
        return self.conn.execute(sql, params)

    def contig_identity(self, name):
        row = self._execute("SELECT Contig_id, Contig_length FROM Contig WHERE Contig_name = ?", [name]).fetchone()
        if row is None:
            raise ValueError(f"Unknown contig: {name}")
        return int(row[0]), int(row[1])

    def contig_identities(self, names):
        if not names:
            return {}
        marks = ",".join("?" for _ in names)
        rows = self._execute(
            f"SELECT Contig_name, Contig_id, Contig_length FROM Contig WHERE Contig_name IN ({marks})",
            list(names),
        ).fetchall()
        return {str(name): (int(contig_id), int(length)) for name, contig_id, length in rows}

    def annotations(self, contig_ids, offsets, start, end, feature_types, plot_isoforms):
        placeholders = ",".join("?" for _ in contig_ids)
        params = list(contig_ids)
        type_sql = ""
        if feature_types:
            type_sql = f' AND ca."Type" IN ({",".join("?" for _ in feature_types)})'
            params.extend(feature_types)
        iso_sql = "" if plot_isoforms else " AND (lq.Value IS NULL OR ca.Main_isoform = true)"
        rows = self._execute(
            'SELECT ca.Annotation_id, ca.Contig_id, ca."Start", ca."End", ca.Strand, ca."Type", '
            "ca.Main_isoform, pq.Value, fq.Value, lq.Value FROM Contig_annotation_core ca "
            "LEFT JOIN Annotation_qualifier pq ON pq.Annotation_id=ca.Annotation_id AND pq.\"Key\"='product' "
            "LEFT JOIN Annotation_qualifier fq ON fq.Annotation_id=ca.Annotation_id AND fq.\"Key\"='function' "
            "LEFT JOIN Annotation_qualifier lq ON lq.Annotation_id=ca.Annotation_id AND lq.\"Key\"='locus_tag' "
            f"WHERE ca.Contig_id IN ({placeholders}){type_sql}{iso_sql}",
            params,
        ).fetchall()
        shifted = []
        for row in rows:
            offset = offsets.get(int(row[1]), 0)
            shifted_row = (row[0], int(row[2]) + offset, int(row[3]) + offset, *row[4:])
            if shifted_row[2] >= start and shifted_row[1] <= end:
                shifted.append(shifted_row)
        return shifted

    def qualifiers(self, annotation_ids, keys):
        if not annotation_ids or not keys:
            return {}
        aid_marks = ",".join("?" for _ in annotation_ids)
        key_marks = ",".join("?" for _ in keys)
        rows = self._execute(
            f'SELECT Annotation_id, "Key", "Value" FROM Annotation_qualifier '
            f'WHERE Annotation_id IN ({aid_marks}) AND "Key" IN ({key_marks})',
            [*annotation_ids, *keys],
        ).fetchall()
        result = {}
        for aid, key, value in rows:
            result.setdefault(int(aid), {})[str(key)] = value
        return result

    def sequence(self, contig_name, local_start, length):
        row = self._execute(
            "SELECT SUBSTR(cs.Sequence, ?, ?) FROM Contig_sequence cs "
            "JOIN Contig c ON c.Contig_id=cs.Contig_id WHERE c.Contig_name=?",
            [local_start, length, contig_name],
        ).fetchone()
        return str(row[0]) if row and row[0] else ""

    def cds(self, contig_name, start, end):
        exists = self._execute(
            "SELECT 1 FROM information_schema.columns WHERE table_name='Contig_annotation' "
            "AND column_name='Protein_sequence'"
        ).fetchone()
        if exists is None:
            return []
        return self._execute(
            'SELECT ca.Start, ca."End", ca.Strand, ca.Nucleotide_sequence, ca.Protein_sequence, ca.Segments '
            "FROM Contig_annotation ca JOIN Contig c ON c.Contig_id=ca.Contig_id "
            "LEFT JOIN Annotation_qualifier lq ON lq.Annotation_id=ca.Annotation_id AND lq.\"Key\"='locus_tag' "
            "WHERE c.Contig_name=? AND ca.Type='CDS' AND ca.Protein_sequence IS NOT NULL "
            'AND ca."End">=? AND ca.Start<=? AND (lq.Value IS NULL OR ca.Main_isoform=true)',
            [contig_name, start, end],
        ).fetchall()

    def codons(self):
        rows = self._execute(
            "SELECT Codon, AminoAcid, AminoAcid_name, AminoAcid_label, Color FROM Codon_table"
        ).fetchall()
        return {str(c).upper(): (aa, name, label, color) for c, aa, name, label, color in rows}
