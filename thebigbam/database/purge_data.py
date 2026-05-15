"""Purge mapping-derived data from a theBIGbam database.

Produces a smaller DB that retains genome-level data (annotations, sequences,
BLAST repeats, GC content) and is ready to be reprocessed with
``thebigbam calculate --extend``.
"""

import os
import tempfile

import duckdb

from thebigbam.database.database_getters import _table_exists, update_database_metadata


COMPRESSION_METADATA_KEYS = [
    "Min_aligned_fraction",
    "Min_coverage_depth",
    "Coverage_percentage",
    "Variation_percentage",
    "Min_occurrences",
]

MAPPING_TABLES = {
    "Sample", "Coverage", "Feature_blob", "Feature_blob_chunk",
    "Phage_mechanisms", "Phage_termini", "Misassembly", "Microdiversity",
    "Side_misassembly", "Topology", "MAG_blob", "MAG_coverage",
}


def purge_mapping_data(input_db, output_db):
    if not os.path.exists(input_db):
        raise FileNotFoundError(f"Input database not found: {input_db}")
    if os.path.abspath(input_db) == os.path.abspath(output_db):
        raise ValueError("Input and output paths must differ")
    if os.path.exists(output_db):
        raise FileExistsError(f"Output path already exists: {output_db}")

    src = duckdb.connect(input_db, read_only=True)
    deleted = {}
    for table in MAPPING_TABLES:
        if _table_exists(src, table):
            count = src.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
            if count > 0:
                deleted[table] = count

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Exporting {input_db} ...")
        src.execute(f"EXPORT DATABASE '{tmpdir}'")
        src.close()

        schema_path = os.path.join(tmpdir, "schema.sql")
        load_path = os.path.join(tmpdir, "load.sql")

        # Remove recreatable views from schema and strip mapping-table data from load
        with open(schema_path) as f:
            schema_lines = f.readlines()
        with open(schema_path, "w") as f:
            for line in schema_lines:
                stripped = line.lstrip()
                if stripped.startswith("CREATE VIEW Explicit_"):
                    continue
                if stripped.startswith("CREATE VIEW Contig_blast_hits_symmetric"):
                    continue
                f.write(line)

        mapping_lower = {t.lower() for t in MAPPING_TABLES}
        with open(load_path) as f:
            load_lines = f.readlines()
        with open(load_path, "w") as f:
            for line in load_lines:
                stripped = line.strip()
                # COPY lines: COPY TableName FROM ... or COPY "TableName" FROM ...
                if stripped.startswith("COPY "):
                    token = stripped.split()[1].strip('"')
                    if token.lower() in mapping_lower:
                        continue
                f.write(line)

        print(f"Building {output_db} ...")
        dst = duckdb.connect(output_db)
        dst.execute(f"IMPORT DATABASE '{tmpdir}'")

        if _table_exists(dst, "Contig"):
            dst.execute("UPDATE Contig SET Number_of_samples = NULL")

        if _table_exists(dst, "Database_metadata"):
            for key in COMPRESSION_METADATA_KEYS:
                dst.execute("DELETE FROM Database_metadata WHERE Key = ?", [key])

        update_database_metadata(dst)
        dst.commit()
        dst.close()

    orig_size = os.path.getsize(input_db) / (1024 * 1024)
    new_size = os.path.getsize(output_db) / (1024 * 1024)
    print(f"\nPurged mapping data:")
    for table, count in deleted.items():
        if count > 0:
            print(f"  {table}: {count} rows removed")
    print(f"\nSize: {orig_size:.1f} MB -> {new_size:.1f} MB ({orig_size - new_size:.1f} MB saved)")
    print(f"\nReady for: thebigbam calculate --extend -b <bams> -o {output_db}")


def add_purge_args(parser):
    parser.add_argument('-d', '--db', required=True, help='Input database to purge')
    parser.add_argument('-o', '--output', required=True, help='Output path for purged database')


def run_purge(args):
    try:
        purge_mapping_data(args.db, args.output)
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 2
