"""Database metadata retrieval for server startup."""

from __future__ import annotations

import duckdb


class DatabaseMetadataRepository:
    def load(self, db_path: str) -> dict[str, str]:
        connection = duckdb.connect(db_path, read_only=True)
        try:
            return dict(connection.execute("SELECT Key, Value FROM Database_metadata").fetchall())
        finally:
            connection.close()
