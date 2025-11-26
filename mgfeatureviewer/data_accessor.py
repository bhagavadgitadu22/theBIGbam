"""Data accessor module for MGFeatureViewer - pure SQLite."""

import sqlite3
from pathlib import Path


class DataAccessor:
    """Unified data accessor for MGFeatureViewer outputs.

    Supports both:
    - Directory format: directory containing metadata.db
    - Single file format: single .db file

    All feature data is stored in Feature_* tables within the SQLite database.
    """

    def __init__(self, db_path: str):
        """Initialize the data accessor.

        Args:
            db_path: Path to either a directory (containing metadata.db) or .db file
        """
        self.db_path = Path(db_path)

        if self.db_path.is_dir():
            # Directory format: look for metadata.db inside
            self.sqlite_path = self.db_path / "metadata.db"
            if not self.sqlite_path.exists():
                raise FileNotFoundError(f"metadata.db not found in {self.db_path}")
        else:
            # Single file format
            self.sqlite_path = self.db_path

        self.sqlite_conn = sqlite3.connect(str(self.sqlite_path))

    def get_sqlite_connection(self) -> sqlite3.Connection:
        """Get the SQLite connection for metadata queries."""
        return self.sqlite_conn

    def get_feature_data(self, feature: str, contig_id: int, sample_id: int,
                         contig_name: str = None, sample_name: str = None):
        """Get feature data for a specific contig and sample.

        Args:
            feature: Feature subplot name (e.g., "Coverage", "Read lengths")
            contig_id: Contig ID from database
            sample_id: Sample ID from database
            contig_name: Unused, kept for API compatibility
            sample_name: Unused, kept for API compatibility

        Returns:
            List of feature dicts with x, y data and rendering info
        """
        cur = self.sqlite_conn.cursor()

        # Get rendering info from Variable table
        cur.execute(
            "SELECT Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name "
            "FROM Variable WHERE Subplot=?",
            (feature,)
        )
        rows = cur.fetchall()

        list_feature_dict = []
        for row in rows:
            type_picked, color, alpha, fill_alpha, size, title, feature_table = row
            feature_dict = {
                "type": type_picked,
                "color": color,
                "alpha": alpha,
                "fill_alpha": fill_alpha,
                "size": size,
                "title": title,
                "x": [],
                "y": []
            }

            # Query Feature_* table directly
            cur.execute(
                f"SELECT Position, Value FROM {feature_table} "
                "WHERE Sample_id=? AND Contig_id=? ORDER BY Position",
                (sample_id, contig_id)
            )
            data_rows = cur.fetchall()
            feature_dict["x"] = [r[0] for r in data_rows]
            feature_dict["y"] = [r[1] for r in data_rows]

            list_feature_dict.append(feature_dict)

        return list_feature_dict

    def close(self):
        """Close database connections."""
        if self.sqlite_conn:
            self.sqlite_conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
