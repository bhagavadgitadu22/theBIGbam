"""Filesystem paths for packaged plotting assets."""

from pathlib import Path

STATIC_DIR = Path(__file__).resolve().parents[2] / "static"


def static_directory() -> str:
    if not STATIC_DIR.is_dir():
        raise FileNotFoundError(f"Packaged static asset directory is missing: {STATIC_DIR}")
    return str(STATIC_DIR)
