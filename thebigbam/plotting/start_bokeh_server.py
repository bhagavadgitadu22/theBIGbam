"""CLI-compatible entry point for the Panel/Bokeh plotting server."""

from __future__ import annotations

import argparse

from .application import PlotApplication, create_application, create_layout, preload_db_data
from .application.server import add_serve_args, run_serve
from .controls.base import build_controls

__all__ = [
    "PlotApplication",
    "add_serve_args",
    "build_controls",
    "create_application",
    "create_layout",
    "preload_db_data",
    "run_serve",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    add_serve_args(parser)
    return run_serve(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
