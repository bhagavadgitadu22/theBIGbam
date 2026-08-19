"""Command-line configuration and startup for the plotting server."""

from __future__ import annotations

import os
import time

import panel as pn

from ..repositories.database_metadata import DatabaseMetadataRepository
from ..settings.persistence import load_settings_document
from ..shared.paths import static_directory
from ..shared.timing import TimingPhase
from . import composition_root as application
from .apply_render_handlers import warm_plot_pipeline_imports


def add_serve_args(parser) -> None:
    parser.add_argument("--db", required=True, help="Path to DuckDB database")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Panel app")
    parser.add_argument(
        "--allow-websocket-origin",
        action="append",
        default=None,
        help="Additional browser origin allowed to connect (for example host.example.org:5006). May be repeated.",
    )
    parser.add_argument(
        "--time",
        action="store_true",
        default=False,
        help="Print timing and memory diagnostics to the terminal",
    )
    parser.add_argument(
        "--json",
        dest="settings_json",
        default=None,
        help="Path to a settings JSON file (from SAVE SETTINGS) to restore on load. "
        "Settings that don't fit this --db are skipped with a logged warning.",
    )


def _print_database_metadata(db_path: str) -> None:
    metadata = DatabaseMetadataRepository().load(db_path)

    keys = (
        "Modules",
        "View_mode",
        "Min_aligned_fraction",
        "Min_coverage_depth",
        "Coverage_percentage",
        "Min_occurrences",
        "Variation_percentage",
    )
    parameters = [f"{key}={metadata[key]}" for key in keys if key in metadata]
    print(
        f"Database '{os.path.basename(db_path)}': "
        f"\n Created on {metadata.get('Date_of_creation', '?')} "
        f"(v{metadata.get('Tool_version_used_for_creation', '?')}) "
        f"\n Last modified on {metadata.get('Date_of_last_modification', '?')} "
        f"(v{metadata.get('Tool_version_used_for_last_modification', '?')})",
        flush=True,
    )
    if parameters:
        parameter_text = "\n ".join(parameters)
        print(f"Calculate parameters used:\n {parameter_text}", flush=True)


def run_serve(args) -> int:
    enable_timing = getattr(args, "time", False)
    application._TIMING = TimingPhase()
    if enable_timing:
        application._start_rss_watchdog()
    _print_database_metadata(args.db)

    settings_path = getattr(args, "settings_json", None)
    initial_settings = load_settings_document(settings_path) if settings_path else None

    print("\nPreloading database data...", flush=True)
    preloaded = application.preload_db_data(args.db, enable_timing=enable_timing)
    warm_started = time.perf_counter()
    warm_plot_pipeline_imports()
    if enable_timing:
        elapsed = time.perf_counter() - warm_started
        print(f"[timing] Composer import warm-up: {elapsed:.3f}s{application._TIMING.tag(elapsed)}", flush=True)
    if enable_timing:
        application._TIMING.summary("Startup")
    print(f"Server ready. Open localhost:{args.port} in your browser.\n", flush=True)

    def create_app():
        return application.create_layout(
            args.db,
            preloaded,
            enable_timing=enable_timing,
            initial_settings=initial_settings,
        )

    allowed_origins = getattr(args, "allow_websocket_origin", None) or [
        f"localhost:{args.port}",
        f"127.0.0.1:{args.port}",
    ]
    pn.serve(
        create_app,
        port=args.port,
        address="0.0.0.0",
        allow_websocket_origin=allowed_origins,
        show=True,
        title="theBIGbam",
        static_dirs={"assets": static_directory()},
        session_token_expiration=3600,
    )
    return 0
