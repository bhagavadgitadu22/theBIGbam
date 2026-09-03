"""Command-line configuration and startup for the plotting server."""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any, Mapping

import panel as pn

from ..repositories.database_metadata import DatabaseMetadataRepository
from ..settings.history import HISTORY_FORMAT, HistoryEntry, entries_from_session_document
from ..settings.persistence import load_settings_document
from ..settings.scenario import ScenarioPathAllocator
from ..shared.paths import static_directory
from ..shared.timing import TimingPhase
from . import composition_root as application
from .apply_handlers import warm_plot_pipeline_imports


def restore_payload(
    document: dict[str, Any] | None,
) -> tuple[dict[str, Any] | None, tuple[HistoryEntry, ...]]:
    """Split either supported --json document into settings and history."""
    if document is None:
        return None, ()
    meta = document.get("_meta")
    if not isinstance(meta, Mapping) or meta.get("format") != HISTORY_FORMAT:
        return document, ()
    entries = entries_from_session_document(document)
    return (entries[-1].settings if entries else None), entries


def add_serve_args(parser) -> None:
    parser.add_argument("--db", required=True, help="Path to DuckDB database")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Panel app")
    parser.add_argument(
        "--json",
        dest="settings_json",
        default=None,
        help="Path to a JSON file from SAVE SETTINGS or SAVE SESSION to restore on load. "
        "A session restores its complete application history and the latest entry's settings. "
        "Settings that don't fit this --db are skipped with a logged warning.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open a browser automatically (useful for headless benchmarks and remote servers)",
    )
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
        help="(for developers) Print timing and memory diagnostics to the terminal",
    )
    parser.add_argument(
        "--scenario",
        default=None,
        help="(for developers) Path where user actions and settings changes are continuously recorded as a "
        "scenario JSON file. The file is kept valid while the server runs, including if the job is terminated.",
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
    initial_document = load_settings_document(settings_path) if settings_path else None
    initial_settings, initial_history = restore_payload(initial_document)
    scenario_path = getattr(args, "scenario", None)
    scenario_paths = None
    if scenario_path:
        requested_path = Path(scenario_path)
        requested_path.parent.mkdir(parents=True, exist_ok=True)
        scenario_paths = ScenarioPathAllocator(requested_path)

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
        allocated_scenario_path = scenario_paths.next_path() if scenario_paths is not None else None
        return application.create_layout(
            args.db,
            preloaded,
            enable_timing=enable_timing,
            initial_settings=initial_settings,
            initial_history=initial_history,
            scenario_path=allocated_scenario_path,
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
        show=not getattr(args, "no_browser", False),
        title="theBIGbam",
        static_dirs={"assets": static_directory()},
        session_token_expiration=3600,
    )
    return 0
