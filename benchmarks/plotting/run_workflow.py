#!/usr/bin/env python3
"""Drive the canonical theBIGbam plotting workflow and record per-step timings."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable

from playwright.sync_api import Page, sync_playwright


@dataclass
class StepResult:
    sequence: int
    action: str
    status: str
    duration_seconds: float
    server_lines: list[str] = field(default_factory=list)
    error: str = ""
    artifacts: dict[str, str] = field(default_factory=dict)


class Server:
    def __init__(self, db: Path, port: int, output: Path, settings: Path) -> None:
        self.db = db
        self.port = port
        self.output = output
        self.settings = settings
        self.process: subprocess.Popen[str] | None = None
        self.lines: list[str] = []
        self._lock = threading.Lock()
        self._reader: threading.Thread | None = None

    def start(self, timeout: float) -> None:
        executable = shutil.which("thebigbam")
        if not executable:
            raise RuntimeError("thebigbam executable is not available in this environment")
        command = [
            executable,
            "serve",
            "--db",
            str(self.db),
            "--port",
            str(self.port),
            "--allow-websocket-origin",
            f"127.0.0.1:{self.port}",
            "--time",
            "--no-browser",
            "--json",
            str(self.settings),
        ]
        self.output.mkdir(parents=True, exist_ok=True)
        log_handle = (self.output / "server.log").open("w", encoding="utf-8")
        self.process = subprocess.Popen(
            command,
            cwd=Path(__file__).resolve().parents[2],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )

        def read_output() -> None:
            assert self.process is not None and self.process.stdout is not None
            for line in self.process.stdout:
                log_handle.write(line)
                log_handle.flush()
                with self._lock:
                    self.lines.append(line.rstrip())
            log_handle.close()

        self._reader = threading.Thread(target=read_output, daemon=True)
        self._reader.start()
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.process.poll() is not None:
                raise RuntimeError("thebigbam serve exited before becoming ready; inspect server.log")
            with self._lock:
                if any("Server ready" in line for line in self.lines):
                    return
            time.sleep(0.2)
        raise TimeoutError("Timed out waiting for theBIGbam server startup")

    def mark(self) -> int:
        with self._lock:
            return len(self.lines)

    def since(self, mark: int) -> list[str]:
        with self._lock:
            return list(self.lines[mark:])

    def stop(self) -> None:
        if self.process is None or self.process.poll() is not None:
            return
        self.process.terminate()
        try:
            self.process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait(timeout=5)


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def set_model(page: Page, name: str, attribute: str, value: Any) -> None:
    page.wait_for_function(
        "name => window.Bokeh && Bokeh.documents.length && Bokeh.documents[0].get_model_by_name(name)",
        name,
    )
    page.evaluate(
        """([name, attribute, value]) => {
            const model = Bokeh.documents[0].get_model_by_name(name);
            model.setv({[attribute]: value});
        }""",
        [name, attribute, value],
    )


def get_model(page: Page, name: str, attribute: str) -> Any:
    return page.evaluate(
        """([name, attribute]) => {
            const model = Bokeh.documents[0].get_model_by_name(name);
            return model ? model[attribute] : null;
        }""",
        [name, attribute],
    )


def click_css(page: Page, selector: str) -> None:
    page.locator(selector).locator("button").first.click()


def wait_sidebar(page: Page, timeout_ms: int) -> None:
    # Let the browser send the event, then wait until any serialized transition
    # has completed. The final quiet window also covers very fast operations for
    # which the busy class appeared and disappeared between polls.
    page.wait_for_timeout(100)
    page.wait_for_function(
        "!document.querySelector('.sidebar-busy')",
        timeout=timeout_ms,
    )
    page.wait_for_timeout(250)


class ScenarioError(ValueError):
    """Raised when a scenario cannot be replayed deterministically."""


def _validate_changes(sequence: int, changes: dict[str, Any]) -> None:
    allowed = {"filtering", "view_mode", "selection", "variables", "contig", "plotting_params"}
    for section in changes:
        if section not in allowed:
            raise ScenarioError(f"Step {sequence} changes unsupported settings section {section!r}")


def load_scenario(path: Path, db: Path | None = None) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ScenarioError(f"Cannot read scenario: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ScenarioError(f"Invalid scenario JSON at line {exc.lineno}, column {exc.colno}") from exc
    if not isinstance(document, dict):
        raise ScenarioError("Scenario root must be a JSON object")
    metadata = document.get("_meta")
    if not isinstance(metadata, dict) or metadata.get("format") != "thebigbam-scenario":
        raise ScenarioError("Not a thebigbam-scenario document")
    if metadata.get("version") != 1:
        raise ScenarioError(f"Unsupported scenario version: {metadata.get('version')!r}")
    if not isinstance(document.get("initial_state"), dict):
        raise ScenarioError("Scenario initial_state must be a JSON object")
    steps = document.get("steps")
    if not isinstance(steps, list) or not steps:
        raise ScenarioError("Scenario must contain at least one step")
    sequences: set[int] = set()
    previous_sequence = 0
    supported = {"state_change", "filter_lookup", "apply_filters", "apply_plot"}
    for position, step in enumerate(steps, 1):
        if not isinstance(step, dict):
            raise ScenarioError(f"Step {position} must be a JSON object")
        sequence = step.get("sequence")
        if not isinstance(sequence, int) or isinstance(sequence, bool) or sequence < 1:
            raise ScenarioError(f"Step {position} has invalid sequence")
        if sequence in sequences:
            raise ScenarioError(f"Duplicate step sequence {sequence}")
        if sequence <= previous_sequence:
            raise ScenarioError(f"Step sequence {sequence} is not in increasing order")
        sequences.add(sequence)
        previous_sequence = sequence
        action = step.get("action")
        if action not in supported:
            raise ScenarioError(f"Step {sequence} has unsupported action {action!r}")
        if action == "state_change" and not isinstance(step.get("changes"), dict):
            raise ScenarioError(f"Step {sequence} state_change requires an object-valued changes field")
        if action == "state_change":
            _validate_changes(sequence, step["changes"])
        if action == "filter_lookup" and not isinstance(step.get("details"), dict):
            raise ScenarioError(f"Step {sequence} filter_lookup requires an object-valued details field")
        if action == "filter_lookup" and "occurrence" in step["details"]:
            occurrence = step["details"]["occurrence"]
            if not isinstance(occurrence, int) or isinstance(occurrence, bool) or occurrence < 1:
                raise ScenarioError(f"Step {sequence} filter_lookup occurrence must be a positive integer")
    source_db = metadata.get("source_db")
    if db is not None and source_db and source_db != db.name:
        raise ScenarioError(f"Scenario expects database {source_db!r}, got {db.name!r}")
    return document


def merge_changes(state: Any, changes: Any) -> Any:
    """Apply the recorder's recursive diff representation."""
    if not isinstance(changes, dict) or not isinstance(state, dict):
        return json.loads(json.dumps(changes))
    merged = json.loads(json.dumps(state))
    for key, value in changes.items():
        if value is None:
            merged.pop(key, None)
        elif isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_changes(merged[key], value)
        else:
            merged[key] = json.loads(json.dumps(value))
    return merged


def apply_state(page: Page, state: dict[str, Any], changed: dict[str, Any], timeout_ms: int) -> None:
    """Project supported SettingsCollector ontology fields into live controls."""
    del state
    nonce = str(time.time_ns())
    request = {"nonce": nonce, "settings": changed}
    set_model(page, "benchmark-scenario-restore", "value", json.dumps(request, sort_keys=True))
    page.wait_for_function(
        """nonce => {
            const model = Bokeh.documents[0].get_model_by_name('benchmark-scenario-restore-status');
            if (!model || !model.value) return false;
            const result = JSON.parse(model.value);
            if (result.nonce !== nonce) return false;
            if (result.status === 'failed') throw new Error(result.error || 'scenario restoration failed');
            return result.status === 'completed';
        }""",
        nonce,
        timeout=timeout_ms,
    )


def choose_filter_lookup_index(candidates: list[int], occurrence: int) -> int:
    """Choose a live duplicate filter target using a one-based occurrence."""
    if not isinstance(occurrence, int) or isinstance(occurrence, bool) or occurrence < 1:
        raise ScenarioError("filter_lookup occurrence must be a positive integer")
    ordered = sorted(set(candidates))
    if occurrence > len(ordered):
        raise ScenarioError(
            f"filter_lookup occurrence {occurrence} does not exist; found {len(ordered)} matching live rows"
        )
    return ordered[occurrence - 1]


def find_filter_lookup_index(page: Page, details: dict[str, Any]) -> int:
    """Resolve a semantic filter target against live, monotonically named rows."""
    candidates = page.evaluate(
        r"""({category, column}) => {
            const document = Bokeh.documents[0];
            const models = document._all_models;
            const values = models && typeof models.values === 'function'
                ? Array.from(models.values())
                : Array.from((models && models.values ? models.values() : []) || []);
            const indices = [];
            for (const model of values) {
                const match = /^benchmark-filter-(\d+)-category$/.exec(model.name || '');
                if (!match || model.value !== category) continue;
                const index = Number(match[1]);
                const metric = document.get_model_by_name(`benchmark-filter-${index}-metric`);
                if (metric && metric.value === column) indices.push(index);
            }
            return indices;
        }""",
        {"category": details.get("category"), "column": details.get("column")},
    )
    return choose_filter_lookup_index(candidates, details.get("occurrence", 1))


def replay_action(
    page: Page,
    step: dict[str, Any],
    state: dict[str, Any],
    timeout_ms: int,
) -> dict[str, Any]:
    action = step["action"]
    if action == "state_change":
        state = merge_changes(state, step["changes"])
        apply_state(page, state, step["changes"], timeout_ms)
    elif action == "filter_lookup":
        details = step["details"]
        index = find_filter_lookup_index(page, details)
        click_css(page, f".benchmark-filter-{index}-lookup")
    elif action == "apply_filters":
        click_css(page, ".benchmark-apply-filters")
    elif action == "apply_plot":
        click_css(page, ".benchmark-apply-plot")
    return state


def run_step(
    page: Page,
    server: Server | None,
    run_name: str,
    output: Path,
    step: dict[str, Any],
    action: Callable[[], None],
    timeout_ms: int,
) -> StepResult:
    mark = server.mark() if server else 0
    started = time.perf_counter()
    error = ""
    ok = True
    try:
        action()
        wait_sidebar(page, timeout_ms)
    except Exception as exc:  # keep later benchmark steps and artifacts
        ok = False
        error = f"{type(exc).__name__}: {exc}"
    elapsed = time.perf_counter() - started
    sequence = step["sequence"]
    action_name = step["action"]
    shot = output / f"{sequence:03d}-{re.sub(r'[^a-z0-9]+', '-', action_name.lower()).strip('-')}.png"
    try:
        page.screenshot(path=str(shot), full_page=True)
    except Exception as exc:
        ok = False
        error = f"{error}; screenshot: {exc}" if error else f"screenshot: {exc}"
    relevant = []
    if server:
        relevant = [
            line
            for line in server.since(mark)
            if line.startswith(
                ("[timing]", "[summary]", "[profile]", "[diagnostic]", "[rss]", "[settings]", "[benchmark]")
            )
        ]
        compatibility_warnings = [line for line in relevant if line.startswith("[settings] Warning")]
        if compatibility_warnings:
            ok = False
            warning = "control compatibility: " + " | ".join(compatibility_warnings)
            error = f"{error}; {warning}" if error else warning
    status = "completed" if ok else "failed"
    print(f"[{run_name}] {sequence}. {action_name}: {elapsed:.3f}s {status.upper()}", flush=True)
    return StepResult(sequence, action_name, status, elapsed, relevant, error, {"screenshot": str(shot)})


def workflow(
    page: Page,
    server: Server | None,
    run_name: str,
    output: Path,
    timeout_ms: int,
    scenario: dict[str, Any],
) -> list[StepResult]:
    results: list[StepResult] = []
    state = json.loads(json.dumps(scenario["initial_state"]))
    initial_mark = server.mark() if server else 0
    apply_state(page, state, state, timeout_ms)
    wait_sidebar(page, timeout_ms)
    if server:
        warnings = [line for line in server.since(initial_mark) if line.startswith("[settings] Warning")]
        if warnings:
            raise ScenarioError("Initial state is incompatible with server controls: " + " | ".join(warnings))
    for scenario_step in scenario["steps"]:

        def action(scenario_step=scenario_step) -> None:
            nonlocal state
            state = replay_action(page, scenario_step, state, timeout_ms)
        results.append(run_step(page, server, run_name, output, scenario_step, action, timeout_ms))
    return results


def augment_scenario(scenario: dict[str, Any], results: list[StepResult], metadata: dict[str, Any]) -> dict[str, Any]:
    document = json.loads(json.dumps(scenario))
    document["_result"] = metadata
    by_sequence = {result.sequence: result for result in results}
    for step in document["steps"]:
        result = by_sequence[step["sequence"]]
        step.update(asdict(result))
    return document


def write_results(
    output: Path,
    scenario: dict[str, Any],
    runs: dict[str, list[StepResult]],
    metadata: dict[str, Any],
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for run_name, results in runs.items():
        result_metadata = {**metadata, "run": run_name}
        (output / f"results.{run_name}.json").write_text(
            json.dumps(augment_scenario(scenario, results, result_metadata), indent=2) + "\n",
            encoding="utf-8",
        )
    with (output / "results.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("run", "sequence", "action", "status", "duration_seconds", "error", "artifacts", "server_lines"),
            delimiter="\t",
        )
        writer.writeheader()
        for run_name, results in runs.items():
            for result in results:
                row = {"run": run_name, **asdict(result)}
                row["server_lines"] = " | ".join(result.server_lines)
                row["artifacts"] = json.dumps(result.artifacts, sort_keys=True)
                writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", required=True, type=Path, help="Recorded thebigbam scenario JSON")
    parser.add_argument("--db", type=Path, help="Database to serve (required unless --url is used)")
    parser.add_argument("--url", help="Use an existing server instead of starting one")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("benchmarks/plotting/results"))
    parser.add_argument("--cold-and-warm", action="store_true")
    parser.add_argument("--headed", action="store_true")
    parser.add_argument("--timeout", type=float, default=180.0, help="Per-step and startup timeout in seconds")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    scenario_path = (root / args.scenario).resolve() if not args.scenario.is_absolute() else args.scenario.resolve()
    if not args.url and args.db is None:
        raise SystemExit("--db is required unless --url is used")
    db = None if args.db is None else ((root / args.db).resolve() if not args.db.is_absolute() else args.db.resolve())
    if not args.url and db is not None and not db.exists():
        raise SystemExit(f"Database does not exist: {db}")
    output = args.output.resolve()
    try:
        scenario = load_scenario(scenario_path, None if args.url else db)
    except ScenarioError as exc:
        raise SystemExit(f"Invalid scenario: {exc}") from exc
    initial_settings = output / ".initial-settings.json"
    output.mkdir(parents=True, exist_ok=True)
    initial_settings.write_text(json.dumps(scenario["initial_state"], indent=2) + "\n", encoding="utf-8")
    port = args.port or free_port()
    server = None if args.url else Server(db, port, output, initial_settings)
    url = args.url or f"http://127.0.0.1:{port}"
    if server:
        server.start(args.timeout)

    runs: dict[str, list[StepResult]] = {}
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=not args.headed)
            run_names = ("cold", "warm") if args.cold_and_warm else ("cold",)
            for run_name in run_names:
                context = browser.new_context(viewport={"width": 1920, "height": 1080})
                page = context.new_page()
                page.goto(url, wait_until="domcontentloaded", timeout=int(args.timeout * 1000))
                page.wait_for_function("window.Bokeh && Bokeh.documents.length", timeout=int(args.timeout * 1000))
                run_output = output / run_name
                run_output.mkdir(parents=True, exist_ok=True)
                runs[run_name] = workflow(
                    page, server, run_name, run_output, int(args.timeout * 1000), scenario
                )
                context.close()
            browser.close()
    finally:
        if server:
            server.stop()

    metadata = {
        "database": str(db) if db is not None else None,
        "scenario": str(scenario_path),
        "url": url,
        "cold_and_warm": args.cold_and_warm,
        "python": sys.executable,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    write_results(output, scenario, runs, metadata)
    all_results = [result for results in runs.values() for result in results]
    failed = [result for result in all_results if result.status == "failed"]
    print(f"Wrote {len(all_results)} step results to {output}; failures={len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
