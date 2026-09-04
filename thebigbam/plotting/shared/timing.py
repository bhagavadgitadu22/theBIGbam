"""Low-overhead server and Bokeh payload diagnostics."""

from __future__ import annotations

import sys
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Callable

from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import TextInput


class PipelineTimings:
    """Accumulate exclusive wall-clock phases for one plot composition."""

    def __init__(self) -> None:
        self.seconds: dict[str, float] = {}

    @contextmanager
    def phase(self, name: str):
        started = time.perf_counter()
        try:
            yield
        finally:
            self.seconds[name] = self.seconds.get(name, 0.0) + time.perf_counter() - started

    def format(self) -> str:
        return " ".join(f"{name}={seconds:.3f}s" for name, seconds in self.seconds.items())


@dataclass
class ApplyProfiler:
    """Structured low-overhead timings for one Apply render request."""

    mode: str
    subject: str
    sample: str | None = None
    enabled: bool = False
    seconds: dict[str, float] = field(default_factory=dict)
    cache_hits: dict[str, int] = field(default_factory=dict)
    cache_misses: dict[str, int] = field(default_factory=dict)
    metrics: dict[str, int | float | str] = field(default_factory=dict)

    @contextmanager
    def phase(self, name: str):
        started = time.perf_counter()
        try:
            yield
        finally:
            self.seconds[name] = self.seconds.get(name, 0.0) + time.perf_counter() - started

    def cache(self, name: str, hit: bool) -> None:
        target = self.cache_hits if hit else self.cache_misses
        target[name] = target.get(name, 0) + 1

    def metric(self, name: str, value: int | float | str) -> None:
        self.metrics[name] = value

    def report(self, cache_stats=None) -> None:
        if not self.enabled:
            return
        phases = " ".join(f"{name}={seconds:.3f}s" for name, seconds in self.seconds.items())
        keys = set(self.cache_hits) | set(self.cache_misses)
        caches = " ".join(
            f"{name}=hit:{self.cache_hits.get(name, 0)}/miss:{self.cache_misses.get(name, 0)}" for name in sorted(keys)
        )
        retained = f" entries={cache_stats.entries}" if cache_stats is not None else ""
        metrics = " ".join(f"{name}={value}" for name, value in sorted(self.metrics.items()))
        print(
            f"[profile] mode={self.mode} subject={self.subject} sample={self.sample or '-'} "
            f"{phases} {metrics} cache[{caches}]{retained}",
            flush=True,
        )


@contextmanager
def profile_phase(profiler: ApplyProfiler | None, name: str):
    if profiler is None:
        yield
    else:
        with profiler.phase(name):
            yield


def rss_mb() -> float:
    try:
        with open("/proc/self/status") as status:
            for line in status:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        pass
    try:
        import resource

        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    except (ImportError, OSError):
        return 0.0


def estimate_grid_data_size(grid) -> tuple[int, int, int]:
    total_bytes = 0
    sources = 0
    references = grid.references()
    for reference in references:
        if isinstance(reference, ColumnDataSource):
            sources += 1
            for column in reference.data.values():
                total_bytes += column.nbytes if hasattr(column, "nbytes") else sys.getsizeof(column)
    return total_bytes, sources, len(references)


class TimingPhase:
    def __init__(self) -> None:
        self._phase_start = time.perf_counter()
        self._phase_timed = 0.0
        self._phase_name = "Startup"

    def tag(self, step_seconds: float) -> str:
        self._phase_timed += step_seconds
        real = time.perf_counter() - self._phase_start
        return f" (timed={self._phase_timed:.3f}s real={real:.3f}s)"

    def summary(self, name: str | None = None) -> None:
        label = name or self._phase_name
        real = time.perf_counter() - self._phase_start
        gap = real - self._phase_timed
        percentage = gap / real * 100 if real > 0 else 0
        print(
            f"[summary] {label} done: timed={self._phase_timed:.3f}s "
            f"real={real:.3f}s gap={gap:.3f}s ({percentage:.1f}%)",
            flush=True,
        )
        self._phase_start = time.perf_counter()
        self._phase_timed = 0.0

    def start_phase(self, name: str) -> None:
        self._phase_name = name
        self._phase_start = time.perf_counter()
        self._phase_timed = 0.0


def start_rss_watchdog(operation: Callable[[], str], interval: int = 5) -> None:
    def watch() -> None:
        while True:
            time.sleep(interval)
            print(f"[rss] watchdog: {rss_mb():.0f} MB  op={operation()}", flush=True)

    threading.Thread(target=watch, daemon=True, name="rss-watchdog").start()


class BrowserTimingRelay:
    """Measure browser paint separately from its acknowledgement round-trip."""

    def __init__(self, enabled: bool) -> None:
        self.enabled = enabled
        self.state: dict[str, object] = {}
        self.ping = TextInput(name="benchmark-render-ping", value="", visible=False)
        self.ack = TextInput(name="benchmark-render-ack", value="", visible=False)
        self.ping.js_on_change(
            "value",
            CustomJS(
                args=dict(ack=self.ack),
                code="""
                try {
                    if (cb_obj.value) {
                        const patchStarted = performance.now();
                        requestAnimationFrame(() => {
                            requestAnimationFrame(() => {
                                var mem = '';
                                if (window.performance && performance.memory) {
                                    mem = '|heap=' + Math.round(performance.memory.usedJSHeapSize / 1024 / 1024);
                                }
                                const painted = performance.now();
                                const paint = '|paint=' + (painted - patchStarted).toFixed(1);
                                const started = window.__thebigbam_view_change_started;
                                const interaction = started == null ? '' : '|interaction=' + (painted - started).toFixed(1);
                                window.__thebigbam_view_change_started = null;
                                ack.value = cb_obj.value + mem + paint + interaction;
                            });
                        });
                    }
                } catch(e) {}
                """,
            ),
        )
        if enabled:
            self.ack.on_change("value", self._on_ack)

    def send(self, label: str, flow_start: float | None = None) -> None:
        token = f"{label}_{time.perf_counter()}"
        self.state.clear()
        self.state.update(token=token, label=label, sent=time.perf_counter())
        if flow_start is not None:
            self.state["flow_start"] = flow_start
        self.ping.value = token

    @staticmethod
    def _field(value: str, name: str, default: str = "unknown") -> str:
        marker = f"|{name}="
        return value.split(marker)[1].split("|", 1)[0] if marker in value else default

    def _on_ack(self, _attr, _old, new: str) -> None:
        token = new.split("|", 1)[0] if new else ""
        if token != self.state.get("token") or "sent" not in self.state:
            return
        elapsed = time.perf_counter() - float(self.state["sent"])
        label = self.state.get("label", "unknown")
        heap = self._field(new, "heap", "")
        heap_text = f" [JS heap: {heap} MB]" if heap else ""
        print(
            f"[timing] Browser paint after patch '{label}': {self._field(new, 'paint')} ms; "
            f"interaction to paint: {self._field(new, 'interaction')} ms; "
            f"acknowledgement round-trip: {elapsed:.3f}s{heap_text}",
            flush=True,
        )
        if "flow_start" in self.state:
            total = time.perf_counter() - float(self.state["flow_start"])
            print(f"[timing] Total flow (query -> send -> render): {total:.3f}s", flush=True)
