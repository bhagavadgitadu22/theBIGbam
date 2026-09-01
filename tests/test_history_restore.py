from thebigbam.plotting.application.history_restore import (
    HistoryRestoreBindings,
    HistoryRestoreCoordinator,
)
from thebigbam.plotting.settings.history import HistoryEntry


def _entry(action="apply_plot"):
    return HistoryEntry("id", 1, action, "time", {"selection": {"contig": "restored"}})


def _coordinator(queue, events, *, apply_plot=lambda: True, install=None, snapshot=None):
    return HistoryRestoreCoordinator(
        HistoryRestoreBindings(
            begin=lambda: events.append("begin") or True,
            end=lambda: events.append("end"),
            set_loading=lambda loading: events.append(("loading", loading)),
            schedule=queue.append,
            snapshot=snapshot or (lambda: {"selection": {"contig": "current"}}),
            install=install or (lambda settings: events.append(("install", settings))),
            apply_plot=apply_plot,
        )
    )


def test_plot_restore_installs_snapshot_regenerates_plot_and_then_completes():
    queue = []
    events = []
    coordinator = _coordinator(queue, events, apply_plot=lambda: events.append("apply_plot") or True)

    assert coordinator.restore(_entry(), lambda: events.append("complete")) is True
    queue.pop(0)()

    assert events == [
        "begin",
        ("loading", True),
        ("install", {"selection": {"contig": "restored"}}),
        "apply_plot",
        ("loading", False),
        "end",
        "complete",
    ]


def test_filter_restore_installs_snapshot_without_regenerating_plot():
    queue = []
    events = []
    coordinator = _coordinator(
        queue, events, apply_plot=lambda: events.append("unexpected apply") or True
    )

    coordinator.restore(_entry("apply_filters"))
    queue.pop(0)()

    assert ("install", {"selection": {"contig": "restored"}}) in events
    assert "unexpected apply" not in events


def test_failed_plot_restore_rolls_back_state_and_plot_before_releasing_lock():
    queue = []
    events = []
    failures = []
    apply_results = iter([False, True])
    coordinator = _coordinator(
        queue,
        events,
        apply_plot=lambda: events.append("apply_plot") or next(apply_results),
    )

    coordinator.restore(_entry(), lambda: events.append("complete"), failures.append)
    queue.pop(0)()

    assert events == [
        "begin",
        ("loading", True),
        ("install", {"selection": {"contig": "restored"}}),
        "apply_plot",
        ("install", {"selection": {"contig": "current"}}),
        "apply_plot",
        ("loading", False),
        "end",
    ]
    assert [str(error) for error in failures] == [
        "previously successful plot snapshot could not be regenerated"
    ]


def test_snapshot_failure_still_releases_lock_and_reports_failure():
    queue = []
    events = []
    failures = []
    coordinator = _coordinator(
        queue,
        events,
        snapshot=lambda: (_ for _ in ()).throw(RuntimeError("snapshot failed")),
    )

    coordinator.restore(_entry(), on_error=failures.append)
    queue.pop(0)()

    assert events == ["begin", ("loading", True), ("loading", False), "end"]
    assert [str(error) for error in failures] == ["snapshot failed"]


def test_restore_rejection_does_not_schedule_partial_work():
    queue = []
    coordinator = HistoryRestoreCoordinator(
        HistoryRestoreBindings(
            begin=lambda: False,
            end=lambda: None,
            set_loading=lambda _loading: None,
            schedule=queue.append,
            snapshot=dict,
            install=lambda _settings: None,
            apply_plot=lambda: True,
        )
    )

    assert coordinator.restore(_entry()) is False
    assert queue == []
