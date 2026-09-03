from types import SimpleNamespace

import pytest

from thebigbam.plotting.application.session_callbacks import SessionCallbacks


def test_apply_is_scheduled_and_uses_attached_controller():
    scheduled = []
    calls = []
    callbacks = SessionCallbacks(scheduled.append)
    placeholder = SimpleNamespace(loading=False)
    callbacks.attach_placeholder(placeholder)
    callbacks.attach_apply(SimpleNamespace(apply=lambda: calls.append("apply")))

    assert callbacks.apply_clicked() is True

    assert placeholder.loading is True
    assert scheduled == [callbacks.do_apply]
    scheduled[0]()
    assert calls == ["apply"]


def test_summary_subject_and_scope_delegation():
    calls = []
    callbacks = SessionCallbacks(lambda callback: None)
    callbacks.attach_summary(SimpleNamespace(show=lambda: calls.append("summary")))
    callbacks.attach_subject(SimpleNamespace(sync_selected_contig_position=lambda: calls.append("sync")))
    callbacks.attach_scope(
        SimpleNamespace(
            variable_callback=lambda group: ("callback", group),
            view_changed=lambda attr, old, new: calls.append((attr, old, new)),
        )
    )

    callbacks.show_summary()
    callbacks.sync_selected_contig_position()
    assert callbacks.variable_callback("variables") == ("callback", "variables")
    callbacks.view_changed("active", 0, 1)

    assert calls == ["summary", "sync", ("active", 0, 1)]


def test_unattached_controller_fails_explicitly():
    callbacks = SessionCallbacks(lambda callback: None)
    with pytest.raises(RuntimeError, match="Summary controller"):
        callbacks.show_summary()


def test_plot_apply_holds_interaction_scope_until_controller_finishes():
    scheduled = []
    events = []
    interactions = SimpleNamespace(
        begin=lambda scope: events.append(("begin", scope)) or True,
        end=lambda: events.append(("end", None)),
    )
    callbacks = SessionCallbacks(scheduled.append, interactions=interactions)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: events.append(("apply", None))))

    callbacks.apply_clicked()
    scheduled[0]()

    assert events == [("begin", "plot"), ("apply", None), ("end", None)]


def test_rejected_plot_apply_reports_that_it_did_not_start():
    placeholder = SimpleNamespace(loading=False)
    callbacks = SessionCallbacks(
        lambda _callback: None,
        interactions=SimpleNamespace(begin=lambda _scope: False, end=lambda: None),
    )
    callbacks.attach_placeholder(placeholder)

    assert callbacks.apply_clicked() is False
    assert placeholder.loading is False


def test_apply_failure_always_clears_loading():
    scheduled = []
    placeholder = SimpleNamespace(loading=False)
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(placeholder)
    callbacks.attach_apply(
        SimpleNamespace(apply=lambda: (_ for _ in ()).throw(RuntimeError("failed before presentation")))
    )

    callbacks.apply_clicked()
    with pytest.raises(RuntimeError, match="failed before presentation"):
        scheduled.pop(0)()

    assert placeholder.loading is False


def test_plot_apply_records_semantic_action():
    scheduled = []
    events = []
    recorder = SimpleNamespace(record_action=lambda action, settings: events.append(("record", action, settings)))
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: events.append(("apply",))))
    callbacks.attach_scenario(recorder, lambda: {"selection": {"sample": "s1"}})

    callbacks.apply_clicked()
    scheduled[0]()

    assert events == [
        ("record", "apply_plot", {"selection": {"sample": "s1"}}),
        ("apply",),
    ]


def test_failed_plot_apply_is_recorded_and_reraised():
    scheduled = []
    recorded = []
    recorder = SimpleNamespace(record_action=lambda action, settings: recorded.append((action, settings)))
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: (_ for _ in ()).throw(RuntimeError("plot failed"))))
    callbacks.attach_scenario(recorder, dict)

    callbacks.apply_clicked()
    with pytest.raises(RuntimeError, match="plot failed"):
        scheduled[0]()

    assert recorded == [("apply_plot", {})]


def test_successful_plot_appends_click_time_applied_settings():
    scheduled = []
    settings = {"selection": {"contig": "c1"}}
    appended = []
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: True))
    callbacks.attach_history(
        lambda: {"selection": dict(settings["selection"])},
        lambda snapshot, _apply_step: appended.append(snapshot),
    )

    callbacks.apply_clicked()
    settings["selection"]["contig"] = "c2"
    scheduled[0]()

    assert appended == [{"selection": {"contig": "c1"}}]


def test_successful_plot_history_retains_originating_scenario_step():
    scheduled = []
    appended = []
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: True))
    callbacks.attach_scenario(
        SimpleNamespace(record_action=lambda _action, _settings: 17),
        lambda: {"selection": {"contig": "c1"}},
    )
    callbacks.attach_history(
        lambda: {"selection": {"contig": "c1"}},
        lambda settings, apply_step: appended.append((settings, apply_step)),
    )

    callbacks.apply_clicked()
    scheduled[0]()

    assert appended == [({"selection": {"contig": "c1"}}, 17)]


def test_failed_plot_does_not_append_history():
    scheduled = []
    appended = []
    callbacks = SessionCallbacks(scheduled.append)
    callbacks.attach_placeholder(SimpleNamespace(loading=False))
    callbacks.attach_apply(SimpleNamespace(apply=lambda: False))
    callbacks.attach_history(lambda: {"selection": {}}, appended.append)

    callbacks.apply_clicked()
    scheduled.pop(0)()

    assert appended == []


def test_atomic_restored_apply_uses_outer_loading_and_interaction_transaction():
    events = []
    placeholder = SimpleNamespace(loading=False)
    callbacks = SessionCallbacks(
        lambda _callback: None,
        interactions=SimpleNamespace(
            begin=lambda _scope: events.append("unexpected begin"),
            end=lambda: events.append("unexpected end"),
        ),
    )
    callbacks.attach_placeholder(placeholder)
    callbacks.attach_apply(SimpleNamespace(apply=lambda: events.append("apply") or True))

    assert callbacks.apply_restored_now() is True
    assert placeholder.loading is False
    assert events == ["apply"]

    callbacks.set_plot_loading(True)
    assert placeholder.loading is True
    callbacks.set_plot_loading(False)
    assert placeholder.loading is False
