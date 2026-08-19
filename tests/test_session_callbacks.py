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

    callbacks.apply_clicked()

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
