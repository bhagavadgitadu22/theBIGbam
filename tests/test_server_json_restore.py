import argparse

from thebigbam.plotting.application.server import add_serve_args, restore_payload
from thebigbam.plotting.settings.history import SessionHistory


def test_settings_json_remains_a_settings_only_restore():
    settings = {"selection": {"contig": "c1"}}

    assert restore_payload(settings) == (settings, ())


def test_session_json_restores_history_and_uses_latest_entry_for_controls():
    history = SessionHistory("example.db")
    first = history.append("apply_filters", {"filtering": []})
    latest = history.append("apply_plot", {"selection": {"contig": "c2"}})

    settings, entries = restore_payload(history.document())

    assert settings == latest.settings
    assert entries == (first, latest)


def test_json_help_documents_both_saved_formats():
    parser = argparse.ArgumentParser()
    add_serve_args(parser)

    help_text = " ".join(parser.format_help().split())

    assert "SAVE SETTINGS or SAVE SESSION" in help_text
    assert "complete application history" in help_text
