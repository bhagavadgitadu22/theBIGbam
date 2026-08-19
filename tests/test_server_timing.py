from types import SimpleNamespace

from thebigbam.plotting.application import server


def _run(monkeypatch, *, enable_timing):
    watchdog_calls = []
    serve_kwargs = {}
    monkeypatch.setattr(server.application, "_start_rss_watchdog", lambda: watchdog_calls.append(True))
    monkeypatch.setattr(server, "_print_database_metadata", lambda _path: None)
    monkeypatch.setattr(server.application, "preload_db_data", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(server, "warm_composer_imports", lambda: None)
    monkeypatch.setattr(server.pn, "serve", lambda *_args, **kwargs: serve_kwargs.update(kwargs))
    result = server.run_serve(
        SimpleNamespace(
            db="example.db",
            port=5006,
            time=enable_timing,
            settings_json=None,
        )
    )
    return result, watchdog_calls, serve_kwargs


def test_normal_server_does_not_start_rss_watchdog(monkeypatch):
    result, watchdog_calls, serve_kwargs = _run(monkeypatch, enable_timing=False)

    assert result == 0
    assert watchdog_calls == []
    assert serve_kwargs["allow_websocket_origin"] == ["localhost:5006", "127.0.0.1:5006"]


def test_timed_server_starts_rss_watchdog(monkeypatch):
    result, watchdog_calls, _serve_kwargs = _run(monkeypatch, enable_timing=True)

    assert result == 0
    assert watchdog_calls == [True]
