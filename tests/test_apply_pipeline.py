from types import SimpleNamespace

import pytest

from thebigbam.plotting.application.apply_pipeline import (
    ApplyDispatcher,
    ApplyMode,
    ApplyRequestBuilder,
    ContigAllApplyRequest,
    ContigOneApplyRequest,
    MagAllApplyRequest,
    MagOneApplyRequest,
    PlotPresenter,
    RequestBuildResult,
    ValidationFailure,
)
from thebigbam.plotting.shared.timing import PipelineTimings


def test_dispatcher_routes_all_four_modes_explicitly():
    calls = []
    handlers = {mode: lambda inputs, started, mode=mode: calls.append(mode) or mode.value for mode in ApplyMode}
    dispatcher = ApplyDispatcher(handlers)
    inputs = SimpleNamespace()

    requests = (
        ContigOneApplyRequest(inputs),
        ContigAllApplyRequest(inputs),
        MagOneApplyRequest(inputs),
        MagAllApplyRequest(inputs),
    )
    for request in requests:
        assert dispatcher.render(request, 1.0) == request.mode.value

    assert calls == list(ApplyMode)


def test_dispatcher_rejects_an_incomplete_mode_table():
    with pytest.raises(ValueError, match="Missing Apply handlers"):
        ApplyDispatcher({ApplyMode.CONTIG_ONE: lambda inputs, started: None})


def test_request_build_result_exposes_validation_state():
    valid = RequestBuildResult(request=ContigOneApplyRequest(SimpleNamespace()))
    invalid = RequestBuildResult(failure=ValidationFailure("bad request"))

    assert valid.is_valid
    assert not invalid.is_valid


def test_request_builder_reports_invalid_positions_without_rendering(monkeypatch):
    inputs = SimpleNamespace(contig="c1", is_all=False)
    bindings = SimpleNamespace(
        widgets={
            "has_mags": False,
            "view_radio": SimpleNamespace(active=1),
            "contig_lengths": {"c1": 100},
        },
        from_position_input=SimpleNamespace(value="bad"),
        to_position_input=SimpleNamespace(value="100"),
    )
    monkeypatch.setattr(
        "thebigbam.plotting.application.apply_pipeline.collect_apply_inputs",
        lambda bindings: inputs,
    )

    result = ApplyRequestBuilder(bindings).from_widgets()

    assert result.failure == ValidationFailure("Invalid position range - positions must be integers.")


def test_presenter_owns_validation_display_and_button_state():
    buttons = [SimpleNamespace(visible=True) for _ in range(4)]
    bindings = SimpleNamespace(
        peruse_button=buttons[0],
        download_metrics_button=buttons[1],
        download_mag_metrics_button=buttons[2],
        download_data_button=buttons[3],
        main_placeholder=SimpleNamespace(objects=[]),
    )

    PlotPresenter(bindings).show_validation(ValidationFailure("Please select a contig."))

    assert all(not button.visible for button in buttons)
    assert "Please select a contig" in bindings.main_placeholder.objects[0].object


def test_presenter_owns_range_callback_registration():
    callbacks = {}
    x_range = SimpleNamespace(on_change=lambda attribute, callback: callbacks.setdefault(attribute, callback))
    state = {"range_callbacks": ()}
    bindings = SimpleNamespace(
        current_plot_state=state,
        from_position_input=SimpleNamespace(value=""),
        to_position_input=SimpleNamespace(value=""),
    )

    PlotPresenter(bindings).bind_range_inputs(x_range)
    callbacks["start"]("start", 1, 12.8)
    callbacks["end"]("end", 20, 42.2)

    assert bindings.from_position_input.value == "12"
    assert bindings.to_position_input.value == "42"
    assert len(state["range_callbacks"]) == 2


def test_presenter_preserves_and_installs_matching_range():
    callbacks = {}
    old_range = SimpleNamespace(start=10, end=20)
    new_range = SimpleNamespace(
        start=1,
        end=100,
        on_change=lambda attribute, callback: callbacks.setdefault(attribute, callback),
    )
    state = {
        "shared_xrange": old_range,
        "contig": "c1",
        "sample": "s1",
        "is_all": False,
        "data_xstart": 1,
        "data_xend": 100,
        "range_callbacks": (),
    }
    bindings = SimpleNamespace(
        current_plot_state=state,
        plot_lifecycle=SimpleNamespace(shared_xrange=lambda grid: new_range),
        from_position_input=SimpleNamespace(value=""),
        to_position_input=SimpleNamespace(value=""),
    )
    presenter = PlotPresenter(bindings)
    snapshot = presenter.capture_range("c1", "s1", False, 1, 100)

    presenter.install_range(object(), snapshot, subject="c1", sample="s1", is_all=False, start=1, end=100)

    assert (new_range.start, new_range.end) == (10, 20)
    assert state["shared_xrange"] is new_range


def test_pipeline_timings_accumulate_named_phases(monkeypatch):
    ticks = iter((1.0, 1.2, 2.0, 2.5, 3.0, 3.4))
    monkeypatch.setattr("thebigbam.plotting.shared.timing.time.perf_counter", lambda: next(ticks))
    timings = PipelineTimings()

    with timings.phase("repository"):
        pass
    with timings.phase("service"):
        pass
    with timings.phase("repository"):
        pass

    assert timings.seconds["repository"] == pytest.approx(0.6)
    assert timings.seconds["service"] == pytest.approx(0.5)
    assert "repository=0.600s" in timings.format()
