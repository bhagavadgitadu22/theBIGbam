from types import SimpleNamespace

import pytest
from bokeh.models import Select

from thebigbam.plotting.application.availability_facade import AvailabilityFacade
from thebigbam.plotting.application.scope_transition import ScopeTransitionBindings, ScopeTransitionController
from thebigbam.plotting.application.state import SampleScope, SubjectScope


class Diagnostics:
    def record(self, *args, **kwargs):
        return None


def _visible():
    return SimpleNamespace(visible=False)


def test_scope_change_revalidates_stale_sample_options_before_mag_options():
    sample = SimpleNamespace(options=["valid", "stale"], value="stale")
    mag = SimpleNamespace(options=["m1"], value="m1")
    widgets = {
        "sample_select": sample,
        "contig_select": SimpleNamespace(options=["c1"], value="c1"),
        "mag_select": mag,
        "has_mags": True,
        "view_radio": SimpleNamespace(active=0),
        "variables_widgets_all": [],
    }
    refresh_order = []

    def update(widget, options):
        widget.options = list(options)
        if widget.value not in widget.options:
            widget.value = ""

    def refresh_availability():
        refresh_order.append("samples")
        update(sample, ["valid"])
        refresh_order.append("mags")
        update(mag, ["m1"])

    controller = ScopeTransitionController(
        ScopeTransitionBindings(
            widgets=widgets,
            sample_scope=SimpleNamespace(active=1),
            sample_section=_visible(),
            variables_one=_visible(),
            variables_all=_visible(),
            sample_parameters_header=_visible(),
            sample_parameters_content=_visible(),
            mag_sort_sample_row=_visible(),
            mag_sort_category=SimpleNamespace(value="Contig"),
            interaction_lock={"locked": False},
            diagnostics=Diagnostics(),
            enable_timing=False,
            timing=None,
            send_timing_ping=lambda *args: None,
            set_operation=lambda operation: None,
            update_completions=update,
            refresh_availability=refresh_availability,
        )
    )

    controller.view_changed("active", 0, 1)

    assert refresh_order == ["samples", "mags"]
    assert sample.options == ["valid"]
    assert sample.value == ""
    assert mag.value == "m1"


@pytest.mark.parametrize("has_mags", [False, True])
@pytest.mark.parametrize("valid_contig", [False, True])
def test_scope_switch_revalidates_search_selection_before_projecting_cached_options(has_mags, valid_contig):
    first_page = [f"c{i}" for i in range(100)]
    contig = Select(options=first_page, value="c0")
    sample = Select(options=["s1"], value="s1")
    mag = Select(options=["m1"], value="m1")
    scope = SimpleNamespace(active=0)
    widgets = {
        "sample_select": sample,
        "contig_select": contig,
        "mag_select": mag,
        "has_mags": has_mags,
        "view_radio": SimpleNamespace(active=1),
        "variables_widgets_all": [],
    }
    refresh_values = []

    def refresh_availability():
        refresh_values.append(contig.value)
        # Model bounded availability that prioritizes a valid current value.
        options = first_page
        if valid_contig and contig.value == "c149":
            options = ["c149", *first_page[:99]]
        AvailabilityFacade.replace_widget(sample, ["s1"])
        if has_mags:
            AvailabilityFacade.replace_widget(mag, ["m1"])
        AvailabilityFacade.replace_widget(contig, options)

    controller = ScopeTransitionController(
        ScopeTransitionBindings(
            widgets=widgets,
            sample_scope=scope,
            sample_section=_visible(),
            variables_one=_visible(),
            variables_all=_visible(),
            sample_parameters_header=_visible(),
            sample_parameters_content=_visible(),
            mag_sort_sample_row=_visible(),
            mag_sort_category=SimpleNamespace(value="Contig"),
            interaction_lock={"locked": False},
            diagnostics=Diagnostics(),
            enable_timing=False,
            timing=None,
            send_timing_ping=lambda *args: None,
            set_operation=lambda operation: None,
            update_completions=AvailabilityFacade.replace_widget,
            refresh_availability=refresh_availability,
        )
    )
    # Search selects a contig absent from both startup scope snapshots.
    contig.options = ["c149"]
    contig.value = "c149"
    changes = []
    contig.on_change("value", lambda attr, old, new: changes.append(new))
    expected = "c149" if valid_contig else ""

    for old, new in [(0, 1), (1, 0), (0, 1)]:
        scope.active = new
        controller.view_changed("active", old, new)
        assert contig.value == expected
        assert controller.plot_controller.state.contig == expected
        assert controller.plot_controller.state.sample == "s1"
        assert controller.plot_controller.state.subject_scope == SubjectScope.CONTIG
        assert controller.plot_controller.state.sample_scope == SampleScope(new)
        assert controller.bindings.interaction_lock["locked"] is False

    assert refresh_values == ["c149", expected, expected]
    # A valid contig must never clear transiently: that can reset positions.
    assert changes == ([] if valid_contig else [""])
