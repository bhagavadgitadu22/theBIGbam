from types import SimpleNamespace

from thebigbam.plotting.application.scope_transition import ScopeTransitionBindings, ScopeTransitionController


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
