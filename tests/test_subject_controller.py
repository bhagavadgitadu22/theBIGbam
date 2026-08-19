from types import SimpleNamespace

from thebigbam.plotting.application.subject import SubjectBindings, SubjectController


def test_subject_controller_projects_selected_contig_into_mag_coordinates():
    widgets = {
        "has_mags": True,
        "view_radio": SimpleNamespace(active=0),
        "mag_select": SimpleNamespace(value="m1"),
        "contig_select": SimpleNamespace(value="c2"),
        "mag_to_contigs": {"m1": ["c1", "c2"]},
        "mag_to_contig_offsets": {},
        "contig_lengths": {"c1": 10, "c2": 20},
    }
    start = SimpleNamespace(value="")
    end = SimpleNamespace(value="")
    controller = SubjectController(
        SubjectBindings(
            widgets=widgets,
            interaction_lock={"locked": False},
            compute_contigs=lambda search: [],
            compute_samples=lambda search: [],
            compute_mags=lambda search: [],
            push_completions=lambda widget, values: None,
            refresh_contigs=lambda: None,
            refresh_samples=lambda: None,
            refresh_mags=lambda: None,
            update_titles=lambda: None,
            from_position=start,
            to_position=end,
            sample_order_category=SimpleNamespace(options=[], value=""),
            sample_contig_categories=["Sample"],
            sample_mag_categories=["MAG coverage"],
            sample_current_categories=["Sample"],
        )
    )

    controller.sync_selected_contig_position()

    assert start.value == "11"
    assert end.value == "30"
    assert widgets["mag_to_contig_offsets"]["m1"] == {"c1": 0, "c2": 10}


def test_invalid_filtered_sample_is_rejected_without_clearing_mag():
    sample = SimpleNamespace(value="stale")
    widgets = {
        "has_mags": True,
        "sample_select": sample,
        "mag_select": SimpleNamespace(value="m1"),
    }
    refreshed_mags = []

    def push(widget, values):
        widget.options = values
        if widget.value not in values:
            widget.value = ""

    controller = SubjectController(
        SubjectBindings(
            widgets=widgets,
            interaction_lock={"locked": False},
            compute_contigs=lambda search: [],
            compute_samples=lambda search: ["valid"],
            compute_mags=lambda search: ["m1"],
            push_completions=push,
            refresh_contigs=lambda: None,
            refresh_samples=lambda: None,
            refresh_mags=lambda: refreshed_mags.append(True),
            update_titles=lambda: None,
            from_position=SimpleNamespace(value=""),
            to_position=SimpleNamespace(value=""),
            sample_order_category=SimpleNamespace(options=[], value=""),
            sample_contig_categories=[],
            sample_mag_categories=[],
            sample_current_categories=[],
        )
    )

    controller.sample_changed(SimpleNamespace(new="stale"))

    assert sample.value == ""
    assert sample.options == ["valid"]
    assert widgets["mag_select"].value == "m1"
    assert refreshed_mags == []
