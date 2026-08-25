from types import SimpleNamespace

import pytest

from thebigbam.plotting.application.availability import AvailabilityBindings, AvailabilityController


class AvailabilityService:
    def contigs_for_sample(self, sample_id, search, preserve, mag_name=None):
        return ("c2",)

    def samples_for_contig(self, contig_id, search):
        return ("s2",)

    def count_contigs_for_sample(self, sample_id, mag_name=None):
        return 1

    def count_samples_for_contig(self, contig_id):
        return 1


class FilteredMagAvailabilityService(AvailabilityService):
    def filtered_contigs(self, sql, parameters, **scope):
        return ("c1",)

    def count_filtered_contigs(self, sql, parameters, **scope):
        return 1

    def filtered_samples(self, sql, parameters, **scope):
        return ("GL34_UP",) if scope.get("mag_name") == "m1" else ("GL34_UP", "other")

    def count_filtered_samples(self, sql, parameters, **scope):
        return 1 if scope.get("mag_name") == "m1" else 2

    def filtered_mags(self, sql, parameters, **scope):
        return ("m1",)

    def count_filtered_mags(self, sql, parameters, **scope):
        return 1


def _select(value="", options=()):
    return SimpleNamespace(value=value, options=list(options))


def test_availability_controller_uses_service_for_one_sample_subject_constraints():
    widgets = {
        "has_mags": False,
        "view_radio": SimpleNamespace(active=1),
        "contig_select": _select(""),
        "sample_select": _select("s1"),
        "mag_select": _select(""),
        "sample_name_to_id": {"s1": 1},
        "contig_name_to_id": {"c1": 1},
    }
    title = SimpleNamespace(text="")
    controller = AvailabilityController(
        AvailabilityBindings(
            availability_service=AvailabilityService(),
            filtering_pairs=lambda: None,
            original_contigs=["c1", "c2"],
            original_samples=["s1", "s2"],
            sample_scope=SimpleNamespace(active=0),
            widgets=widgets,
            sort_sample_select=_select(),
            update_completions=lambda widget, values: setattr(widget, "options", list(values)),
            total_coverage_count=2,
            filtering_title=title,
            contig_title=SimpleNamespace(text=""),
            sample_title=SimpleNamespace(text=""),
            mag_title=SimpleNamespace(text=""),
        )
    )

    assert controller.compute_contigs() == ["c2"]
    assert controller.compute_samples() == ["s1", "s2"]
    assert controller.sort_samples_for_mag(["s1"]) == ["s1"]


def test_selected_sample_does_not_collapse_available_sample_count_to_one():
    widgets = {
        "has_mags": False,
        "view_radio": SimpleNamespace(active=1),
        "contig_select": _select(""),
        "sample_select": _select("s1", ("s1", "s2")),
        "mag_select": _select(""),
        "sample_name_to_id": {"s1": 1, "s2": 2},
        "contig_name_to_id": {},
    }
    sample_title = SimpleNamespace(text="")
    controller = AvailabilityController(
        AvailabilityBindings(
            availability_service=AvailabilityService(),
            filtering_pairs=lambda: None,
            original_contigs=["c1", "c2"],
            original_samples=["s1", "s2"],
            sample_scope=SimpleNamespace(active=0),
            widgets=widgets,
            sort_sample_select=_select(),
            update_completions=lambda widget, values: setattr(widget, "options", list(values)),
            total_coverage_count=2,
            filtering_title=SimpleNamespace(text=""),
            contig_title=SimpleNamespace(text=""),
            sample_title=sample_title,
            mag_title=SimpleNamespace(text=""),
        )
    )

    controller.update_titles()

    assert "2 available" in sample_title.text


def test_contig_outside_selected_mag_does_not_override_mag_sample_scope():
    widgets = {
        "has_mags": True,
        "view_radio": SimpleNamespace(active=1),
        "contig_select": _select("c1"),
        "sample_select": _select("s1"),
        "mag_select": _select("m1"),
        "sample_name_to_id": {"s1": 1},
        "contig_name_to_id": {"c1": 1},
        "mag_to_contigs": {"m1": ["mag-child"]},
        "mag_to_sample_ids": {"m1": {1}},
    }
    controller = AvailabilityController(
        AvailabilityBindings(
            availability_service=AvailabilityService(),
            filtering_pairs=lambda: None,
            original_contigs=["c1", "c2"],
            original_samples=["s1", "s2"],
            sample_scope=SimpleNamespace(active=0),
            widgets=widgets,
            sort_sample_select=_select(),
            update_completions=lambda widget, values: setattr(widget, "options", list(values)),
            total_coverage_count=2,
            filtering_title=SimpleNamespace(text=""),
            contig_title=SimpleNamespace(text=""),
            sample_title=SimpleNamespace(text=""),
            mag_title=SimpleNamespace(text=""),
        )
    )

    assert controller.compute_contigs() == ["c2"]
    assert controller.compute_samples() == ["s1"]


@pytest.mark.parametrize("subject_view", [0, 1])
def test_sample_options_and_title_share_selected_mag_scope_without_contig(subject_view):
    widgets = {
        "has_mags": True,
        "view_radio": SimpleNamespace(active=subject_view),
        "contig_select": _select(""),
        "sample_select": _select("", ("GL34_UP", "other")),
        "mag_select": _select("m1", ("m1",)),
        "sample_name_to_id": {"GL34_UP": 1, "other": 2},
        "contig_name_to_id": {"c1": 1},
        "mag_to_contigs": {"m1": ["c1"]},
        "mag_to_sample_ids": {"m1": {1}},
    }
    sample_title = SimpleNamespace(text="")
    controller = AvailabilityController(
        AvailabilityBindings(
            availability_service=FilteredMagAvailabilityService(),
            filtering_pairs=lambda: {
                "sql": "SELECT 1",
                "params": [],
                "count_pairs": 1,
                "count_mag_pairs": 1,
            },
            original_contigs=["c1"],
            original_samples=["GL34_UP", "other"],
            sample_scope=SimpleNamespace(active=0),
            widgets=widgets,
            sort_sample_select=_select(),
            update_completions=lambda widget, values: setattr(widget, "options", list(values)),
            total_coverage_count=2,
            filtering_title=SimpleNamespace(text=""),
            contig_title=SimpleNamespace(text=""),
            sample_title=sample_title,
            mag_title=SimpleNamespace(text=""),
        )
    )

    controller.refresh_samples()
    controller.update_titles()

    assert widgets["sample_select"].options == ["GL34_UP"]
    assert "1 available" in sample_title.text
