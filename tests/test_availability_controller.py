from types import SimpleNamespace

from thebigbam.plotting.application.availability import AvailabilityBindings, AvailabilityController


class AvailabilityService:
    def contigs_for_sample(self, sample_id, search, preserve):
        return ("c2",)

    def samples_for_contig(self, contig_id, search):
        return ("s2",)

    def count_contigs_for_sample(self, sample_id):
        return 1

    def count_samples_for_contig(self, contig_id):
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


def test_mag_capable_database_uses_contig_constraints_in_contig_view():
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
    assert controller.compute_samples() == ["s2"]
