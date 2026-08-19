from types import SimpleNamespace

from thebigbam.plotting.downloads.inspect_command import InspectCommandBindings, InspectCommandController


def _group(labels, active):
    return SimpleNamespace(labels=labels, active=active)


def _controller(*, all_samples=False, mag_view=False):
    widgets = {
        "has_samples": True,
        "has_mags": True,
        "view_radio": SimpleNamespace(active=0 if mag_view else 1),
        "mag_select": SimpleNamespace(value="MAG-1"),
        "contig_select": SimpleNamespace(value="contig-1"),
        "sample_select": SimpleNamespace(value="sample-a"),
        "contig_lengths": {"contig-1": 1000, "contig-2": 500},
        "mag_to_contigs": {"MAG-1": ["contig-1", "contig-2"]},
        "variables_widgets_one": [_group(["Coverage"], [0])],
        "variables_widgets_all": [_group(["Coverage"], [0])],
    }
    return InspectCommandController(
        InspectCommandBindings(
            db_path="example.db",
            widgets=widgets,
            sample_scope=SimpleNamespace(active=1 if all_samples else 0),
            combined_features=None,
            subplot_variables={"Coverage": ["Mean_coverage"]},
            from_position=SimpleNamespace(value="101"),
            to_position=SimpleNamespace(value="900"),
            filtered_samples=lambda contig: ["sample-a", "sample-b"],
            feature_classifier=lambda feature: False,
        )
    )


def test_contig_command_uses_filtered_all_sample_set_and_region():
    commands = _controller(all_samples=True).commands()
    assert commands == [
        "thebigbam inspect -d example.db --contig contig-1 --sample sample-a,sample-b "
        "--feature Mean_coverage --region 101-900 > output.tsv"
    ]


def test_mag_command_uses_mag_length_for_region_clamping():
    controller = _controller(mag_view=True)
    controller.bindings.to_position.value = "2000"
    commands = controller.commands()
    assert commands[0].startswith("thebigbam inspect -d example.db --mag MAG-1 --sample sample-a")
    assert "--region 101-1500" in commands[0]
