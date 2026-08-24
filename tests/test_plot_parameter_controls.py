from types import SimpleNamespace

from thebigbam.plotting.controls.parameter_options import ParameterOptionCatalog
from thebigbam.plotting.controls.plot_parameters import build_plot_parameter_controls


def test_plot_parameter_controls_build_from_catalog_without_database_access():
    metadata = {
        "Sample": {"source": "Sample", "columns": {"Read_count": {"type": "numeric"}}},
        "Contig": {
            "source": "Contig",
            "columns": {"Contig_length": {"type": "numeric"}, "Topology": {"type": "text"}},
        },
        "Coverage": {"source": "Coverage", "columns": {"Mean_coverage": {"type": "numeric"}}},
    }
    catalog = ParameterOptionCatalog.from_filtering_metadata(metadata)

    controls = build_plot_parameter_controls(
        widgets={"has_mags": True, "has_samples": True},
        sample_scope=SimpleNamespace(active=0),
        original_samples=["sample-a"],
        parameter_options=catalog,
        toggle_stylesheet="",
        stylesheet="",
        make_toggle_callback=lambda button, content: lambda event: None,
        compute_sample_completions=lambda query: ["sample-a"],
        push_search_completions=lambda widget, values: None,
        sort_sample_completions=lambda values: values,
        interaction_lock={"locked": False},
    )

    assert controls.mag_category.value == "Contig"
    assert controls.mag_metric.value == "Contig_length"
    assert controls.max_genemap_window.value > 0
    assert controls.max_genemap_window.high >= 10_000_000
    controls.max_genemap_window.value = 10_000_000
    assert controls.max_genemap_window.value == 10_000_000
    assert controls.sample_order_metric.value == "Read_count"
    assert controls.mag_category_sources["Coverage"] == "Coverage"


def test_sample_ordering_includes_text_filtering_metrics():
    metadata = {
        "Sample": {
            "source": "Sample",
            "columns": {"Sample_name": {"type": "text"}, "Group": {"type": "text"}},
        },
        "Coverage": {
            "source": "Coverage",
            "columns": {"Mean_coverage": {"type": "numeric"}, "Coverage_class": {"type": "text"}},
        },
    }
    catalog = ParameterOptionCatalog.from_filtering_metadata(metadata)

    controls = build_plot_parameter_controls(
        widgets={"has_mags": False, "has_samples": True},
        sample_scope=SimpleNamespace(active=1),
        original_samples=["sample-a"],
        parameter_options=catalog,
        toggle_stylesheet="",
        stylesheet="",
        make_toggle_callback=lambda button, content: lambda event: None,
        compute_sample_completions=lambda query: ["sample-a"],
        push_search_completions=lambda widget, values: None,
        sort_sample_completions=lambda values: values,
        interaction_lock={"locked": False},
    )

    controls.sample_order_category.value = "Coverage"

    assert controls.sample_order_metric.options == [
        ("Mean_coverage", "Mean coverage"),
        ("Coverage_class", "Coverage class"),
    ]
