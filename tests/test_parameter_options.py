from thebigbam.plotting.controls.parameter_options import ParameterOptionCatalog, format_column


def test_parameter_option_catalog_derives_stable_cached_options():
    metadata = {
        "Sample": {
            "source": "Sample",
            "columns": {"Read_count": {"type": "numeric"}, "Group": {"type": "text"}},
        },
        "Contig": {
            "source": "Contig",
            "columns": {"Contig_length": {"type": "numeric"}, "Topology": {"type": "text"}},
        },
        "Coverage": {"source": "Coverage", "columns": {"Breadth_percentage": {"type": "numeric"}}},
        "MAG coverage": {"source": "MAG_coverage", "columns": {"Mean": {"type": "numeric"}}},
        "Annotations": {"source": "Contig_annotation", "columns": {}},
    }

    catalog = ParameterOptionCatalog.from_filtering_metadata(metadata)

    assert format_column("Breadth_percentage") == ("Breadth_percentage", "Breadth (%)")
    assert catalog.sample_metrics["Sample"] == catalog.columns["Sample"]
    assert catalog.sample_metrics["Coverage"] == catalog.columns["Coverage"]
    assert catalog.mag_metrics["Contig"] == [("Contig_length", "Contig length"), ("Topology", "Topology")]
    assert catalog.mag_categories == ["Contig", "Coverage"]
    assert catalog.sample_contig_categories == ["Sample", "Coverage", "MAG coverage"]
    assert catalog.sample_mag_categories == ["Sample", "MAG coverage"]


def test_ordering_metrics_match_filtering_columns_for_shared_categories():
    metadata = {
        "Sample": {
            "source": "Sample",
            "columns": {"Sample_name": {"type": "text"}, "Group": {"type": "text"}},
        },
        "Coverage": {
            "source": "Explicit_coverage",
            "columns": {
                "Mean_coverage": {"type": "numeric"},
                "Coverage_class": {"type": "text"},
            },
        },
        "MAG coverage": {
            "source": "Explicit_coverage_per_MAG",
            "columns": {
                "Mean_coverage": {"type": "numeric"},
                "Coverage_class": {"type": "text"},
            },
        },
    }

    catalog = ParameterOptionCatalog.from_filtering_metadata(metadata)

    for category in catalog.sample_contig_categories:
        assert catalog.sample_metrics[category] == catalog.columns[category]
    for category in catalog.sample_mag_categories:
        assert catalog.sample_metrics[category] == catalog.columns[category]
