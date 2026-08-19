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
    assert catalog.sample_metrics["Sample"][0][0] == "Sample_name"
    assert catalog.mag_metrics["Contig"] == [("Contig_length", "Contig length"), ("Topology", "Topology")]
    assert catalog.mag_categories == ["Contig", "Coverage"]
    assert catalog.sample_contig_categories == ["Sample", "Coverage"]
    assert catalog.sample_mag_categories == ["Sample", "MAG coverage"]
