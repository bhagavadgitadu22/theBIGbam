from types import SimpleNamespace

from thebigbam.plotting.controls.filter_rows import FilterRowFactory


class Visualizations:
    def build_numeric_histogram(self, *args, **kwargs):
        return None

    def build_text_treemap(self, *args, **kwargs):
        return None


def test_filter_row_factory_builds_row_from_cached_column_options():
    metadata = {
        "Contig": {
            "columns": {"Contig_length": {"type": "numeric", "is_bool": False}},
        }
    }
    factory = FilterRowFactory(
        metadata_service=SimpleNamespace(),
        filtering_metadata=metadata,
        columns={"Contig": [("Contig_length", "Contig length")]},
        raw_columns={"Contig": ["Contig_length"]},
        visualizations=Visualizations(),
        refresh=lambda: None,
        stylesheet="",
        enable_timing=False,
    )
    factory.attach_controller(SimpleNamespace(count_rows=lambda: 1, sections=[]))

    result = factory.create_row({}, initial_category="Contig", initial_column="Contig_length", initial_operator=">")

    assert result["category_select"].value == "Contig"
    assert result["subcategory_select"].value == "Contig_length"
    assert result["comparison_select"].value == ">"
    assert result["input_ref"]["widget"].value is None
