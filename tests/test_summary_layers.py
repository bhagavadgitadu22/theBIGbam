import duckdb

from thebigbam.plotting.repositories.summary import SummaryRepository
from thebigbam.plotting.services.summary import column_scales, decode_map, filter_encode, round_significant


def test_summary_scale_retrieval_and_decoding_are_layered():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Column_scales(Feature_name VARCHAR, Column_name VARCHAR, Scale DOUBLE)")
    connection.execute("INSERT INTO Column_scales VALUES ('Coverage', 'Depth', 10), ('Other', 'Ignored', 2)")

    repository = SummaryRepository(connection)
    assert repository.column_scales() == [("Depth", 10.0)]
    assert repository.query_count == 1
    assert column_scales(connection) == {"Depth": 10.0}
    assert decode_map(connection)["Depth"](123) == 12.3
    assert filter_encode(connection) == {"Depth": 10.0}


def test_summary_numeric_formatting_is_database_free():
    assert round_significant(1234.5, 2) == 1200.0
    assert round_significant(None) is None
