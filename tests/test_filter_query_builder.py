import duckdb
import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.repositories.filtering import FilteringRepository
from thebigbam.plotting.services.filter_query import FilterQueryBuilder
from thebigbam.plotting.services.filtering import FilterExpressionService


def _database():
    connection = duckdb.connect(":memory:")
    connection.execute("CREATE TABLE Contig (Contig_id INTEGER, Contig_length INTEGER, Contig_name VARCHAR)")
    connection.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    connection.execute("CREATE TABLE Coverage (Contig_id INTEGER, Sample_id INTEGER)")
    connection.execute(
        'CREATE TABLE Contig_annotation_core (Annotation_id INTEGER, Contig_id INTEGER, "Start" INTEGER, "End" INTEGER)'
    )
    connection.execute("CREATE TABLE Annotation_qualifier (Annotation_id INTEGER, Key VARCHAR, Value VARCHAR)")
    connection.execute("INSERT INTO Contig VALUES (1, 1000, 'c1'), (2, 2000, 'c2')")
    connection.execute("INSERT INTO Sample VALUES (1, 's1'), (2, 's2')")
    connection.execute("INSERT INTO Coverage VALUES (1, 1), (2, 1), (2, 2)")
    connection.execute("INSERT INTO Contig_annotation_core VALUES (1, 1, 10, 20), (2, 2, 10, 20)")
    connection.execute("INSERT INTO Annotation_qualifier VALUES (1, 'product', 'capsid^tail'), (2, 'product', 'tail')")
    return connection


def _builder(has_samples=True):
    metadata = {
        "Contig": {"columns": {"Contig_length": {"type": "numeric"}}},
        "Annotations": {
            "columns": {
                "product": {
                    "type": "text",
                    "source": "Annotation_qualifier",
                    "qualifier_key": "product",
                }
            }
        },
    }
    return FilterQueryBuilder(metadata, {"Contig_length": 10}, has_samples=has_samples)


def test_compiler_preserves_encoding_qualifier_tokens_and_intersection():
    expression = FilterExpression(
        (
            FilterSection(
                (
                    FilterPredicate("Contig", "Contig_length", ">", "150"),
                    FilterPredicate("Annotations", "product", "=", "tail"),
                ),
                ("AND",),
            ),
        )
    )
    compiled = _builder().compile(expression)
    assert compiled is not None
    assert compiled.parameters == (1500, "product", "tail")

    connection = _database()
    rows = connection.execute(compiled.sql, compiled.parameters).fetchall()
    assert sorted(rows) == [(2, 1), (2, 2)]
    connection.close()


def test_sections_use_union_and_no_sample_database_projects_null():
    expression = FilterExpression(
        (
            FilterSection((FilterPredicate("Contig", "Contig_length", "=", 100),)),
            FilterSection((FilterPredicate("Contig", "Contig_length", "=", 200),)),
        ),
        ("OR",),
    )
    compiled = _builder(has_samples=False).compile(expression)
    assert compiled is not None
    assert "UNION" in compiled.sql
    assert "NULL" in compiled.sql
    assert compiled.parameters == (1000, 2000)


def test_expression_service_caches_compilation_and_count_queries():
    connection = _database()
    repository = FilteringRepository(connection)
    service = FilterExpressionService(repository, _builder(), has_mags=False)
    expression = FilterExpression((FilterSection((FilterPredicate("Contig", "Contig_length", ">", 100),)),))

    first = service.evaluate(expression)
    second = service.evaluate(expression)

    assert first == second
    assert first is not None
    assert (first.pair_count, first.contig_count, first.sample_count) == (2, 1, 2)
    assert first.compiled.sql.startswith('SELECT Contig_id, Sample_id FROM "_thebigbam_filter_')
    assert first.compiled.parameters == ()
    assert repository.query_count == 1
    service.invalidate()
    service.evaluate(expression)
    assert repository.query_count == 1
    connection.close()


def test_expression_service_reuses_an_older_filter_snapshot():
    connection = _database()
    repository = FilteringRepository(connection)
    service = FilterExpressionService(repository, _builder(), has_mags=False)
    narrow = FilterExpression(
        (FilterSection((FilterPredicate("Contig", "Contig_length", ">", 100),)),)
    )
    broad = FilterExpression(
        (FilterSection((FilterPredicate("Contig", "Contig_length", ">", 50),)),)
    )

    first = service.evaluate(narrow)
    service.invalidate()
    service.evaluate(broad)
    service.invalidate()
    restored = service.evaluate(narrow)

    assert restored == first
    assert repository.query_count == 2
    connection.close()


def test_scaled_metric_filter_uses_compact_raw_table():
    builder = FilterQueryBuilder(
        {"Coverage": {"columns": {"Coverage_mean": {"type": "numeric"}}}},
        {"Coverage_mean": 10000},
        has_samples=True,
    )

    compiled = builder.compile(
        FilterExpression(
            (FilterSection((FilterPredicate("Coverage", "Coverage_mean", ">", 12.5),)),)
        )
    )

    assert compiled is not None
    assert "FROM Coverage v" in compiled.sql
    assert "Explicit_coverage" not in compiled.sql
    assert compiled.parameters == (125000,)


@pytest.mark.parametrize(
    ("category", "column", "metadata", "sql_fragment"),
    (
        ("Sample", "Sample_name", {"type": "text"}, "FROM Sample s"),
        ("MAG", "MAG_name", {"type": "text"}, "FROM MAG mg"),
        ("Coverage", "Coverage_mean", {"type": "numeric"}, "FROM Coverage v"),
        ("MAG coverage", "Mean_coverage", {"type": "numeric"}, "FROM Explicit_coverage_per_MAG v"),
        (
            "Annotations",
            "taxonomy",
            {"type": "text", "source": "Contig_qualifier", "qualifier_key": "taxonomy"},
            "FROM Contig_qualifier cq",
        ),
    ),
)
def test_all_filter_source_families_compile(category, column, metadata, sql_fragment):
    builder = FilterQueryBuilder({category: {"columns": {column: metadata}}}, {}, has_samples=True)
    compiled = builder.compile(FilterExpression((FilterSection((FilterPredicate(category, column, "=", "1"),)),)))
    assert compiled is not None
    assert sql_fragment in compiled.sql


def test_invalid_numeric_predicate_is_skipped_without_breaking_following_connector():
    expression = FilterExpression(
        (
            FilterSection(
                (
                    FilterPredicate("Contig", "Contig_length", ">", "invalid"),
                    FilterPredicate("Contig", "Contig_length", ">", "100"),
                    FilterPredicate("Contig", "Contig_length", "<", "300"),
                ),
                ("OR", "AND"),
            ),
        )
    )
    compiled = _builder().compile(expression)
    assert compiled is not None
    assert "IS NOT DISTINCT FROM" in compiled.sql
    assert "UNION" not in compiled.sql


def test_standalone_contig_scoped_filter_expands_to_available_samples_at_boundary():
    compiled = _builder().compile(
        FilterExpression(
            (FilterSection((FilterPredicate("Annotations", "product", "=", "tail"),)),)
        )
    )

    assert compiled is not None
    assert "LEFT JOIN Coverage _coverage" in compiled.sql
    assert sorted(_database().execute(compiled.sql, compiled.parameters).fetchall()) == [
        (1, 1),
        (2, 1),
        (2, 2),
    ]
