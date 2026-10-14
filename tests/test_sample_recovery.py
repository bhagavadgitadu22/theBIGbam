import duckdb
from thebigbam.database.calculating_data import _retryable_sample_names


def test_explicit_failed_and_interrupted_states_are_retryable():
    conn = duckdb.connect()
    conn.execute("CREATE TABLE Sample (Sample_id INTEGER, Sample_name VARCHAR)")
    conn.execute("INSERT INTO Sample VALUES (1,'complete'),(2,'failed'),(3,'interrupted'),(4,'legacy')")
    assert _retryable_sample_names(conn) == set()
    conn.execute("CREATE TABLE Sample_processing (Sample_id INTEGER, State VARCHAR)")
    conn.execute("INSERT INTO Sample_processing VALUES (1,'complete'),(2,'failed'),(3,'processing')")
    assert _retryable_sample_names(conn) == {'failed', 'interrupted'}
    conn.execute("UPDATE Sample_processing SET State='complete' WHERE Sample_id=3")
    assert _retryable_sample_names(conn) == {'failed'}
