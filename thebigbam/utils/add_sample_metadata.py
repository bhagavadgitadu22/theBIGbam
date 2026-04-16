"""Add sample metadata from CSV as new columns in Sample table.

This module provides functionality to extend the Sample table with
additional metadata columns from a user-provided CSV file.
"""
import csv
import duckdb

from thebigbam.database.database_getters import update_database_metadata
from thebigbam.utils.metadata_helpers import infer_column_type as _infer_column_type, convert_value as _convert_value


def add_add_sample_metadata_args(parser):
    """Define command-line arguments for add-sample-metadata."""
    parser.add_argument('--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--csv', dest='csv_file', required=True, help='CSV file with Sample column and metadata columns. Header should be like: Sample,Var1,Var2,... followed by one row per sample containing the values per variable and per sample.')
    parser.add_argument('--force', action='store_true', default=False, help='Drop existing columns that conflict with CSV columns before adding them')


def run_add_sample_metadata(args):
    """Execute the add-sample-metadata command.
    
    Reads a CSV file and adds its columns (except 'Sample') as new columns
    to the Sample table in the database.
    
    Args:
        args: Namespace with 'db' and 'csv_file' attributes
        
    Returns:
        0 on success, 1 on error
    """
    db_path = args.db
    csv_file = args.csv_file
    
    try:
        # Connect to database (writable)
        conn = duckdb.connect(db_path)
        
        # Check Sample table exists (genbank-only mode has no Sample table)
        if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone() is None:
            print("No Sample table in database (genbank-only mode). Cannot add sample metadata.")
            conn.close()
            return 1

        # Get existing Sample table columns
        result = conn.execute("PRAGMA table_info(Sample)").fetchall()
        existing_columns = {row[1] for row in result}  # Column name is at index 1
        
        # Get existing samples from database
        sample_rows = conn.execute("SELECT Sample_name FROM Sample").fetchall()
        db_samples = {row[0] for row in sample_rows}
        
        # Read CSV file
        with open(csv_file, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            
            # Validate 'Sample' column exists (case-sensitive)
            if 'Sample' not in reader.fieldnames:
                print(f"Error: CSV must have a 'Sample' column (case-sensitive). Found columns: {reader.fieldnames}")
                conn.close()
                return 1
            
            # Get new column names (all except 'Sample')
            new_columns = [col for col in reader.fieldnames if col != 'Sample']
            
            if not new_columns:
                print("Error: CSV has no metadata columns (only 'Sample' column found)")
                conn.close()
                return 1
            
            # Check for column name conflicts
            conflicts = [col for col in new_columns if col in existing_columns]
            if conflicts:
                if args.force:
                    print(f"Overwriting existing columns (--force): {', '.join(conflicts)}")
                else:
                    print(f"Error: The following columns already exist in Sample table: {', '.join(conflicts)}")
                    print("Use --force to overwrite them, or rename these columns in your CSV file.")
                    conn.close()
                    return 1
            
            # Read all rows to infer types and prepare data
            rows = list(reader)
        
        if not rows:
            print("Warning: CSV file has no data rows")
            conn.close()
            return 0
        
        # Collect values per column for type inference
        column_values = {col: [] for col in new_columns}
        for row in rows:
            for col in new_columns:
                column_values[col].append(row.get(col, ''))
        
        # Infer types for each column
        column_types = {col: _infer_column_type(values) for col, values in column_values.items()}
        
        # Add new columns to Sample table (skip columns that already exist)
        for col in new_columns:
            if col in existing_columns:
                continue
            col_type = column_types[col]
            # Quote column name to handle special characters
            conn.execute(f'ALTER TABLE Sample ADD COLUMN "{col}" {col_type}')
            print(f"Added column '{col}' ({col_type})")
        
        # Update rows - match by Sample name
        updated_count = 0
        skipped_samples = []
        
        for row in rows:
            sample_name = row['Sample']
            
            # Check if sample exists in database
            if sample_name not in db_samples:
                skipped_samples.append(sample_name)
                continue
            
            # Build UPDATE statement for this row
            set_clauses = []
            values = []
            for col in new_columns:
                raw_value = row.get(col, '')
                converted = _convert_value(raw_value, column_types[col])
                set_clauses.append(f'"{col}" = ?')
                values.append(converted)
            
            values.append(sample_name)  # For WHERE clause
            
            update_sql = f"UPDATE Sample SET {', '.join(set_clauses)} WHERE Sample_name = ?"
            conn.execute(update_sql, values)
            updated_count += 1
        
        update_database_metadata(conn)
        conn.close()

        # Report results
        print(f"\nUpdated {updated_count} sample(s) with {len(new_columns)} new column(s)")
        
        if skipped_samples:
            print(f"\nWarning: {len(skipped_samples)} sample(s) from CSV not found in database:")
            for s in skipped_samples[:10]:  # Show first 10
                print(f"  - {s}")
            if len(skipped_samples) > 10:
                print(f"  ... and {len(skipped_samples) - 10} more")
        
        return 0
        
    except FileNotFoundError:
        print(f"Error: CSV file not found: {csv_file}")
        return 1
    except duckdb.Error as e:
        print(f"Database error: {e}")
        return 1
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
