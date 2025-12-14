//! SQLite database operations.
//!
//! Database schema:
//! - `Contig`: Contig metadata (name, length, annotation tool)
//! - `Sample`: Sample names from BAM files
//! - `Presences`: Coverage percentage per contig per sample
//! - `Contig_annotation`: Gene annotations from GenBank
//! - `Variable`: Feature metadata (name, type, table name)
//! - Feature tables: One table per feature type (coverage, tau, etc.)
//!
//! Each sample writes to a temporary database during parallel processing.
//! Temp databases are merged into the main database sequentially after processing.

use anyhow::{Context, Result};
use rusqlite::{params, Connection};
use std::collections::HashMap;
use std::path::Path;

use crate::types::{feature_table_name, ContigInfo, FeatureAnnotation, FeaturePoint, PresenceData, VARIABLES};

/// Create the main SQLite database with schema and initial data.
/// If `create_indexes` is false, skip creating indexes (for Python compatibility).
pub fn create_metadata_db(
    db_path: &Path,
    contigs: &[ContigInfo],
    annotations: &[FeatureAnnotation],
    create_indexes: bool,
) -> Result<()> {
    let conn = Connection::open(db_path)
        .with_context(|| format!("Failed to create database: {}", db_path.display()))?;

    // Enable WAL mode for better concurrent read performance
    conn.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;")
        .context("Failed to set database pragmas")?;

    // Create core tables
    create_core_tables(&conn, create_indexes)?;

    // Insert contigs
    insert_contigs(&conn, contigs)?;

    // Insert annotations
    insert_annotations(&conn, annotations)?;

    // Insert variables and create feature tables
    create_variable_tables(&conn, create_indexes)?;

    Ok(())
}

/// Create core database tables.
fn create_core_tables(conn: &Connection, create_indexes: bool) -> Result<()> {
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Annotation_tool TEXT
        )",
        [],
    )
    .context("Failed to create Contig table")?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Sample_name TEXT UNIQUE
        )",
        [],
    )
    .context("Failed to create Sample table")?;

    conn.execute(
        "CREATE TABLE Presences (
            Presence_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Sample_id INTEGER,
            Coverage_percentage REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        )",
        [],
    )
    .context("Failed to create Presences table")?;

    if create_indexes {
        conn.execute(
            "CREATE INDEX idx_presences_contig_sample ON Presences(Contig_id, Sample_id)",
            [],
        )
        .context("Failed to create Presences index")?;
    }

    conn.execute(
        "CREATE TABLE Contig_annotation (
            Contig_annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Start INTEGER,
            End INTEGER,
            Strand INTEGER,
            Type TEXT,
            Product TEXT,
            Function TEXT,
            Phrog INTEGER,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id)
        )",
        [],
    )
    .context("Failed to create Contig_annotation table")?;

    conn.execute(
        "CREATE TABLE Variable (
            Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Variable_name TEXT UNIQUE,
            Subplot TEXT,
            Module TEXT,
            Type TEXT,
            Color TEXT,
            Alpha REAL,
            Fill_alpha REAL,
            Size REAL,
            Title TEXT,
            Help TEXT,
            Feature_table_name TEXT
        )",
        [],
    )
    .context("Failed to create Variable table")?;

    conn.execute(
        "CREATE TABLE Summary (
            Summary_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER NOT NULL,
            Sample_id INTEGER NOT NULL,
            Variable_id INTEGER NOT NULL,
            Row_count INTEGER NOT NULL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id),
            FOREIGN KEY(Variable_id) REFERENCES Variable(Variable_id),
            UNIQUE(Contig_id, Sample_id, Variable_id)
        )",
        [],
    )
    .context("Failed to create Summary table")?;

    if create_indexes {
        conn.execute(
            "CREATE INDEX idx_summary_lookup ON Summary(Contig_id, Sample_id, Variable_id)",
            [],
        )
        .context("Failed to create Summary index")?;
    }

    Ok(())
}

/// Insert contigs into the database.
fn insert_contigs(conn: &Connection, contigs: &[ContigInfo]) -> Result<()> {
    let mut stmt = conn
        .prepare("INSERT INTO Contig (Contig_name, Contig_length, Annotation_tool) VALUES (?1, ?2, ?3)")
        .context("Failed to prepare contig insert")?;

    for contig in contigs {
        stmt.execute(params![
            &contig.name,
            contig.length as i64,
            &contig.annotation_tool
        ])
        .with_context(|| format!("Failed to insert contig: {}", contig.name))?;
    }

    Ok(())
}

/// Insert annotations into the database.
fn insert_annotations(conn: &Connection, annotations: &[FeatureAnnotation]) -> Result<()> {
    conn.execute("BEGIN TRANSACTION", [])
        .context("Failed to begin annotation transaction")?;

    let mut stmt = conn
        .prepare(
            "INSERT INTO Contig_annotation (Contig_id, Start, End, Strand, Type, Product, Function, Phrog)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
        )
        .context("Failed to prepare annotation insert")?;

    for ann in annotations {
        stmt.execute(params![
            ann.contig_id,
            ann.start,
            ann.end,
            ann.strand,
            &ann.feature_type,
            &ann.product,
            &ann.function,
            &ann.phrog,
        ])
        .context("Failed to insert annotation")?;
    }

    conn.execute("COMMIT", [])
        .context("Failed to commit annotation transaction")?;

    Ok(())
}

/// Create Variable entries and Feature_* tables.
fn create_variable_tables(conn: &Connection, create_indexes: bool) -> Result<()> {
    for v in VARIABLES {
        let table_name = feature_table_name(v.name);

        conn.execute(
            "INSERT INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            params![
                v.name,
                v.subplot,
                v.module,
                v.plot_type.as_str(),
                v.color,
                v.alpha,
                v.fill_alpha,
                v.size,
                v.title,
                v.help,
                &table_name
            ],
        )
        .with_context(|| format!("Failed to insert variable: {}", v.name))?;

        // Create the Feature_* table
        conn.execute(
            &format!(
                "CREATE TABLE {} (
                    Feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    Contig_id INTEGER,
                    Sample_id INTEGER,
                    First_position INTEGER,
                    Last_position INTEGER,
                    Value REAL,
                    FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
                    FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
                )",
                table_name
            ),
            [],
        )
        .with_context(|| format!("Failed to create table: {}", table_name))?;

        if create_indexes {
            conn.execute(
                &format!(
                    "CREATE INDEX idx_{}_lookup ON {}(Sample_id, Contig_id)",
                    v.name, table_name
                ),
                [],
            )
            .with_context(|| format!("Failed to create index for: {}", table_name))?;
        }
    }

    Ok(())
}

/// Create a temporary per-sample database to avoid lock contention.
pub fn create_temp_sample_db(temp_db_path: &Path) -> Result<Connection> {
    // Remove existing temp DB if it exists
    let _ = std::fs::remove_file(temp_db_path);
    
    let conn = Connection::open(temp_db_path)
        .with_context(|| format!("Failed to create temp DB: {}", temp_db_path.display()))?;

    conn.execute_batch("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;")
        .context("Failed to set temp DB pragmas")?;

    conn.execute(
        "CREATE TABLE TempPresences (
            Contig_name TEXT,
            Sample_name TEXT,
            Coverage_percentage REAL
        )",
        [],
    )
    .context("Failed to create TempPresences table")?;

    conn.execute(
        "CREATE TABLE TempFeatures (
            Variable_name TEXT,
            Contig_name TEXT,
            First_position INTEGER,
            Last_position INTEGER,
            Value REAL
        )",
        [],
    )
    .context("Failed to create TempFeatures table")?;

    Ok(conn)
}

/// Write features to a temp database.
pub fn write_features_to_temp_db(conn: &Connection, features: &[FeaturePoint]) -> Result<()> {
    conn.execute("BEGIN TRANSACTION", [])
        .context("Failed to begin feature transaction")?;

    let mut stmt = conn
        .prepare("INSERT INTO TempFeatures (Variable_name, Contig_name, First_position, Last_position, Value) VALUES (?1, ?2, ?3, ?4, ?5)")
        .context("Failed to prepare feature insert")?;

    for f in features {
        stmt.execute(params![&f.feature, &f.contig_name, f.start_pos, f.end_pos, f.value])
            .context("Failed to insert feature")?;
    }

    conn.execute("COMMIT", [])
        .context("Failed to commit feature transaction")?;

    Ok(())
}

/// Write presences to a temp database.
pub fn write_presences_to_temp_db(
    conn: &Connection,
    sample_name: &str,
    presences: &[PresenceData],
) -> Result<()> {
    let mut stmt = conn
        .prepare("INSERT INTO TempPresences (Contig_name, Sample_name, Coverage_percentage) VALUES (?1, ?2, ?3)")
        .context("Failed to prepare presence insert")?;

    for p in presences {
        stmt.execute(params![&p.contig_name, sample_name, p.coverage_pct])
            .context("Failed to insert presence")?;
    }

    Ok(())
}

/// Merge a temp sample DB into the main database.
pub fn merge_temp_db_into_main(
    main_db_path: &Path,
    temp_db_path: &Path,
    contigs: &[ContigInfo],
) -> Result<()> {
    let conn = Connection::open(main_db_path)
        .with_context(|| format!("Failed to open main DB: {}", main_db_path.display()))?;

    // Build contig name -> id mapping
    let contig_name_to_id: HashMap<String, i64> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
        .collect();

    // Attach temp DB
    let temp_path_str = temp_db_path
        .to_str()
        .context("Invalid temp DB path encoding")?;
    conn.execute("ATTACH DATABASE ?1 AS src", [temp_path_str])
        .context("Failed to attach temp DB")?;

    // Get sample names from temp DB
    let sample_names: Vec<String> = {
        let mut stmt = conn
            .prepare("SELECT DISTINCT Sample_name FROM src.TempPresences")
            .context("Failed to query sample names")?;
        let rows = stmt.query_map([], |row| row.get(0))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    conn.execute("BEGIN TRANSACTION", [])
        .context("Failed to begin merge transaction")?;

    // Insert samples
    for sample_name in &sample_names {
        conn.execute(
            "INSERT OR IGNORE INTO Sample (Sample_name) VALUES (?1)",
            [sample_name],
        )?;
    }

    // Get sample_id mapping
    let sample_name_to_id: HashMap<String, i64> = {
        let mut stmt = conn.prepare("SELECT Sample_name, Sample_id FROM Sample")?;
        let rows = stmt.query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    // Insert presences
    merge_presences(&conn, &contig_name_to_id, &sample_name_to_id)?;

    // Insert features
    merge_features(&conn, &contig_name_to_id, &sample_name_to_id, &sample_names)?;

    conn.execute("COMMIT", [])
        .context("Failed to commit merge transaction")?;
    conn.execute("DETACH DATABASE src", [])
        .context("Failed to detach temp DB")?;

    Ok(())
}

/// Merge presences from temp DB.
fn merge_presences(
    conn: &Connection,
    contig_name_to_id: &HashMap<String, i64>,
    sample_name_to_id: &HashMap<String, i64>,
) -> Result<()> {
    let mut insert_stmt = conn.prepare(
        "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?1, ?2, ?3)",
    )?;

    let presences: Vec<(String, String, f32)> = {
        let mut stmt =
            conn.prepare("SELECT Contig_name, Sample_name, Coverage_percentage FROM src.TempPresences")?;
        let rows = stmt.query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    for (contig_name, sample_name, coverage_pct) in presences {
        if let (Some(&contig_id), Some(&sample_id)) = (
            contig_name_to_id.get(&contig_name),
            sample_name_to_id.get(&sample_name),
        ) {
            insert_stmt.execute(params![contig_id, sample_id, coverage_pct])?;
        }
    }

    Ok(())
}

/// Merge features from temp DB into feature tables.
fn merge_features(
    conn: &Connection,
    contig_name_to_id: &HashMap<String, i64>,
    sample_name_to_id: &HashMap<String, i64>,
    sample_names: &[String],
) -> Result<()> {
    let sample_id = sample_names
        .first()
        .and_then(|name| sample_name_to_id.get(name).copied())
        .unwrap_or(1);

    // Get Variable_id mapping
    let variable_name_to_id: HashMap<String, i64> = {
        let mut stmt = conn.prepare("SELECT Variable_name, Variable_id FROM Variable")?;
        let rows = stmt.query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    let mut summary_stmt = conn.prepare(
        "INSERT INTO Summary (Contig_id, Sample_id, Variable_id, Row_count) VALUES (?1, ?2, ?3, ?4)"
    )?;

    for v in VARIABLES {
        let table_name = feature_table_name(v.name);

        let features: Vec<(String, i32, i32, f32)> = {
            let mut stmt = conn.prepare(
                "SELECT Contig_name, First_position, Last_position, Value FROM src.TempFeatures WHERE Variable_name = ?1",
            )?;
            let rows = stmt.query_map([v.name], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        if features.is_empty() {
            continue;
        }

        // Track row counts per contig for summary
        let mut contig_row_counts: HashMap<i64, usize> = HashMap::new();

        let mut stmt = conn.prepare(&format!(
            "INSERT INTO {} (Contig_id, Sample_id, First_position, Last_position, Value) VALUES (?1, ?2, ?3, ?4, ?5)",
            table_name
        ))?;

        for (contig_name, first_pos, last_pos, value) in features {
            if let Some(&contig_id) = contig_name_to_id.get(&contig_name) {
                stmt.execute(params![contig_id, sample_id, first_pos, last_pos, value])?;
                *contig_row_counts.entry(contig_id).or_insert(0) += 1;
            }
        }

        // Insert summary data for this variable
        if let Some(&variable_id) = variable_name_to_id.get(v.name) {
            for (contig_id, row_count) in contig_row_counts {
                summary_stmt.execute(params![contig_id, sample_id, variable_id, row_count as i64])?;
            }
        }
    }

    Ok(())
}

/// Finalize the database after all samples are processed.
pub fn finalize_db(db_path: &Path) -> Result<()> {
    let conn = Connection::open(db_path)
        .with_context(|| format!("Failed to open DB for finalization: {}", db_path.display()))?;

    conn.execute_batch("PRAGMA wal_checkpoint(TRUNCATE);")
        .context("Failed to checkpoint WAL")?;

    Ok(())
}
