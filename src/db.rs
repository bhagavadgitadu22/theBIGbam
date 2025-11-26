//! SQLite database operations.
//!
//! This module handles creating and updating the SQLite database.
//! Uses WAL mode and temp-DB-per-sample pattern to avoid lock contention.

use anyhow::Result;
use rusqlite::{Connection, params};
use std::collections::HashMap;
use std::path::Path;

use crate::types::{ContigInfo, FeatureAnnotation, FeaturePoint, PlotType, PresenceData, VARIABLES};

/// Create the main SQLite database with schema and initial data.
pub fn create_metadata_db(
    db_path: &Path,
    contigs: &[ContigInfo],
    annotations: &[FeatureAnnotation],
) -> Result<()> {
    let conn = Connection::open(db_path)?;

    // Enable WAL mode for better concurrent read performance
    conn.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;")?;

    // Create core tables
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Annotation_tool TEXT
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Sample_name TEXT UNIQUE
        )",
        [],
    )?;

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
    )?;

    conn.execute(
        "CREATE INDEX idx_presences_contig_sample ON Presences(Contig_id, Sample_id)",
        [],
    )?;

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
    )?;

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
    )?;

    // Insert contigs
    {
        let mut stmt = conn.prepare(
            "INSERT INTO Contig (Contig_name, Contig_length, Annotation_tool) VALUES (?1, ?2, ?3)"
        )?;
        for contig in contigs {
            stmt.execute(params![&contig.name, contig.length as i64, &contig.annotation_tool])?;
        }
    }

    // Insert annotations in a transaction
    conn.execute("BEGIN TRANSACTION", [])?;
    {
        let mut stmt = conn.prepare(
            "INSERT INTO Contig_annotation (Contig_id, Start, End, Strand, Type, Product, Function, Phrog)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)"
        )?;
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
            ])?;
        }
    }
    conn.execute("COMMIT", [])?;

    // Insert variables and create Feature_* tables
    for v in VARIABLES {
        let vtype = match v.plot_type { PlotType::Curve => "curve", PlotType::Bars => "bars" };
        let table_name = format!("Feature_{}", v.name);

        conn.execute(
            "INSERT INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            params![v.name, v.subplot, v.module, vtype, v.color, v.alpha, v.fill_alpha, v.size, v.title, "", &table_name],
        )?;

        // Create the Feature_* table for this variable
        conn.execute(
            &format!(
                "CREATE TABLE {} (
                    Feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    Contig_id INTEGER,
                    Sample_id INTEGER,
                    Position INTEGER,
                    Value REAL,
                    FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
                    FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
                )", table_name
            ),
            [],
        )?;

        // Create index for fast lookups
        conn.execute(
            &format!(
                "CREATE INDEX idx_{}_lookup ON {}(Sample_id, Contig_id)",
                v.name, table_name
            ),
            [],
        )?;
    }

    Ok(())
}

/// Create a temporary per-sample database to avoid lock contention.
/// Each parallel worker writes to its own temp DB.
pub fn create_temp_sample_db(temp_db_path: &Path) -> Result<Connection> {
    let conn = Connection::open(temp_db_path)?;

    conn.execute_batch("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;")?;

    conn.execute(
        "CREATE TABLE TempPresences (
            Contig_name TEXT,
            Sample_name TEXT,
            Coverage_percentage REAL
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE TempFeatures (
            Variable_name TEXT,
            Contig_name TEXT,
            Position INTEGER,
            Value REAL
        )",
        [],
    )?;

    Ok(conn)
}

/// Write features to a temp database (no locking issues).
pub fn write_features_to_temp_db(
    conn: &Connection,
    features: &[FeaturePoint],
) -> Result<()> {
    conn.execute("BEGIN TRANSACTION", [])?;
    {
        let mut stmt = conn.prepare(
            "INSERT INTO TempFeatures (Variable_name, Contig_name, Position, Value) VALUES (?1, ?2, ?3, ?4)"
        )?;
        for f in features {
            stmt.execute(params![&f.feature, &f.contig_name, f.position, f.value])?;
        }
    }
    conn.execute("COMMIT", [])?;
    Ok(())
}

/// Write presences to a temp database.
pub fn write_presences_to_temp_db(
    conn: &Connection,
    sample_name: &str,
    presences: &[PresenceData],
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO TempPresences (Contig_name, Sample_name, Coverage_percentage) VALUES (?1, ?2, ?3)"
    )?;
    for p in presences {
        stmt.execute(params![&p.contig_name, sample_name, p.coverage_pct])?;
    }
    Ok(())
}

/// Merge a temp sample DB into the main database.
/// This is called sequentially after all parallel workers complete.
pub fn merge_temp_db_into_main(
    main_db_path: &Path,
    temp_db_path: &Path,
    contigs: &[ContigInfo],
) -> Result<()> {
    let conn = Connection::open(main_db_path)?;

    // Build contig name -> id mapping
    let contig_name_to_id: HashMap<String, i64> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
        .collect();

    // Attach temp DB
    conn.execute("ATTACH DATABASE ?1 AS src", [temp_db_path.to_str().unwrap()])?;

    // Get distinct sample names from temp DB and insert into main
    let sample_names: Vec<String> = {
        let mut stmt = conn.prepare("SELECT DISTINCT Sample_name FROM src.TempPresences")?;
        let rows = stmt.query_map([], |row| row.get(0))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    conn.execute("BEGIN TRANSACTION", [])?;

    for sample_name in &sample_names {
        // Insert sample (get or create)
        conn.execute("INSERT OR IGNORE INTO Sample (Sample_name) VALUES (?1)", [sample_name])?;
    }

    // Get sample_id mapping
    let sample_name_to_id: HashMap<String, i64> = {
        let mut stmt = conn.prepare("SELECT Sample_name, Sample_id FROM Sample")?;
        let rows = stmt.query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    // Insert presences
    {
        let mut stmt = conn.prepare(
            "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?1, ?2, ?3)"
        )?;

        let presences: Vec<(String, String, f32)> = {
            let mut pstmt = conn.prepare("SELECT Contig_name, Sample_name, Coverage_percentage FROM src.TempPresences")?;
            let rows = pstmt.query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        for (contig_name, sample_name, coverage_pct) in presences {
            if let (Some(&contig_id), Some(&sample_id)) = (
                contig_name_to_id.get(&contig_name),
                sample_name_to_id.get(&sample_name),
            ) {
                stmt.execute(params![contig_id, sample_id, coverage_pct])?;
            }
        }
    }

    // Insert features into their respective Feature_* tables
    for v in VARIABLES {
        let table_name = format!("Feature_{}", v.name);

        // Get features for this variable from temp DB
        let features: Vec<(String, i32, f32)> = {
            let mut stmt = conn.prepare(
                "SELECT Contig_name, Position, Value FROM src.TempFeatures WHERE Variable_name = ?1"
            )?;
            let rows = stmt.query_map([v.name], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        if features.is_empty() {
            continue;
        }

        // Insert into feature table - need sample_id from the presences
        // All features in this temp DB belong to the same sample
        let sample_id = sample_name_to_id.get(&sample_names[0]).copied().unwrap_or(1);

        let mut stmt = conn.prepare(
            &format!("INSERT INTO {} (Contig_id, Sample_id, Position, Value) VALUES (?1, ?2, ?3, ?4)", table_name)
        )?;

        for (contig_name, position, value) in features {
            if let Some(&contig_id) = contig_name_to_id.get(&contig_name) {
                stmt.execute(params![contig_id, sample_id, position, value])?;
            }
        }
    }

    conn.execute("COMMIT", [])?;
    conn.execute("DETACH DATABASE src", [])?;

    Ok(())
}

/// Finalize the database after all samples are processed.
/// Runs VACUUM and switches back to normal journal mode for reads.
pub fn finalize_db(db_path: &Path) -> Result<()> {
    let conn = Connection::open(db_path)?;
    conn.execute_batch("PRAGMA wal_checkpoint(TRUNCATE);")?;
    Ok(())
}
