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

use crate::types::{feature_table_name, ContigInfo, FeatureAnnotation, FeaturePoint, PackagingData, PresenceData, VARIABLES};

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
            Coverage_variation REAL,
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
        "CREATE TABLE PhageMechanisms (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Phage_packaging_mechanism TEXT,
            Phage_left_terminus INTEGER,
            Phage_right_terminus INTEGER,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        )",
        [],
    )
    .context("Failed to create PhageMechanisms table")?;

    conn.execute(
        "CREATE TABLE Completeness (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Prevalence_completeness_left REAL,
            Distance_contaminated_left INTEGER,
            Min_missing_left INTEGER,
            Prevalence_completeness_right REAL,
            Distance_contaminated_right INTEGER,
            Min_missing_right INTEGER,
            Total_mismatches REAL,
            Total_deletions REAL,
            Total_insertions REAL,
            Total_reads_clipped REAL,
            Total_reference_clipped REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        )",
        [],
    )
    .context("Failed to create Completeness table")?;

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
    // Features that need statistics columns
    let features_with_stats = ["left_clippings", "right_clippings", "insertions"];

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

        // Skip creating physical table for primary_reads - it will be a VIEW instead
        if v.name == "primary_reads" {
            continue;
        }

        // Create the Feature_* table with or without statistics columns
        let table_sql = if features_with_stats.contains(&v.name) {
            format!(
                "CREATE TABLE {} (
                    Feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    Contig_id INTEGER,
                    Sample_id INTEGER,
                    First_position INTEGER,
                    Last_position INTEGER,
                    Value REAL,
                    Mean REAL,
                    Median REAL,
                    Std REAL,
                    FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
                    FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
                )",
                table_name
            )
        } else {
            format!(
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
            )
        };

        conn.execute(&table_sql, [])
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
            
            // Create additional indexes for strand-specific tables used in VIEW
            if v.name == "primary_reads_plus_only" || v.name == "primary_reads_minus_only" {
                conn.execute(
                    &format!(
                        "CREATE INDEX idx_{}_intervals ON {}(Contig_id, Sample_id, First_position, Last_position)",
                        v.name, table_name
                    ),
                    [],
                )
                .with_context(|| format!("Failed to create intervals index for: {}", table_name))?;
            }
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
            Coverage_percentage REAL,
            Coverage_variation REAL
        )",
        [],
    )
    .context("Failed to create TempPresences table")?;

    // New: TempPhageMechanisms table for per-sample phage mechanism info
    conn.execute(
        "CREATE TABLE TempPhageMechanisms (
            Contig_name TEXT,
            Sample_name TEXT,
            Phage_packaging_mechanism TEXT,
            Phage_left_terminus INTEGER,
            Phage_right_terminus INTEGER
        )",
        [],
    )
    .context("Failed to create TempPhageMechanisms table")?;

    // TempCompleteness table for per-sample completeness statistics
    conn.execute(
        "CREATE TABLE TempCompleteness (
            Contig_name TEXT,
            Sample_name TEXT,
            Prevalence_completeness_left REAL,
            Distance_contaminated_left INTEGER,
            Min_missing_left INTEGER,
            Prevalence_completeness_right REAL,
            Distance_contaminated_right INTEGER,
            Min_missing_right INTEGER,
            Total_mismatches REAL,
            Total_deletions REAL,
            Total_insertions REAL,
            Total_reads_clipped REAL,
            Total_reference_clipped REAL
        )",
        [],
    )
    .context("Failed to create TempCompleteness table")?;

    conn.execute(
        "CREATE TABLE TempFeatures (
            Variable_name TEXT,
            Contig_name TEXT,
            First_position INTEGER,
            Last_position INTEGER,
            Value REAL,
            Mean REAL,
            Median REAL,
            Std REAL
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
        .prepare("INSERT INTO TempFeatures (Variable_name, Contig_name, First_position, Last_position, Value, Mean, Median, Std) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)")
        .context("Failed to prepare feature insert")?;

    for f in features {
        stmt.execute(params![
            &f.feature,
            &f.contig_name,
            f.start_pos,
            f.end_pos,
            f.value,
            f.mean,
            f.median,
            f.std
        ])
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
        .prepare("INSERT INTO TempPresences (Contig_name, Sample_name, Coverage_percentage, Coverage_variation) VALUES (?1, ?2, ?3, ?4)")
        .context("Failed to prepare presence insert")?;

    for p in presences {
        stmt.execute(params![
            &p.contig_name,
            sample_name,
            p.coverage_pct,
            p.coverage_variation
        ])
        .context("Failed to insert presence")?;
    }

    Ok(())
}

/// Write phage packaging data to a temp database.
/// Only called with entries where packaging mechanism was detected (not "No_packaging").
pub fn write_packaging_to_temp_db(
    conn: &rusqlite::Connection,
    sample_name: &str,
    packaging_data: &[PackagingData],
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO TempPhageMechanisms (Contig_name, Sample_name, Phage_packaging_mechanism, Phage_left_terminus, Phage_right_terminus) VALUES (?1, ?2, ?3, ?4, ?5)"
    ).context("Failed to prepare phage packaging insert")?;

    for pkg in packaging_data {
        stmt.execute(params![&pkg.contig_name, sample_name, &pkg.mechanism, pkg.left_terminus, pkg.right_terminus])
            .context("Failed to insert phage packaging")?;
    }

    Ok(())
}

/// Completeness data for a contig.
/// Contains assembly completeness statistics computed from clipping events at reference ends.
/// Individual score components are stored; score_completeness and score_contamination are computed in VIEW.
#[derive(Clone, Debug)]
pub struct CompletenessData {
    pub contig_name: String,
    pub prevalence_left: Option<f64>,
    pub distance_left: Option<i32>,
    pub min_missing_left: Option<i32>,
    pub prevalence_right: Option<f64>,
    pub distance_right: Option<i32>,
    pub min_missing_right: Option<i32>,
    /// Weighted mismatches: Σ(count/coverage) - contributes to both completeness and contamination
    pub total_mismatches: Option<f64>,
    /// Weighted deletions: Σ(count/coverage * median_length) - contributes to contamination
    pub total_deletions: Option<f64>,
    /// Weighted insertions: Σ(count/coverage * median_length) - contributes to completeness
    pub total_insertions: Option<f64>,
    /// Weighted clippings: Σ(count/coverage * median_length) for left+right - contributes to completeness
    pub total_reads_clipped: Option<f64>,
    /// Weighted reference gaps: Σ(avg_prevalence * distance) for paired clips - contributes to contamination
    pub total_reference_clipped: Option<f64>,
}

impl CompletenessData {
    /// Returns true if there is any completeness data.
    pub fn has_data(&self) -> bool {
        self.prevalence_left.is_some() || self.prevalence_right.is_some()
            || self.total_mismatches.is_some() || self.total_deletions.is_some()
            || self.total_insertions.is_some() || self.total_reads_clipped.is_some()
            || self.total_reference_clipped.is_some()
    }
}

/// Write completeness data to a temp database.
/// Only stores entries that have at least one side with data.
pub fn write_completeness_to_temp_db(
    conn: &rusqlite::Connection,
    sample_name: &str,
    completeness_data: &[CompletenessData],
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO TempCompleteness (Contig_name, Sample_name, Prevalence_completeness_left, Distance_contaminated_left, Min_missing_left, Prevalence_completeness_right, Distance_contaminated_right, Min_missing_right, Total_mismatches, Total_deletions, Total_insertions, Total_reads_clipped, Total_reference_clipped) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)"
    ).context("Failed to prepare completeness insert")?;

    for data in completeness_data {
        if data.has_data() {
            stmt.execute(params![
                &data.contig_name,
                sample_name,
                data.prevalence_left,
                data.distance_left,
                data.min_missing_left,
                data.prevalence_right,
                data.distance_right,
                data.min_missing_right,
                data.total_mismatches,
                data.total_deletions,
                data.total_insertions,
                data.total_reads_clipped,
                data.total_reference_clipped,
            ])
            .context("Failed to insert completeness data")?;
        }
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

    // Insert phage mechanisms
    merge_phage_mechanisms(&conn, &contig_name_to_id, &sample_name_to_id)?;

    // Insert completeness data
    merge_completeness(&conn, &contig_name_to_id, &sample_name_to_id)?;

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
        "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage, Coverage_variation) VALUES (?1, ?2, ?3, ?4)"
    )?;

    let presences: Vec<(String, String, f32, f32)> = {
        let mut stmt =
            conn.prepare("SELECT Contig_name, Sample_name, Coverage_percentage, Coverage_variation FROM src.TempPresences")?;
        let rows = stmt.query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    for (contig_name, sample_name, coverage_pct, coverage_variation) in presences {
        if let (Some(&contig_id), Some(&sample_id)) = (
            contig_name_to_id.get(&contig_name),
            sample_name_to_id.get(&sample_name),
        ) {
            insert_stmt.execute(params![contig_id, sample_id, coverage_pct, coverage_variation])?;
        }
    }

    Ok(())
}

/// Merge phage mechanisms from temp DB into PhageMechanisms table.
fn merge_phage_mechanisms(
    conn: &rusqlite::Connection,
    contig_name_to_id: &std::collections::HashMap<String, i64>,
    sample_name_to_id: &std::collections::HashMap<String, i64>,
) -> Result<()> {
    let mut insert_stmt = conn.prepare(
        "INSERT OR IGNORE INTO PhageMechanisms (Contig_id, Sample_id, Phage_packaging_mechanism, Phage_left_terminus, Phage_right_terminus) VALUES (?1, ?2, ?3, ?4, ?5)"
    )?;

    let mechanisms: Vec<(String, String, String, Option<i64>, Option<i64>)> = {
        let mut stmt = conn.prepare(
            "SELECT Contig_name, Sample_name, Phage_packaging_mechanism, Phage_left_terminus, Phage_right_terminus FROM src.TempPhageMechanisms"
        )?;
        let rows = stmt.query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?, row.get(4)?)))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    for (contig_name, sample_name, mechanism, left, right) in mechanisms {
        if let (Some(&contig_id), Some(&sample_id)) = (
            contig_name_to_id.get(&contig_name),
            sample_name_to_id.get(&sample_name),
        ) {
            insert_stmt.execute(params![contig_id, sample_id, mechanism, left, right])?;
        }
    }

    Ok(())
}

/// Merge completeness data from temp DB into Completeness table.
fn merge_completeness(
    conn: &rusqlite::Connection,
    contig_name_to_id: &std::collections::HashMap<String, i64>,
    sample_name_to_id: &std::collections::HashMap<String, i64>,
) -> Result<()> {
    let mut insert_stmt = conn.prepare(
        "INSERT OR IGNORE INTO Completeness (Contig_id, Sample_id, Prevalence_completeness_left, Distance_contaminated_left, Min_missing_left, Prevalence_completeness_right, Distance_contaminated_right, Min_missing_right, Total_mismatches, Total_deletions, Total_insertions, Total_reads_clipped, Total_reference_clipped) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)"
    )?;

    let completeness_rows: Vec<(String, String, Option<f64>, Option<i64>, Option<i64>, Option<f64>, Option<i64>, Option<i64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>)> = {
        let mut stmt = conn.prepare(
            "SELECT Contig_name, Sample_name, Prevalence_completeness_left, Distance_contaminated_left, Min_missing_left, Prevalence_completeness_right, Distance_contaminated_right, Min_missing_right, Total_mismatches, Total_deletions, Total_insertions, Total_reads_clipped, Total_reference_clipped FROM src.TempCompleteness"
        )?;
        let rows = stmt.query_map([], |row| Ok((
            row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?,
            row.get(4)?, row.get(5)?, row.get(6)?, row.get(7)?,
            row.get(8)?, row.get(9)?, row.get(10)?, row.get(11)?, row.get(12)?
        )))?;
        rows.filter_map(|r| r.ok()).collect()
    };

    for (contig_name, sample_name, prev_left, dist_left, min_left, prev_right, dist_right, min_right, total_mm, total_del, total_ins, total_reads_clip, total_ref_clip) in completeness_rows {
        if let (Some(&contig_id), Some(&sample_id)) = (
            contig_name_to_id.get(&contig_name),
            sample_name_to_id.get(&sample_name),
        ) {
            insert_stmt.execute(params![contig_id, sample_id, prev_left, dist_left, min_left, prev_right, dist_right, min_right, total_mm, total_del, total_ins, total_reads_clip, total_ref_clip])?;
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

    // Features that have statistics columns
    let features_with_stats = ["left_clippings", "right_clippings", "insertions"];

    for v in VARIABLES {
        let table_name = feature_table_name(v.name);
        let has_stats = features_with_stats.contains(&v.name);

        let features: Vec<(String, i32, i32, f32, Option<f32>, Option<f32>, Option<f32>)> = {
            let mut stmt = conn.prepare(
                "SELECT Contig_name, First_position, Last_position, Value, Mean, Median, Std FROM src.TempFeatures WHERE Variable_name = ?1",
            )?;
            let rows = stmt.query_map([v.name], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?, row.get(4)?, row.get(5)?, row.get(6)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        if features.is_empty() {
            continue;
        }

        // Track row counts per contig for summary
        let mut contig_row_counts: HashMap<i64, usize> = HashMap::new();

        // Use appropriate INSERT statement based on whether feature has statistics
        if has_stats {
            let mut stmt = conn.prepare(&format!(
                "INSERT INTO {} (Contig_id, Sample_id, First_position, Last_position, Value, Mean, Median, Std) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
                table_name
            ))?;

            for (contig_name, first_pos, last_pos, value, mean, median, std) in features {
                if let Some(&contig_id) = contig_name_to_id.get(&contig_name) {
                    stmt.execute(params![contig_id, sample_id, first_pos, last_pos, value, mean, median, std])?;
                    *contig_row_counts.entry(contig_id).or_insert(0) += 1;
                }
            }
        } else {
            let mut stmt = conn.prepare(&format!(
                "INSERT INTO {} (Contig_id, Sample_id, First_position, Last_position, Value) VALUES (?1, ?2, ?3, ?4, ?5)",
                table_name
            ))?;

            for (contig_name, first_pos, last_pos, value, _, _, _) in features {
                if let Some(&contig_id) = contig_name_to_id.get(&contig_name) {
                    stmt.execute(params![contig_id, sample_id, first_pos, last_pos, value])?;
                    *contig_row_counts.entry(contig_id).or_insert(0) += 1;
                }
            }
        }

        // Insert summary data for this variable
        if let Some(&variable_id) = variable_name_to_id.get(v.name) {
            for (contig_id, row_count) in &contig_row_counts {
                summary_stmt.execute(params![contig_id, sample_id, variable_id, *row_count as i64])?;
            }
        }

        // Special case: also populate Summary for primary_reads VIEW when saving plus_only
        if v.name == "primary_reads_plus_only" {
            if let Some(&primary_reads_variable_id) = variable_name_to_id.get("primary_reads") {
                for (contig_id, row_count) in &contig_row_counts {
                    summary_stmt.execute(params![contig_id, sample_id, primary_reads_variable_id, *row_count as i64])?;
                }
            }
        }
    }

    Ok(())
}

/// Finalize the database after all samples are processed.
pub fn finalize_db(db_path: &Path) -> Result<()> {
    let conn = Connection::open(db_path)
        .with_context(|| format!("Failed to open DB for finalization: {}", db_path.display()))?;

    // Create a VIEW that computes primary_reads as sum of plus and minus strands
    // This properly handles overlapping ranges by collecting boundaries, building minimal segments,
    // and summing values that cover each segment
    conn.execute(
        "CREATE VIEW IF NOT EXISTS Feature_primary_reads AS
         WITH
         boundaries AS (
             SELECT Contig_id, Sample_id, First_position AS pos
             FROM Feature_primary_reads_plus_only
             UNION
             SELECT Contig_id, Sample_id, Last_position + 1
             FROM Feature_primary_reads_plus_only
             UNION
             SELECT Contig_id, Sample_id, First_position
             FROM Feature_primary_reads_minus_only
             UNION
             SELECT Contig_id, Sample_id, Last_position + 1
             FROM Feature_primary_reads_minus_only
         ),
         segments AS (
             SELECT
                 Contig_id,
                 Sample_id,
                 pos AS First_position,
                 LEAD(pos) OVER (
                     PARTITION BY Contig_id, Sample_id
                     ORDER BY pos
                 ) - 1 AS Last_position
             FROM boundaries
         ),
         summed AS (
             SELECT
                 s.Contig_id,
                 s.Sample_id,
                 s.First_position,
                 s.Last_position,
                 COALESCE(SUM(p.Value), 0) +
                 COALESCE(SUM(m.Value), 0) AS Value
             FROM segments s
             LEFT JOIN Feature_primary_reads_plus_only p
                 ON p.Contig_id = s.Contig_id
                AND p.Sample_id = s.Sample_id
                AND p.First_position <= s.First_position
                AND p.Last_position >= s.Last_position
             LEFT JOIN Feature_primary_reads_minus_only m
                 ON m.Contig_id = s.Contig_id
                AND m.Sample_id = s.Sample_id
                AND m.First_position <= s.First_position
                AND m.Last_position >= s.Last_position
             WHERE s.Last_position >= s.First_position
             GROUP BY
                 s.Contig_id,
                 s.Sample_id,
                 s.First_position,
                 s.Last_position
         )
         SELECT
             ROW_NUMBER() OVER (
                 ORDER BY Contig_id, Sample_id, First_position
             ) AS Feature_id,
             Contig_id,
             Sample_id,
             First_position,
             Last_position,
             Value
         FROM summed
         WHERE Value != 0",
        [],
    )
    .context("Failed to create primary_reads VIEW")?;

    // Create a view joining Contig, Sample, and Presences for explicit presence data
    conn.execute(
        "CREATE VIEW IF NOT EXISTS Explicit_presences AS
         SELECT
             p.Presence_id,
             c.Contig_name,
             s.Sample_name,
             p.Coverage_percentage,
             p.Coverage_variation
         FROM Presences p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_presences VIEW")?;

    // Create a view joining Contig, Sample, Presences, and PhageMechanisms for explicit metadata
    conn.execute(
        "CREATE VIEW IF NOT EXISTS Explicit_phage_mechanisms AS
         SELECT
             p.Presence_id,
             c.Contig_name,
             s.Sample_name,
             p.Coverage_percentage,
             m.Phage_packaging_mechanism,
             m.Phage_left_terminus,
             m.Phage_right_terminus
         FROM Presences p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id
         LEFT JOIN PhageMechanisms m ON p.Contig_id = m.Contig_id AND p.Sample_id = m.Sample_id",
        [],
    )
    .context("Failed to create Explicit_phage_mechanisms VIEW")?;

    // Create a view joining Contig, Sample, and Completeness for explicit completeness data
    // Score_completeness = mismatches + insertions + reads_clipped (sequence in reads not in reference)
    // Score_contamination = mismatches + deletions + reference_clipped (sequence in reference not in reads)
    conn.execute(
        "CREATE VIEW IF NOT EXISTS Explicit_completeness AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             comp.Prevalence_completeness_left,
             comp.Distance_contaminated_left,
             comp.Min_missing_left,
             comp.Prevalence_completeness_right,
             comp.Distance_contaminated_right,
             comp.Min_missing_right,
             comp.Total_mismatches,
             comp.Total_deletions,
             comp.Total_insertions,
             comp.Total_reads_clipped,
             comp.Total_reference_clipped,
             COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_insertions, 0) + COALESCE(comp.Total_reads_clipped, 0) AS Score_completeness,
             CASE WHEN c.Contig_length > 0 THEN 100 - ((COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_insertions, 0) + COALESCE(comp.Total_reads_clipped, 0)) * 100.0 / c.Contig_length) ELSE NULL END AS Percentage_completeness,
             COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_deletions, 0) + COALESCE(comp.Total_reference_clipped, 0) AS Score_contamination,
             CASE WHEN c.Contig_length > 0 THEN (COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_deletions, 0) + COALESCE(comp.Total_reference_clipped, 0)) * 100.0 / c.Contig_length ELSE NULL END AS Percentage_contamination
         FROM Completeness comp
         JOIN Contig c ON comp.Contig_id = c.Contig_id
         JOIN Sample s ON comp.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_completeness VIEW")?;

    conn.execute_batch("PRAGMA wal_checkpoint(TRUNCATE);")
        .context("Failed to checkpoint WAL")?;

    Ok(())
}
