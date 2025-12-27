//! DuckDB database operations.
//!
//! Database schema:
//! - `Contig`: Contig metadata (name, length, annotation tool)
//! - `Sample`: Sample names from BAM files
//! - `Presences`: Coverage percentage per contig per sample
//! - `Contig_annotation`: Gene annotations from GenBank
//! - `Variable`: Feature metadata (name, type, table name)
//! - Feature tables: One table per feature type (coverage, tau, etc.)
//!
//! DuckDB advantages over SQLite:
//! - Columnar storage with automatic compression
//! - No need for explicit indexes (uses zone maps)
//! - Better write performance with batch inserts
//! - Simpler concurrency model (single writer, multiple readers)

use anyhow::{Context, Result};
use duckdb::{params, Connection};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Mutex;

use crate::types::{feature_table_name, ContigInfo, FeatureAnnotation, FeaturePoint, PackagingData, PresenceData, VARIABLES};

/// Thread-safe database connection wrapper for sequential writes.
pub struct DbWriter {
    conn: Mutex<Connection>,
    contig_name_to_id: HashMap<String, i64>,
    sample_name_to_id: Mutex<HashMap<String, i64>>,
    next_sample_id: Mutex<i64>,
    variable_name_to_id: HashMap<String, i64>,
}

impl DbWriter {
    /// Create a new database and return a writer for sequential sample insertion.
    pub fn create(
        db_path: &Path,
        contigs: &[ContigInfo],
        annotations: &[FeatureAnnotation],
    ) -> Result<Self> {
        // Remove existing database if present
        let _ = std::fs::remove_file(db_path);

        let conn = Connection::open(db_path)
            .with_context(|| format!("Failed to create database: {}", db_path.display()))?;

        // Create tables
        create_core_tables(&conn)?;
        create_variable_tables(&conn)?;

        // Insert contigs
        insert_contigs(&conn, contigs)?;

        // Insert annotations
        insert_annotations(&conn, annotations)?;

        // Build contig name -> id mapping
        let contig_name_to_id: HashMap<String, i64> = contigs
            .iter()
            .enumerate()
            .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
            .collect();

        // Build variable name -> id mapping
        let variable_name_to_id: HashMap<String, i64> = {
            let mut stmt = conn.prepare("SELECT Variable_name, Variable_id FROM Variable")?;
            let rows = stmt.query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        Ok(Self {
            conn: Mutex::new(conn),
            contig_name_to_id,
            sample_name_to_id: Mutex::new(HashMap::new()),
            next_sample_id: Mutex::new(1),
            variable_name_to_id,
        })
    }

    /// Insert a sample and return its ID.
    pub fn insert_sample(&self, sample_name: &str) -> Result<i64> {
        // Get next sample ID
        let sample_id = {
            let mut next_id = self.next_sample_id.lock().unwrap();
            let id = *next_id;
            *next_id += 1;
            id
        };

        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO Sample (Sample_id, Sample_name) VALUES (?1, ?2)",
            params![sample_id, sample_name],
        )?;
        drop(conn);

        let mut sample_map = self.sample_name_to_id.lock().unwrap();
        sample_map.insert(sample_name.to_string(), sample_id);

        Ok(sample_id)
    }

    /// Write all data for a sample (presences, packaging, completeness, features).
    pub fn write_sample_data(
        &self,
        sample_name: &str,
        presences: &[PresenceData],
        packaging: &[PackagingData],
        completeness: &[CompletenessData],
        features: &[FeaturePoint],
    ) -> Result<()> {
        let sample_id = {
            let sample_map = self.sample_name_to_id.lock().unwrap();
            *sample_map.get(sample_name).context("Sample not found")?
        };

        let conn = self.conn.lock().unwrap();

        // Appenders handle their own transactions - no explicit BEGIN/COMMIT needed
        // Insert presences
        self.write_presences(&conn, sample_id, presences)?;

        // Insert packaging mechanisms
        self.write_packaging(&conn, sample_id, packaging)?;

        // Insert completeness data
        self.write_completeness(&conn, sample_id, completeness)?;

        // Insert features
        self.write_features(&conn, sample_id, features)?;

        Ok(())
    }

    fn write_presences(&self, conn: &Connection, sample_id: i64, presences: &[PresenceData]) -> Result<()> {
        let mut appender = conn.appender("Presences")
            .context("Failed to create Presences appender")?;

        for p in presences {
            if let Some(&contig_id) = self.contig_name_to_id.get(&p.contig_name) {
                // Store coverage_percentage and coverage_variation as integers (×100)
                let cov_pct = (p.coverage_pct * 100.0).round() as i32;
                let cov_var = (p.coverage_variation * 100.0).round() as i32;
                appender.append_row(params![contig_id, sample_id, cov_pct, cov_var])?;
            }
        }

        appender.flush().context("Failed to flush Presences appender")?;
        Ok(())
    }

    fn write_packaging(&self, conn: &Connection, sample_id: i64, packaging: &[PackagingData]) -> Result<()> {
        let mut appender = conn.appender("PhageMechanisms")
            .context("Failed to create PhageMechanisms appender")?;

        for pkg in packaging {
            if let Some(&contig_id) = self.contig_name_to_id.get(&pkg.contig_name) {
                appender.append_row(params![contig_id, sample_id, &pkg.mechanism, pkg.left_terminus, pkg.right_terminus])?;
            }
        }

        appender.flush().context("Failed to flush PhageMechanisms appender")?;
        Ok(())
    }

    fn write_completeness(&self, conn: &Connection, sample_id: i64, completeness: &[CompletenessData]) -> Result<()> {
        let mut appender = conn.appender("Completeness")
            .context("Failed to create Completeness appender")?;

        for data in completeness {
            if data.has_data() {
                if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                    // Store prevalences as integers (×100)
                    let prev_left = data.prevalence_left.map(|v| (v * 100.0).round() as i32);
                    let prev_right = data.prevalence_right.map(|v| (v * 100.0).round() as i32);
                    // Store totals as integers (already counts or small values)
                    let total_mm = data.total_mismatches.map(|v| v.round() as i32);
                    let total_del = data.total_deletions.map(|v| v.round() as i32);
                    let total_ins = data.total_insertions.map(|v| v.round() as i32);
                    let total_rc = data.total_reads_clipped.map(|v| v.round() as i32);
                    let total_ref = data.total_reference_clipped.map(|v| v.round() as i32);

                    appender.append_row(params![
                        contig_id, sample_id,
                        prev_left, data.distance_left, data.min_missing_left,
                        prev_right, data.distance_right, data.min_missing_right,
                        total_mm, total_del, total_ins, total_rc, total_ref
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Completeness appender")?;
        Ok(())
    }

    fn write_features(&self, conn: &Connection, sample_id: i64, features: &[FeaturePoint]) -> Result<()> {
        // Features that have statistics columns
        let features_with_stats = ["left_clippings", "right_clippings", "insertions"];

        // Group features by variable name for efficient batch inserts
        let mut by_variable: HashMap<&str, Vec<&FeaturePoint>> = HashMap::new();
        for f in features {
            by_variable.entry(&f.feature).or_default().push(f);
        }

        // Track row counts for summary table
        let mut summary_data: Vec<(i64, i64, i64)> = Vec::new(); // (contig_id, variable_id, count)

        for v in VARIABLES {
            // Skip primary_reads - it's a VIEW
            if v.name == "primary_reads" {
                continue;
            }

            let Some(feature_points) = by_variable.get(v.name) else {
                continue;
            };

            let table_name = feature_table_name(v.name);
            let has_stats = features_with_stats.contains(&v.name);

            // Determine if this variable stores scaled integers
            let is_scaled = v.name == "tau" || v.name == "mapq";

            // Count rows per contig
            let mut contig_row_counts: HashMap<i64, i64> = HashMap::new();

            // Use Appender for bulk loading
            let mut appender = conn.appender(&table_name)
                .with_context(|| format!("Failed to create appender for {}", table_name))?;

            if has_stats {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let mean = f.mean.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let median = f.median.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let std = f.std.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        appender.append_row(params![contig_id, sample_id, f.start_pos, f.end_pos, value, mean, median, std])?;
                        *contig_row_counts.entry(contig_id).or_insert(0) += 1;
                    }
                }
            } else {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, f.end_pos, value])?;
                        *contig_row_counts.entry(contig_id).or_insert(0) += 1;
                    }
                }
            }

            appender.flush()
                .with_context(|| format!("Failed to flush appender for {}", table_name))?;

            // Collect summary data
            if let Some(&variable_id) = self.variable_name_to_id.get(v.name) {
                for (&contig_id, &count) in &contig_row_counts {
                    summary_data.push((contig_id, variable_id, count));
                }

                // For primary_reads VIEW, use plus_only counts
                if v.name == "primary_reads_plus_only" {
                    if let Some(&primary_reads_var_id) = self.variable_name_to_id.get("primary_reads") {
                        for (&contig_id, &count) in &contig_row_counts {
                            summary_data.push((contig_id, primary_reads_var_id, count));
                        }
                    }
                }
            }
        }

        // Insert summary data using appender
        let mut summary_appender = conn.appender("Summary")
            .context("Failed to create Summary appender")?;
        for (contig_id, variable_id, count) in summary_data {
            summary_appender.append_row(params![contig_id, sample_id, variable_id, count])?;
        }
        summary_appender.flush().context("Failed to flush Summary appender")?;

        Ok(())
    }

    /// Finalize the database after all samples are processed.
    pub fn finalize(self) -> Result<()> {
        let conn = self.conn.into_inner().unwrap();
        create_views(&conn)?;

        // Force checkpoint to compress data and write to disk
        conn.execute("CHECKPOINT", [])
            .context("Failed to checkpoint database")?;

        Ok(())
    }

    /// Get the contig ID for a contig name.
    pub fn get_contig_id(&self, name: &str) -> Option<i64> {
        self.contig_name_to_id.get(name).copied()
    }
}

/// Create core database tables (no indexes - DuckDB uses zone maps).
fn create_core_tables(conn: &Connection) -> Result<()> {
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Annotation_tool TEXT
        )",
        [],
    )
    .context("Failed to create Contig table")?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY,
            Sample_name TEXT UNIQUE
        )",
        [],
    )
    .context("Failed to create Sample table")?;

    // Coverage_percentage and Coverage_variation stored as INTEGER (×100)
    // No Presence_id needed - (Contig_id, Sample_id) is the natural key
    conn.execute(
        "CREATE TABLE Presences (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Coverage_percentage INTEGER,
            Coverage_variation INTEGER,
            PRIMARY KEY (Contig_id, Sample_id)
        )",
        [],
    )
    .context("Failed to create Presences table")?;

    conn.execute(
        "CREATE TABLE PhageMechanisms (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Phage_packaging_mechanism TEXT,
            Phage_left_terminus INTEGER,
            Phage_right_terminus INTEGER
        )",
        [],
    )
    .context("Failed to create PhageMechanisms table")?;

    // Completeness table - prevalences stored as INTEGER (×100), totals as INTEGER
    conn.execute(
        "CREATE TABLE Completeness (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Prevalence_completeness_left INTEGER,
            Distance_contaminated_left INTEGER,
            Min_missing_left INTEGER,
            Prevalence_completeness_right INTEGER,
            Distance_contaminated_right INTEGER,
            Min_missing_right INTEGER,
            Total_mismatches INTEGER,
            Total_deletions INTEGER,
            Total_insertions INTEGER,
            Total_reads_clipped INTEGER,
            Total_reference_clipped INTEGER
        )",
        [],
    )
    .context("Failed to create Completeness table")?;

    // No auto-increment ID needed
    conn.execute(
        "CREATE TABLE Contig_annotation (
            Contig_id INTEGER,
            \"Start\" INTEGER,
            \"End\" INTEGER,
            Strand INTEGER,
            \"Type\" TEXT,
            Product TEXT,
            \"Function\" TEXT,
            Phrog INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_annotation table")?;

    conn.execute(
        "CREATE TABLE Variable (
            Variable_id INTEGER PRIMARY KEY,
            Variable_name TEXT UNIQUE,
            Subplot TEXT,
            Module TEXT,
            \"Type\" TEXT,
            Color TEXT,
            Alpha REAL,
            Fill_alpha REAL,
            \"Size\" REAL,
            Title TEXT,
            Help TEXT,
            Feature_table_name TEXT
        )",
        [],
    )
    .context("Failed to create Variable table")?;

    // (Contig_id, Sample_id, Variable_id) is the natural key
    conn.execute(
        "CREATE TABLE Summary (
            Contig_id INTEGER NOT NULL,
            Sample_id INTEGER NOT NULL,
            Variable_id INTEGER NOT NULL,
            Row_count INTEGER NOT NULL,
            PRIMARY KEY (Contig_id, Sample_id, Variable_id)
        )",
        [],
    )
    .context("Failed to create Summary table")?;

    Ok(())
}

/// Insert contigs into the database using Appender for bulk loading.
fn insert_contigs(conn: &Connection, contigs: &[ContigInfo]) -> Result<()> {
    let mut appender = conn.appender("Contig")
        .context("Failed to create Contig appender")?;

    for (i, contig) in contigs.iter().enumerate() {
        appender.append_row(params![
            (i + 1) as i64,
            &contig.name,
            contig.length as i64,
            &contig.annotation_tool
        ])
        .with_context(|| format!("Failed to append contig: {}", contig.name))?;
    }

    appender.flush().context("Failed to flush Contig appender")?;
    Ok(())
}

/// Insert annotations into the database using Appender for bulk loading.
fn insert_annotations(conn: &Connection, annotations: &[FeatureAnnotation]) -> Result<()> {
    let mut appender = conn.appender("Contig_annotation")
        .context("Failed to create Contig_annotation appender")?;

    for ann in annotations {
        appender.append_row(params![
            ann.contig_id,
            ann.start,
            ann.end,
            ann.strand,
            &ann.feature_type,
            &ann.product,
            &ann.function,
            &ann.phrog,
        ])
        .with_context(|| format!(
            "Failed to append annotation: contig_id={}, start={}, end={}, type={}",
            ann.contig_id, ann.start, ann.end, ann.feature_type
        ))?;
    }

    appender.flush().context("Failed to flush Contig_annotation appender")?;
    Ok(())
}

/// Create Variable entries and Feature_* tables.
fn create_variable_tables(conn: &Connection) -> Result<()> {
    // Features that need statistics columns
    let features_with_stats = ["left_clippings", "right_clippings", "insertions"];

    for (i, v) in VARIABLES.iter().enumerate() {
        let table_name = feature_table_name(v.name);

        conn.execute(
            "INSERT INTO Variable (Variable_id, Variable_name, Subplot, Module, \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12)",
            params![
                (i + 1) as i64,
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

        // All feature values stored as INTEGER (tau and mapq are ×100)
        // No Feature_id needed - data is accessed by (Contig_id, Sample_id, position)
        let table_sql = if features_with_stats.contains(&v.name) {
            format!(
                "CREATE TABLE {} (
                    Contig_id INTEGER,
                    Sample_id INTEGER,
                    First_position INTEGER,
                    Last_position INTEGER,
                    Value INTEGER,
                    Mean INTEGER,
                    Median INTEGER,
                    Std INTEGER
                )",
                table_name
            )
        } else {
            format!(
                "CREATE TABLE {} (
                    Contig_id INTEGER,
                    Sample_id INTEGER,
                    First_position INTEGER,
                    Last_position INTEGER,
                    Value INTEGER
                )",
                table_name
            )
        };

        conn.execute(&table_sql, [])
            .with_context(|| format!("Failed to create table: {}", table_name))?;
    }

    Ok(())
}

/// Create views after all data is inserted.
fn create_views(conn: &Connection) -> Result<()> {
    // Create a VIEW that computes primary_reads as sum of plus and minus strands
    conn.execute(
        "CREATE VIEW Feature_primary_reads AS
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

    // Explicit_presences VIEW - note: divide by 100.0 to get original values
    conn.execute(
        "CREATE VIEW Explicit_presences AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             p.Coverage_percentage / 100.0 AS Coverage_percentage,
             p.Coverage_variation / 100.0 AS Coverage_variation
         FROM Presences p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_presences VIEW")?;

    // Explicit_phage_mechanisms VIEW
    conn.execute(
        "CREATE VIEW Explicit_phage_mechanisms AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             p.Coverage_percentage / 100.0 AS Coverage_percentage,
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

    // Explicit_completeness VIEW - divide prevalences by 100.0
    conn.execute(
        "CREATE VIEW Explicit_completeness AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             comp.Prevalence_completeness_left / 100.0 AS Prevalence_completeness_left,
             comp.Distance_contaminated_left,
             comp.Min_missing_left,
             comp.Prevalence_completeness_right / 100.0 AS Prevalence_completeness_right,
             comp.Distance_contaminated_right,
             comp.Min_missing_right,
             comp.Total_mismatches,
             comp.Total_deletions,
             comp.Total_insertions,
             comp.Total_reads_clipped,
             comp.Total_reference_clipped,
             COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_insertions, 0) + COALESCE(comp.Total_reads_clipped, 0) AS Score_completeness,
             CASE WHEN c.Contig_length > 0 THEN 100.0 - ((COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_insertions, 0) + COALESCE(comp.Total_reads_clipped, 0)) * 100.0 / c.Contig_length) ELSE NULL END AS Percentage_completeness,
             COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_deletions, 0) + COALESCE(comp.Total_reference_clipped, 0) AS Score_contamination,
             CASE WHEN c.Contig_length > 0 THEN (COALESCE(comp.Total_mismatches, 0) + COALESCE(comp.Total_deletions, 0) + COALESCE(comp.Total_reference_clipped, 0)) * 100.0 / c.Contig_length ELSE NULL END AS Percentage_contamination
         FROM Completeness comp
         JOIN Contig c ON comp.Contig_id = c.Contig_id
         JOIN Sample s ON comp.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_completeness VIEW")?;

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

// ============================================================================
// Legacy API - for backward compatibility during migration
// ============================================================================

/// Create the main database with schema and initial data.
/// This is the legacy API - prefer using DbWriter::create() for new code.
pub fn create_metadata_db(
    db_path: &Path,
    contigs: &[ContigInfo],
    annotations: &[FeatureAnnotation],
    _create_indexes: bool, // Ignored - DuckDB doesn't need explicit indexes
) -> Result<()> {
    // Just create an empty database with schema
    let writer = DbWriter::create(db_path, contigs, annotations)?;
    writer.finalize()?;
    Ok(())
}
