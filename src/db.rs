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
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;

use crate::gc_content::{GCContentRun, GCSkewRun, GCSkewStats, GCStats};
use crate::types::{feature_table_name, ContigInfo, FeatureAnnotation, FeaturePoint, PackagingData, PresenceData, VARIABLES};

/// Thread-safe database connection wrapper for sequential writes.
pub struct DbWriter {
    conn: Mutex<Connection>,
    contig_name_to_id: HashMap<String, i64>,
    sample_name_to_id: Mutex<HashMap<String, i64>>,
    next_sample_id: Mutex<i64>,
    variable_name_to_id: HashMap<String, i64>,
    /// Tracks which feature tables have been created (lazy creation)
    created_tables: Mutex<HashSet<String>>,
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

        // Populate Annotated_types table with distinct Type values from Contig_annotation
        populate_annotated_types(&conn)?;

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
            created_tables: Mutex::new(HashSet::new()),
        })
    }

    /// Insert a sample and return its ID.
    pub fn insert_sample(&self, sample_name: &str, sequencing_type: &str, total_reads: u64, mapped_reads: u64) -> Result<i64> {
        // Get next sample ID
        let sample_id = {
            let mut next_id = self.next_sample_id.lock().unwrap();
            let id = *next_id;
            *next_id += 1;
            id
        };

        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO Sample (Sample_id, Sample_name, Sequencing_type, Number_of_reads, Number_of_mapped_reads) VALUES (?1, ?2, ?3, ?4, ?5)",
            params![sample_id, sample_name, sequencing_type, total_reads as i64, mapped_reads as i64],
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
        circular: bool,
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
        self.write_packaging(&conn, sample_id, packaging, circular)?;

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
                // Store coverage_percentage, coverage_variation, and coverage_sd as integers (scaled)
                // Store coverage_mean and coverage_median as integers (rounded)
                let cov_pct = p.coverage_pct.round() as i32;
                let cov_var = p.coverage_variation.round() as i32;
                let cov_sd = p.coverage_sd.round() as i32;
                let cov_mean = p.coverage_mean.round() as i32;
                let cov_median = p.coverage_median.round() as i32;
                appender.append_row(params![contig_id, sample_id, cov_pct, cov_var, cov_sd, cov_mean, cov_median])?;
            }
        }

        appender.flush().context("Failed to flush Presences appender")?;
        Ok(())
    }

    fn write_packaging(&self, conn: &Connection, sample_id: i64, packaging: &[PackagingData], circular: bool) -> Result<()> {
        // Track next packaging_id and terminus_id
        // Query current max to ensure unique IDs
        let mut packaging_id: i64 = conn.query_row(
            "SELECT COALESCE(MAX(Packaging_id), 0) FROM PhageMechanisms",
            [],
            |row| row.get(0),
        ).unwrap_or(0) + 1;

        let mut terminus_id: i64 = conn.query_row(
            "SELECT COALESCE(MAX(Terminus_id), 0) FROM PhageTermini",
            [],
            |row| row.get(0),
        ).unwrap_or(0) + 1;

        let mut mechanism_appender = conn.appender("PhageMechanisms")
            .context("Failed to create PhageMechanisms appender")?;
        let mut termini_appender = conn.appender("PhageTermini")
            .context("Failed to create PhageTermini appender")?;

        for pkg in packaging {
            if let Some(&contig_id) = self.contig_name_to_id.get(&pkg.contig_name) {
                // Collect kept terminus center positions for terminase distance calculation
                let kept_terminus_positions: Vec<i32> = pkg.left_termini.iter()
                    .chain(pkg.right_termini.iter())
                    .filter(|t| t.kept)
                    .map(|t| t.center_pos)
                    .collect();

                // Get genome length for circular distance calculation
                let genome_length: i32 = conn.query_row(
                    "SELECT Contig_length FROM Contig WHERE Contig_id = ?",
                    params![contig_id],
                    |row| row.get(0),
                ).unwrap_or(0);

                // Calculate terminase_distance: minimal distance from any kept terminus
                // center to any terminase gene annotation
                let terminase_distance = calculate_terminase_distance(
                    conn, contig_id, &kept_terminus_positions, genome_length, circular,
                );

                // Insert into PhageMechanisms
                mechanism_appender.append_row(params![
                    packaging_id,
                    contig_id,
                    sample_id,
                    &pkg.mechanism,
                    pkg.duplication,
                    pkg.repeat_length,
                    terminase_distance
                ])?;

                // Insert left termini (status = "start")
                for terminus in &pkg.left_termini {
                    // Store tau as integer (×100)
                    let tau_int = (terminus.tau * 100.0).round() as i32;
                    let kept_str = if terminus.kept { "yes" } else { "no" };
                    termini_appender.append_row(params![
                        terminus_id,
                        packaging_id,
                        terminus.start_pos,
                        terminus.end_pos,
                        terminus.center_pos,
                        "start",
                        terminus.total_spc as i32,
                        terminus.coverage as i32,
                        tau_int,
                        terminus.number_peaks as i32,
                        kept_str,
                        terminus.sum_clippings as i64,
                        (terminus.clipped_ratio * 100.0).round() as i32,
                        terminus.expected_clippings.round() as i64
                    ])?;
                    terminus_id += 1;
                }

                // Insert right termini (status = "end")
                for terminus in &pkg.right_termini {
                    let tau_int = (terminus.tau * 100.0).round() as i32;
                    let kept_str = if terminus.kept { "yes" } else { "no" };
                    termini_appender.append_row(params![
                        terminus_id,
                        packaging_id,
                        terminus.start_pos,
                        terminus.end_pos,
                        terminus.center_pos,
                        "end",
                        terminus.total_spc as i32,
                        terminus.coverage as i32,
                        tau_int,
                        terminus.number_peaks as i32,
                        kept_str,
                        terminus.sum_clippings as i64,
                        (terminus.clipped_ratio * 100.0).round() as i32,
                        terminus.expected_clippings.round() as i64
                    ])?;
                    terminus_id += 1;
                }

                packaging_id += 1;
            }
        }

        mechanism_appender.flush().context("Failed to flush PhageMechanisms appender")?;
        termini_appender.flush().context("Failed to flush PhageTermini appender")?;
        Ok(())
    }

    fn write_completeness(&self, conn: &Connection, sample_id: i64, completeness: &[CompletenessData]) -> Result<()> {
        let mut appender = conn.appender("Completeness")
            .context("Failed to create Completeness appender")?;

        for data in completeness {
            if data.has_data() {
                if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                    // Store prevalences as integers (already percentages)
                    let prev_left = data.prevalence_left.map(|v| v.round() as i32);
                    let prev_right = data.prevalence_right.map(|v| v.round() as i32);
                    // Store totals as integers (already counts or small values)
                    let total_mm = data.total_mismatches.map(|v| v.round() as i32);
                    let total_del = data.total_deletions.map(|v| v.round() as i32);
                    let total_ins = data.total_insertions.map(|v| v.round() as i32);
                    let total_rc = data.total_reads_clipped.map(|v| v.round() as i32);
                    let total_ref = data.total_reference_clipped.map(|v| v.round() as i32);
                    // Circularising reads
                    let circ_reads = data.circularising_reads.map(|v| v as i64);
                    let circ_pct = data.circularising_reads_percentage;

                    appender.append_row(params![
                        contig_id, sample_id,
                        prev_left, data.left_contamination_length, data.left_missing_length,
                        prev_right, data.right_contamination_length, data.right_missing_length,
                        total_mm, total_del, total_ins, total_rc, total_ref,
                        circ_reads, circ_pct
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

        // Lock created_tables for lazy table creation
        let mut created_tables = self.created_tables.lock().unwrap();

        for v in VARIABLES {
            // Skip primary_reads - it's a VIEW
            // Skip contig-level variables stored in Contig_* tables (not per-sample)
            if v.name == "primary_reads" || v.name == "direct_repeats" || v.name == "inverted_repeats" || v.name == "gc_content" {
                continue;
            }

            let Some(feature_points) = by_variable.get(v.name) else {
                continue;
            };

            let table_name = feature_table_name(v.name);
            let has_stats = features_with_stats.contains(&v.name);

            // Create table lazily if it doesn't exist yet
            create_feature_table_if_needed(conn, &table_name, has_stats, &mut created_tables)?;

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

    /// Write repeats data to the database (both direct and inverted).
    /// This is called once during database creation (not per-sample).
    pub fn write_repeats(&self, repeats: &[RepeatsData]) -> Result<()> {
        if repeats.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().unwrap();
        let mut direct_appender = conn.appender("Contig_directRepeats")
            .context("Failed to create Contig_directRepeats appender")?;
        let mut inverted_appender = conn.appender("Contig_invertedRepeats")
            .context("Failed to create Contig_invertedRepeats appender")?;

        let mut direct_count = 0;
        let mut inverted_count = 0;

        for dup in repeats {
            let contig_id_opt = self.contig_name_to_id.get(&dup.contig_name)
                .or_else(|| {
                    // Try without version suffix (e.g., "NC_003071.7" -> "NC_003071")
                    dup.contig_name.rfind('.')
                        .filter(|&dot_pos| dup.contig_name[dot_pos+1..].chars().all(|c| c.is_ascii_digit()))
                        .and_then(|dot_pos| self.contig_name_to_id.get(&dup.contig_name[..dot_pos]))
                });
            if let Some(&contig_id) = contig_id_opt {
                // Store pident as INTEGER (×100)
                let pident_int = (dup.pident * 100.0).round() as i32;

                if dup.is_direct {
                    direct_appender.append_row(params![
                        contig_id,
                        dup.position1,
                        dup.position2,
                        dup.position1prime,
                        dup.position2prime,
                        pident_int
                    ])?;
                    direct_count += 1;
                } else {
                    inverted_appender.append_row(params![
                        contig_id,
                        dup.position1,
                        dup.position2,
                        dup.position1prime,
                        dup.position2prime,
                        pident_int
                    ])?;
                    inverted_count += 1;
                }
            }
        }

        direct_appender.flush().context("Failed to flush Contig_directRepeats appender")?;
        inverted_appender.flush().context("Failed to flush Contig_invertedRepeats appender")?;
        eprintln!("Wrote {} direct repeats and {} inverted repeats to database", direct_count, inverted_count);

        // Warn about contig names from autoblast that couldn't be matched
        let unmatched: std::collections::HashSet<&str> = repeats.iter()
            .filter(|dup| {
                self.contig_name_to_id.get(&dup.contig_name).is_none()
                    && dup.contig_name.rfind('.')
                        .filter(|&p| dup.contig_name[p+1..].chars().all(|c| c.is_ascii_digit()))
                        .and_then(|p| self.contig_name_to_id.get(&dup.contig_name[..p]))
                        .is_none()
            })
            .map(|dup| dup.contig_name.as_str())
            .collect();
        if !unmatched.is_empty() {
            eprintln!("WARNING: {} contig(s) from autoblast not found in database (check naming):", unmatched.len());
            for name in unmatched.iter().take(5) {
                eprintln!("  - '{}'", name);
            }
        }

        // Calculate and update duplication percentages for all contigs
        drop(direct_appender);
        drop(inverted_appender);
        self.update_duplication_percentages(&conn, repeats)?;

        Ok(())
    }

    /// Calculate and update Duplication_percentage for all contigs based on repeat data.
    /// For each contig, merges overlapping repeat intervals (both copies) and calculates
    /// the percentage of the contig covered by repeats.
    fn update_duplication_percentages(&self, conn: &Connection, repeats: &[RepeatsData]) -> Result<()> {
        // Build a map of contig_id -> contig_length from the database
        let contig_lengths: HashMap<i64, i64> = {
            let mut stmt = conn.prepare("SELECT Contig_id, Contig_length FROM Contig")?;
            let rows = stmt.query_map([], |row| Ok((row.get::<_, i64>(0)?, row.get::<_, i64>(1)?)))?;
            rows.filter_map(|r| r.ok()).collect()
        };

        // Group repeat intervals by contig
        let mut intervals_by_contig: HashMap<i64, Vec<(i32, i32)>> = HashMap::new();
        for rep in repeats {
            if let Some(&contig_id) = self.contig_name_to_id.get(&rep.contig_name) {
                // Add both copies of the repeat (first copy and second copy)
                // Normalize so start <= end
                let (start1, end1) = (rep.position1.min(rep.position2), rep.position1.max(rep.position2));
                let (start2, end2) = (rep.position1prime.min(rep.position2prime), rep.position1prime.max(rep.position2prime));

                intervals_by_contig.entry(contig_id).or_default().push((start1, end1));
                intervals_by_contig.entry(contig_id).or_default().push((start2, end2));
            }
        }

        // Set Duplication_percentage to 0 for all contigs without repeats
        let contig_ids_with_repeats: HashSet<i64> = intervals_by_contig.keys().copied().collect();
        for &contig_id in contig_lengths.keys() {
            if !contig_ids_with_repeats.contains(&contig_id) {
                conn.execute(
                    "UPDATE Contig SET Duplication_percentage = 0 WHERE Contig_id = ?1",
                    params![contig_id],
                )?;
            }
        }

        // Calculate and update duplication percentage for each contig with repeats
        for (contig_id, intervals) in intervals_by_contig {
            let contig_length = *contig_lengths
                .get(&contig_id)
                .expect("contig_id not found");

            // Merge overlapping intervals
            let merged = merge_intervals(intervals);

            // Sum the total base pairs covered by merged intervals
            let total_bp: i64 = merged.iter().map(|(start, end)| (*end - *start + 1) as i64).sum();

            // Calculate percentage (rounded to integer)
            let percentage = ((total_bp as f64 / contig_length as f64) * 100.0).round() as i32;

            // Update the Contig table
            conn.execute(
                "UPDATE Contig SET Duplication_percentage = ?1 WHERE Contig_id = ?2",
                params![percentage, contig_id],
            )?;
        }

        Ok(())
    }

    /// Write GC content data to the database.
    /// This is called once during database creation (not per-sample).
    pub fn write_gc_content(&self, gc_data: &[GCContentData]) -> Result<()> {
        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().unwrap();
        let mut appender = conn.appender("Contig_GCContent")
            .context("Failed to create Contig_GCContent appender")?;

        let mut total_runs = 0;

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                for run in &data.runs {
                    appender.append_row(params![
                        contig_id,
                        run.start_pos,
                        run.end_pos,
                        run.gc_percentage as i32
                    ])?;
                    total_runs += 1;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_GCContent appender")?;
        eprintln!("Wrote {} GC content runs to database", total_runs);

        Ok(())
    }

    /// Write GC skew data to the database.
    /// This is called once during database creation (not per-sample).
    pub fn write_gc_skew(&self, gc_data: &[GCContentData]) -> Result<()> {
        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().unwrap();
        let mut appender = conn.appender("Contig_GCSkew")
            .context("Failed to create Contig_GCSkew appender")?;

        let mut total_runs = 0;

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                for run in &data.skew_runs {
                    appender.append_row(params![
                        contig_id,
                        run.start_pos,
                        run.end_pos,
                        run.gc_skew as i32
                    ])?;
                    total_runs += 1;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_GCSkew appender")?;
        eprintln!("Wrote {} GC skew runs to database", total_runs);

        Ok(())
    }

    /// Update Contig table with GC statistics (average, sd) and GC skew stats.
    /// Called after computing GC content and GC skew for all contigs.
    pub fn update_contig_gc_stats(&self, gc_data: &[GCContentData]) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "UPDATE Contig SET GC_mean = ?, GC_sd = ?, GC_skew_amplitude = ?, Positive_GC_skew_windows_percentage = ? WHERE Contig_name = ?"
        )?;

        for data in gc_data {
            stmt.execute(params![
                data.stats.average.round() as i32,                    // GC_mean as int (0-100)
                (data.stats.sd * 100.0).round() as i32,               // GC_sd * 100
                (data.skew_stats.amplitude * 100.0).round() as i32,   // GC_skew_amplitude * 100
                data.skew_stats.percent_positive.round() as i32,      // Positive_GC_skew_windows_percentage as int (0-100)
                &data.contig_name
            ])?;
        }

        Ok(())
    }

    /// Finalize the database after all samples are processed.
    pub fn finalize(self) -> Result<()> {
        let conn = self.conn.into_inner().unwrap();
        let created_tables = self.created_tables.into_inner().unwrap();

        // Only create views if underlying tables exist
        create_views(&conn, &created_tables)?;

        // Delete Variable entries for features without tables
        cleanup_unused_variables(&conn, &created_tables)?;

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

/// Calculate the minimal distance between any kept terminus center position
/// and any terminase gene annotation for a given contig.
/// Returns None if no terminase genes are found or no terminus positions are provided.
/// When circular=true, uses circular distance (min of linear or wrap-around).
fn calculate_terminase_distance(
    conn: &Connection,
    contig_id: i64,
    terminus_positions: &[i32],
    genome_length: i32,
    circular: bool,
) -> Option<i32> {
    if terminus_positions.is_empty() {
        return None;
    }

    // Query terminase gene annotations for this contig
    let mut stmt = conn.prepare(
        "SELECT \"Start\", \"End\" FROM Contig_annotation WHERE Contig_id = ? AND Product LIKE '%terminase%'"
    ).ok()?;

    let terminase_genes: Vec<(i32, i32)> = stmt
        .query_map(params![contig_id], |row| {
            Ok((row.get::<_, i32>(0)?, row.get::<_, i32>(1)?))
        })
        .ok()?
        .filter_map(|r| r.ok())
        .collect();

    if terminase_genes.is_empty() {
        return None;
    }

    // Calculate minimum distance from any terminus to any terminase gene edge
    let mut min_dist = i32::MAX;
    for &term_pos in terminus_positions {
        for &(gene_start, gene_end) in &terminase_genes {
            // Distance to nearest edge of the gene (linear)
            let linear_dist = if term_pos < gene_start {
                gene_start - term_pos
            } else if term_pos > gene_end {
                term_pos - gene_end
            } else {
                0 // terminus is within the gene
            };
            // For circular genomes, consider wrap-around distance
            let dist = if circular && linear_dist > 0 {
                linear_dist.min(genome_length - linear_dist)
            } else {
                linear_dist
            };
            min_dist = min_dist.min(dist);
        }
    }

    Some(min_dist)
}

/// Create core database tables (no indexes - DuckDB uses zone maps).
fn create_core_tables(conn: &Connection) -> Result<()> {
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Duplication_percentage INTEGER,
            GC_mean INTEGER,
            GC_sd INTEGER,
            GC_skew_amplitude INTEGER,
            Positive_GC_skew_windows_percentage INTEGER
        )",
        [],
    )
    .context("Failed to create Contig table")?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY,
            Sample_name TEXT UNIQUE,
            Sequencing_type TEXT,
            Number_of_reads INTEGER,
            Number_of_mapped_reads INTEGER
        )",
        [],
    )
    .context("Failed to create Sample table")?;

    // Aligned_fraction_percentage, Coverage_variation, and Coverage_sd stored as INTEGER (scaled)
    // Coverage_mean and Coverage_median stored as INTEGER (rounded)
    // No Presence_id needed - (Contig_id, Sample_id) is the natural key
    conn.execute(
        "CREATE TABLE Presences (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Aligned_fraction_percentage INTEGER,
            Coverage_variation INTEGER,
            Coverage_sd INTEGER,
            Coverage_mean INTEGER,
            Coverage_median INTEGER,
            PRIMARY KEY (Contig_id, Sample_id)
        )",
        [],
    )
    .context("Failed to create Presences table")?;

    conn.execute(
        "CREATE TABLE PhageMechanisms (
            Packaging_id INTEGER PRIMARY KEY,
            Contig_id INTEGER NOT NULL,
            Sample_id INTEGER NOT NULL,
            Packaging_mechanism TEXT NOT NULL,
            Duplication BOOLEAN,
            Repeat_length INTEGER,
            Terminase_distance INTEGER
        )",
        [],
    )
    .context("Failed to create PhageMechanisms table")?;

    // PhageTermini table - stores individual terminus areas with full metadata
    // One row per terminus area, linked to PhageMechanisms via Packaging_id
    // Includes filtering diagnostics: NumberPeaks, Kept, Clippings, Clipping_excess, Max_clippings_allowed
    conn.execute(
        "CREATE TABLE PhageTermini (
            Terminus_id INTEGER PRIMARY KEY,
            Packaging_id INTEGER NOT NULL,
            \"Start\" INTEGER NOT NULL,
            \"End\" INTEGER NOT NULL,
            Center INTEGER NOT NULL,
            Status TEXT NOT NULL,
            SPC INTEGER NOT NULL,
            Coverage INTEGER NOT NULL,
            Tau INTEGER NOT NULL,
            NumberPeaks INTEGER NOT NULL,
            Kept TEXT NOT NULL,
            Clippings INTEGER NOT NULL,
            Clipping_excess INTEGER NOT NULL,
            Max_clippings_allowed INTEGER NOT NULL
        )",
        [],
    )
    .context("Failed to create PhageTermini table")?;

    // Contig_directRepeats table - stores direct repeats (same orientation)
    // Pident stored as INTEGER (×100)
    conn.execute(
        "CREATE TABLE Contig_directRepeats (
            Contig_id INTEGER,
            Position1 INTEGER,
            Position2 INTEGER,
            Position1prime INTEGER,
            Position2prime INTEGER,
            Pident INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_directRepeats table")?;

    // Contig_invertedRepeats table - stores inverted repeats (opposite orientation)
    // Pident stored as INTEGER (×100)
    conn.execute(
        "CREATE TABLE Contig_invertedRepeats (
            Contig_id INTEGER,
            Position1 INTEGER,
            Position2 INTEGER,
            Position1prime INTEGER,
            Position2prime INTEGER,
            Pident INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_invertedRepeats table")?;

    // Contig_GCContent table - stores GC content (contig-level, sample-independent)
    // Value stored as INTEGER (0-100 percentage)
    conn.execute(
        "CREATE TABLE Contig_GCContent (
            Contig_id INTEGER,
            First_position INTEGER,
            Last_position INTEGER,
            Value INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_GCContent table")?;

    // Contig_GCSkew table - stores GC skew (contig-level, sample-independent)
    // Value stored as INTEGER × 100 (range: -100 to +100)
    conn.execute(
        "CREATE TABLE Contig_GCSkew (
            Contig_id INTEGER,
            First_position INTEGER,
            Last_position INTEGER,
            Value INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_GCSkew table")?;

    // Completeness table - clipping prevalences stored as INTEGER (×100), totals as INTEGER
    // Note: Left/Right_clippings_percentage stores clipping_count/coverage (raw clipping prevalence)
    // The Explicit_completeness VIEW derives completeness as: 1 - clipping_prevalence
    conn.execute(
        "CREATE TABLE Completeness (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Left_clippings_percentage INTEGER,
            Left_contamination_length INTEGER,
            Left_missing_length INTEGER,
            Right_clippings_percentage INTEGER,
            Right_contamination_length INTEGER,
            Right_missing_length INTEGER,
            Mismatch_frequency INTEGER,
            Deletion_frequency INTEGER,
            Insertion_frequency INTEGER,
            Read_based_clipping_frequency INTEGER,
            Reference_based_clippings_frequency INTEGER,
            Circularising_reads INTEGER,
            Circularising_reads_percentage INTEGER
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
            Phrog INTEGER,
            Locus_tag TEXT,
            Longest_isoform BOOLEAN
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
            Module_order INTEGER,
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

    // Annotated_types table - stores distinct Type values from Contig_annotation ordered by frequency
    conn.execute(
        "CREATE TABLE Annotated_types (
            Type_id INTEGER PRIMARY KEY,
            Type_name TEXT UNIQUE,
            Frequency INTEGER
        )",
        [],
    )
    .context("Failed to create Annotated_types table")?;

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

    // Constants table - stores metadata flags about database content
    // Pre-computed at database creation time (not per-request)
    conn.execute(
        "CREATE TABLE Constants (
            Constant TEXT PRIMARY KEY,
            Status BOOLEAN
        )",
        [],
    )
    .context("Failed to create Constants table")?;

    Ok(())
}

/// Insert contigs into the database using Appender for bulk loading.
fn insert_contigs(conn: &Connection, contigs: &[ContigInfo]) -> Result<()> {
    let mut appender = conn.appender("Contig")
        .context("Failed to create Contig appender")?;

    for (i, contig) in contigs.iter().enumerate() {
        let null_int: Option<i32> = None;
        appender.append_row(params![
            (i + 1) as i64,
            &contig.name,
            contig.length as i64,
            null_int,  // Duplication_percentage - set later from autoblast
            null_int,  // GC_mean - set later from GC content computation (stored as int 0-100)
            null_int,  // GC_sd - set later from GC content computation (stored as int, ×100)
            null_int,  // GC_skew_amplitude - set later from GC skew computation (stored as int, ×100)
            null_int,  // Positive_GC_skew_windows_percentage - set later from GC skew computation (stored as int 0-100)
        ])
        .with_context(|| format!("Failed to append contig: {}", contig.name))?;
    }

    appender.flush().context("Failed to flush Contig appender")?;
    Ok(())
}

/// Insert annotations into the database using Appender for bulk loading.
fn insert_annotations(conn: &Connection, annotations: &[FeatureAnnotation]) -> Result<()> {
    // PHAROKKA function categories (must match plotting_data_per_sample.py)
    let pharokka_keys = vec![
        "vfdb_card",
        "unknown function",
        "other",
        "tail",
        "transcription regulation",
        "dna, rna and nucleotide metabolism",
        "lysis",
        "moron, auxiliary metabolic gene and host takeover",
        "integration and excision",
        "head and packaging",
        "connector",
    ];

    // First pass: group annotations by (locus_tag, feature_type) to find longest per group
    // Only group features that HAVE a locus_tag (features without locus_tag are always shown)
    let mut group_max_length: std::collections::HashMap<(String, String), usize> = std::collections::HashMap::new();
    
    // Also check if any CDS feature Function contains PHAROKKA keys
    let mut has_pharokka_function = false;
    
    for ann in annotations {
        let length = (ann.end - ann.start + 1) as usize;
        
        // Only group annotations that have a locus_tag
        if let Some(ref locus_tag) = ann.locus_tag {
            let key = (locus_tag.clone(), ann.feature_type.clone());
            group_max_length.entry(key)
                .and_modify(|max| { if length > *max { *max = length; } })
                .or_insert(length);
        }
        
        // Check if this is a CDS with a PHAROKKA function
        if ann.feature_type == "CDS" {
            if let Some(function) = &ann.function {
                let func_lower = function.to_lowercase();
                for pharokka_key in &pharokka_keys {
                    if func_lower.contains(pharokka_key) {
                        has_pharokka_function = true;
                        break;
                    }
                }
                if has_pharokka_function {
                    break;  // Early exit once we find one
                }
            }
        }
    }

    let mut appender = conn.appender("Contig_annotation")
        .context("Failed to create Contig_annotation appender")?;

    // Track which (locus_tag, Type) groups have already had their longest assigned
    // This ensures only the FIRST occurrence gets TRUE when multiple isoforms have equal max length
    let mut longest_assigned: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();

    for ann in annotations {
        let length = (ann.end - ann.start + 1) as usize;
        
        // Determine Longest_isoform value:
        // - NULL if no locus_tag (always displayed)
        // - TRUE if this is the longest in its (locus_tag, Type) group AND is the first to be marked
        // - FALSE if it's not the longest in its (locus_tag, Type) group OR already marked
        let longest_isoform: Option<bool> = if let Some(ref locus_tag) = ann.locus_tag {
            let key = (locus_tag.clone(), ann.feature_type.clone());
            let max_length = *group_max_length.get(&key).unwrap_or(&0);
            
            // Only mark as longest if:
            // 1. Length matches the max for this group
            // 2. This group hasn't been assigned a longest yet (tie-breaking: first wins)
            let is_longest = if length == max_length && !longest_assigned.contains(&key) {
                longest_assigned.insert(key);
                true
            } else {
                false
            };
            Some(is_longest)
        } else {
            None  // NULL for features without locus_tag
        };
        
        appender.append_row(params![
            ann.contig_id,
            ann.start,
            ann.end,
            ann.strand,
            &ann.feature_type,
            &ann.product,
            &ann.function,
            &ann.phrog,
            &ann.locus_tag,
            longest_isoform,
        ])
        .with_context(|| format!(
            "Failed to append annotation: contig_id={}, start={}, end={}, type={}",
            ann.contig_id, ann.start, ann.end, ann.feature_type
        ))?;
    }

    appender.flush().context("Failed to flush Contig_annotation appender")?;
    
    // Check if any locus_tag appears more than once (isoforms present)
    let mut locus_tag_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for ann in annotations {
        if let Some(ref locus_tag) = ann.locus_tag {
            *locus_tag_counts.entry(locus_tag.clone()).or_insert(0) += 1;
        }
    }
    let has_isoforms = locus_tag_counts.values().any(|&count| count > 1);
    
    // Insert constants into Constants table
    conn.execute(
        "INSERT INTO Constants (Constant, Status) VALUES ('pharokka', ?), ('isoforms', ?)",
        params![has_pharokka_function, has_isoforms],
    )
    .context("Failed to insert constants")?;
    
    Ok(())
}

/// Populate Annotated_types table with distinct Type values from Contig_annotation ordered by frequency.
fn populate_annotated_types(conn: &Connection) -> Result<()> {
    conn.execute(
        "INSERT INTO Annotated_types (Type_id, Type_name, Frequency)
         SELECT
             ROW_NUMBER() OVER (ORDER BY COUNT(*) DESC) AS Type_id,
             \"Type\" AS Type_name,
             COUNT(*) AS Frequency
         FROM Contig_annotation
         GROUP BY \"Type\"
         ORDER BY COUNT(*) DESC",
        [],
    )
    .context("Failed to populate Annotated_types table")?;
    Ok(())
}

/// Create Variable metadata entries only. Feature tables are created lazily.
fn create_variable_tables(conn: &Connection) -> Result<()> {
    for (i, v) in VARIABLES.iter().enumerate() {
        // Special case: contig-level data is stored in separate Contig_* tables
        let table_name = match v.name {
            "direct_repeats" => "Contig_directRepeats".to_string(),
            "inverted_repeats" => "Contig_invertedRepeats".to_string(),
            "gc_content" => "Contig_GCContent".to_string(),
            "gc_skew" => "Contig_GCSkew".to_string(),
            _ => feature_table_name(v.name),
        };

        conn.execute(
            "INSERT INTO Variable (Variable_id, Variable_name, Subplot, Module, Module_order, \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)",
            params![
                (i + 1) as i64,
                v.name,
                v.subplot,
                v.module,
                v.module_order,
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
    }

    Ok(())
}

/// Create a feature table if it doesn't exist yet.
fn create_feature_table_if_needed(
    conn: &Connection,
    table_name: &str,
    has_stats: bool,
    created_tables: &mut HashSet<String>,
) -> Result<()> {
    if created_tables.contains(table_name) {
        return Ok(());
    }

    let table_sql = if has_stats {
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

    created_tables.insert(table_name.to_string());
    Ok(())
}

/// Create views after all data is inserted.
fn create_views(conn: &Connection, created_tables: &HashSet<String>) -> Result<()> {
    // Only create primary_reads VIEW if both underlying tables exist
    let plus_table = "Feature_primary_reads_plus_only";
    let minus_table = "Feature_primary_reads_minus_only";
    if created_tables.contains(plus_table) && created_tables.contains(minus_table) {
        // Create a VIEW that computes primary_reads as sum of plus and minus strands
        // NOTE: The Bokeh server now bypasses this VIEW by computing in Python
        // (merging Feature_primary_reads_plus_only + Feature_primary_reads_minus_only)
        // to avoid OOM issues. This VIEW is kept for backward compatibility.
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
    }

    // Explicit_presences VIEW with correction ratios for normalizing across samples
    conn.execute(
        "CREATE VIEW Explicit_presences AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             p.Aligned_fraction_percentage,
             p.Coverage_variation / 1000000.0 AS Coverage_variation,
             p.Coverage_sd / 1000000.0 AS Coverage_sd,
             p.Coverage_mean,
             p.Coverage_median,
             CAST(s.Number_of_reads AS REAL) / (SELECT MIN(Number_of_reads) FROM Sample WHERE Number_of_reads > 0) AS Read_number_correction_ratio,
             CAST(s.Number_of_mapped_reads AS REAL) / (SELECT MIN(Number_of_mapped_reads) FROM Sample WHERE Number_of_mapped_reads > 0) AS Read_mapped_correction_ratio,
             p.Coverage_mean / (CAST(s.Number_of_reads AS REAL) / (SELECT MIN(Number_of_reads) FROM Sample WHERE Number_of_reads > 0)) AS Coverage_mean_corrected_by_number_of_reads,
             p.Coverage_mean / (CAST(s.Number_of_mapped_reads AS REAL) / (SELECT MIN(Number_of_mapped_reads) FROM Sample WHERE Number_of_mapped_reads > 0)) AS Coverage_mean_corrected_by_number_of_mapped_reads,
             p.Coverage_median / (CAST(s.Number_of_reads AS REAL) / (SELECT MIN(Number_of_reads) FROM Sample WHERE Number_of_reads > 0)) AS Coverage_median_corrected_by_number_of_reads,
             p.Coverage_median / (CAST(s.Number_of_mapped_reads AS REAL) / (SELECT MIN(Number_of_mapped_reads) FROM Sample WHERE Number_of_mapped_reads > 0)) AS Coverage_median_corrected_by_number_of_mapped_reads
         FROM Presences p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_presences VIEW")?;

    // Explicit_completeness VIEW - derives completeness from clipping prevalence
    // Completeness = 1 - clipping_prevalence (higher clipping = lower completeness)
    conn.execute(
        "CREATE VIEW Explicit_completeness AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             COALESCE(100 - comp.Left_clippings_percentage, 100) AS Left_completeness_percentage,
             COALESCE(comp.Left_contamination_length, 0) AS Left_contamination_length,
             COALESCE(comp.Left_missing_length, 0) AS Left_missing_length,
             COALESCE(100 - comp.Right_clippings_percentage, 100) AS Right_completeness_percentage,
             COALESCE(comp.Right_contamination_length, 0) AS Right_contamination_length,
             COALESCE(comp.Right_missing_length, 0) AS Right_missing_length,
             CASE WHEN c.Contig_length > 0 THEN COALESCE(comp.Mismatch_frequency, 0.0) / c.Contig_length ELSE 0 END AS Mismatch_frequency,
             CASE WHEN c.Contig_length > 0 THEN COALESCE(comp.Deletion_frequency, 0.0) / c.Contig_length ELSE 0 END AS Deletion_frequency,
             CASE WHEN c.Contig_length > 0 THEN COALESCE(comp.Insertion_frequency, 0.0) / c.Contig_length ELSE 0 END AS Insertion_frequency,
             CASE WHEN c.Contig_length > 0 THEN COALESCE(comp.Read_based_clipping_frequency, 0.0) / c.Contig_length ELSE 0 END AS Read_based_clipping_frequency,
             CASE WHEN c.Contig_length > 0 THEN COALESCE(comp.Reference_based_clippings_frequency, 0.0) / c.Contig_length ELSE 0 END AS Reference_based_clippings_frequency,
             COALESCE(comp.Mismatch_frequency, 0) + COALESCE(comp.Insertion_frequency, 0) + COALESCE(comp.Read_based_clipping_frequency, 0) AS Score_completeness,
             CASE WHEN c.Contig_length > 0 THEN GREATEST(0, ROUND(
                100 - (
                    (COALESCE(comp.Mismatch_frequency, 0) + COALESCE(comp.Insertion_frequency, 0) + COALESCE(comp.Read_based_clipping_frequency, 0)) * 100.0 / c.Contig_length
                )
             ))::INTEGER ELSE 100 END AS Completeness_percentage,
             COALESCE(comp.Mismatch_frequency, 0) + COALESCE(comp.Deletion_frequency, 0) + COALESCE(comp.Reference_based_clippings_frequency, 0) AS Score_contamination,
             CASE WHEN c.Contig_length > 0 THEN LEAST(100, ROUND(
                (COALESCE(comp.Mismatch_frequency, 0) + COALESCE(comp.Deletion_frequency, 0) + COALESCE(comp.Reference_based_clippings_frequency, 0)) * 100.0 / c.Contig_length
             ))::INTEGER ELSE 0 END AS Contamination_percentage,
             COALESCE(comp.Circularising_reads, 0) AS Circularising_reads,
             COALESCE(comp.Circularising_reads_percentage, 0) AS Circularising_reads_percentage
         FROM Completeness comp
         JOIN Contig c ON comp.Contig_id = c.Contig_id
         JOIN Sample s ON comp.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_completeness VIEW")?;

    // Explicit_phage_mechanisms VIEW - backwards compatible with comma-separated termini
    // Includes aggregated diagnostics columns from PhageTermini
    conn.execute(
        "CREATE VIEW Explicit_phage_mechanisms AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             m.Packaging_mechanism,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM PhageTermini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'start' AND pt.Kept = 'yes'),
                 ''
             ) AS Left_termini,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM PhageTermini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'end' AND pt.Kept = 'yes'),
                 ''
             ) AS Right_termini,
             CASE WHEN m.Duplication = true THEN 'DTR'
                  WHEN m.Duplication = false THEN 'ITR'
                  ELSE NULL END AS Duplication,
             COALESCE(
                 (SELECT COUNT(*)
                  FROM PhageTermini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Kept = 'yes'),
                 0
             ) AS Total_peaks,
             m.Repeat_length,
             m.Terminase_distance,
             CASE WHEN m.Terminase_distance IS NOT NULL AND c.Contig_length > 0
                  THEN CAST(ROUND(CAST(m.Terminase_distance AS REAL) / c.Contig_length * 100) AS INTEGER)
                  ELSE NULL END AS Terminase_percentage
         FROM Presences p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id
         LEFT JOIN PhageMechanisms m ON p.Contig_id = m.Contig_id AND p.Sample_id = m.Sample_id",
        [],
    )
    .context("Failed to create Explicit_phage_mechanisms VIEW")?;

    // Explicit_phage_termini VIEW
    conn.execute_batch(
        "CREATE VIEW Explicit_phage_termini AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             m.Packaging_mechanism,
             CASE WHEN m.Duplication = true THEN 'DTR'
                  WHEN m.Duplication = false THEN 'ITR'
                  ELSE NULL END AS Duplication,
             m.Repeat_length,
             m.Terminase_distance,
             CASE WHEN m.Terminase_distance IS NOT NULL AND c.Contig_length > 0
                  THEN CAST(ROUND(CAST(m.Terminase_distance AS REAL) / c.Contig_length * 100) AS INTEGER)
                  ELSE NULL END AS Terminase_percentage,
             t.\"Start\",
             t.\"End\",
             t.Center,
             t.Status,
             t.SPC,
             t.Coverage,
             t.Tau,
             t.NumberPeaks,
             t.Kept,
             t.Clippings,
             t.Clipping_excess,
             t.Max_clippings_allowed
         FROM PhageTermini t
         JOIN PhageMechanisms m ON t.Packaging_id = m.Packaging_id
         JOIN Contig c ON m.Contig_id = c.Contig_id
         JOIN Sample s ON m.Sample_id = s.Sample_id",
    )
    .context("Failed to create Explicit_phage_termini VIEW")?;

    Ok(())
}

/// Delete Variable entries for features that don't have data (no table created).
fn cleanup_unused_variables(conn: &Connection, created_tables: &HashSet<String>) -> Result<()> {
    // Get all variable names and their table names
    let mut stmt = conn.prepare("SELECT Variable_name, Feature_table_name FROM Variable")?;
    let rows: Vec<(String, String)> = stmt
        .query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?)))?
        .filter_map(|r| r.ok())
        .collect();

    let mut deleted_count = 0;
    for (var_name, table_name) in rows {
        // Skip special cases:
        // - direct_repeats/inverted_repeats/gc_content/gc_skew: data stored in Contig_* tables (always exist)
        // - primary_reads: is a VIEW, not a table (depends on plus/minus tables)
        if var_name == "direct_repeats" || var_name == "inverted_repeats" || var_name == "gc_content" || var_name == "gc_skew" {
            continue;
        }

        // For primary_reads VIEW, check if the underlying tables exist
        if var_name == "primary_reads" {
            let plus_exists = created_tables.contains("Feature_primary_reads_plus_only");
            let minus_exists = created_tables.contains("Feature_primary_reads_minus_only");
            if !plus_exists || !minus_exists {
                conn.execute("DELETE FROM Variable WHERE Variable_name = ?1", params![var_name])?;
                deleted_count += 1;
            }
            continue;
        }

        // For regular features, check if the table was created
        if !created_tables.contains(&table_name) {
            conn.execute("DELETE FROM Variable WHERE Variable_name = ?1", params![var_name])?;
            deleted_count += 1;
        }
    }

    if deleted_count > 0 {
        eprintln!("Removed {} unused variable entries from database", deleted_count);
    }

    Ok(())
}

/// Merge overlapping intervals into non-overlapping intervals.
/// Takes a vector of (start, end) intervals and returns merged intervals.
/// Assumes start <= end for each interval.
fn merge_intervals(mut intervals: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    // Sort by start position
    intervals.sort_by_key(|&(start, _)| start);

    let mut merged: Vec<(i32, i32)> = Vec::new();
    let mut current = intervals[0];

    for &(start, end) in &intervals[1..] {
        if start <= current.1 + 1 {
            // Overlapping or adjacent: extend the current interval
            current.1 = current.1.max(end);
        } else {
            // No overlap: save current and start a new interval
            merged.push(current);
            current = (start, end);
        }
    }
    merged.push(current);

    merged
}

/// Repeats data from self-BLAST results.
/// Represents a repeated region within a contig.
#[derive(Clone, Debug)]
pub struct RepeatsData {
    pub contig_name: String,
    /// Start position of the first copy (query start)
    pub position1: i32,
    /// End position of the first copy (query end)
    pub position2: i32,
    /// Start position of the second copy (subject start) - can be > position2prime if inverted
    pub position1prime: i32,
    /// End position of the second copy (subject end)
    pub position2prime: i32,
    /// Percentage identity
    pub pident: f64,
    /// True if direct repeat (same orientation), false if inverted repeat
    pub is_direct: bool,
}

/// GC content data for a contig.
/// Contains pre-computed GC content runs with RLE compression and statistics.
#[derive(Clone, Debug)]
pub struct GCContentData {
    pub contig_name: String,
    pub runs: Vec<GCContentRun>,
    pub skew_runs: Vec<GCSkewRun>,
    pub stats: GCStats,
    pub skew_stats: GCSkewStats,
}

/// Completeness data for a contig.
/// Contains assembly completeness statistics computed from clipping events at reference ends.
/// Individual score components are stored; score_completeness and score_contamination are computed in VIEW.
#[derive(Clone, Debug)]
pub struct CompletenessData {
    pub contig_name: String,
    pub prevalence_left: Option<f64>,
    pub left_contamination_length: Option<i32>,
    pub left_missing_length: Option<i32>,
    pub prevalence_right: Option<f64>,
    pub right_contamination_length: Option<i32>,
    pub right_missing_length: Option<i32>,
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
    /// Count of reads supporting genome circularity
    pub circularising_reads: Option<u64>,
    /// Percentage of primary reads that are circularising
    pub circularising_reads_percentage: Option<i32>,
}

impl CompletenessData {
    /// Returns true if there is any completeness data.
    pub fn has_data(&self) -> bool {
        self.prevalence_left.is_some() || self.prevalence_right.is_some()
            || self.total_mismatches.is_some() || self.total_deletions.is_some()
            || self.total_insertions.is_some() || self.total_reads_clipped.is_some()
            || self.total_reference_clipped.is_some() || self.circularising_reads.is_some()
    }
}

/// Parse BLAST outfmt 6 file and return repeats data.
///
/// BLAST outfmt 6 columns:
/// 1. qseqid - query sequence id (contig name)
/// 2. sseqid - subject sequence id (same as query for self-blast)
/// 3. pident - percentage identity
/// 4. length - alignment length
/// 5. mismatch - number of mismatches
/// 6. gapopen - number of gap openings
/// 7. qstart - query start
/// 8. qend - query end
/// 9. sstart - subject start
/// 10. send - subject end
/// 11. evalue - e-value
/// 12. bitscore - bit score
pub fn parse_autoblast_file(path: &Path) -> Result<Vec<RepeatsData>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    if !path.exists() || path.as_os_str().is_empty() {
        return Ok(Vec::new());
    }

    let file = File::open(path).with_context(|| format!("Failed to open autoblast file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut repeats = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.with_context(|| format!("Failed to read line {} of autoblast file", line_num + 1))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            eprintln!("Warning: Skipping line {} with {} fields (expected 12)", line_num + 1, fields.len());
            continue;
        }

        let contig_name = fields[0].to_string();
        let pident: f64 = fields[2].parse().unwrap_or(0.0);
        let qstart: i32 = fields[6].parse().unwrap_or(0);
        let qend: i32 = fields[7].parse().unwrap_or(0);
        let sstart: i32 = fields[8].parse().unwrap_or(0);
        let send: i32 = fields[9].parse().unwrap_or(0);

        // Skip self-hits (exact same region)
        if qstart == sstart && qend == send {
            continue;
        }

        // Determine if direct (same orientation) or inverted (opposite orientation)
        let is_direct = (qstart < qend && sstart < send) || (qstart > qend && sstart > send);

        repeats.push(RepeatsData {
            contig_name,
            position1: qstart,
            position2: qend,
            position1prime: sstart,
            position2prime: send,
            pident,
            is_direct,
        });
    }

    Ok(repeats)
}
