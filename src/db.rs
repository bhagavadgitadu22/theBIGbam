//! DuckDB database operations.
//!
//! Database schema:
//! - `Contig`: Contig metadata (name, length, annotation tool)
//! - `Sample`: Sample names from BAM files
//! - `Coverage`: Coverage metrics per contig per sample
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
// Re-export new metric data structs (defined below in this file)

/// Thread-safe database connection wrapper for sequential writes.
pub struct DbWriter {
    conn: Mutex<Connection>,
    has_bam: bool,
    is_extend: bool,
    contig_name_to_id: HashMap<String, i64>,
    sample_name_to_id: Mutex<HashMap<String, i64>>,
    next_sample_id: Mutex<i64>,
    /// Tracks which feature tables have been created (lazy creation)
    created_tables: Mutex<HashSet<String>>,
}

impl DbWriter {
    /// Create a new database and return a writer for sequential sample insertion.
    pub fn create(
        db_path: &Path,
        contigs: &[ContigInfo],
        annotations: &[FeatureAnnotation],
        has_bam: bool,
    ) -> Result<Self> {
        // Remove existing database if present
        let _ = std::fs::remove_file(db_path);

        let conn = Connection::open(db_path)
            .with_context(|| format!("Failed to create database: {}", db_path.display()))?;

        // Create tables
        create_core_tables(&conn, has_bam)?;
        create_variable_tables(&conn)?;

        // Insert contigs
        insert_contigs(&conn, contigs)?;

        // Insert contig sequences (from GenBank or FASTA)
        insert_contig_sequences(&conn, contigs)?;

        // Insert codon table (standard genetic code)
        insert_codon_table(&conn)?;

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

        Ok(Self {
            conn: Mutex::new(conn),
            has_bam,
            is_extend: false,
            contig_name_to_id,
            sample_name_to_id: Mutex::new(HashMap::new()),
            next_sample_id: Mutex::new(1),
            created_tables: Mutex::new(HashSet::new()),
        })
    }

    /// Open an existing database for extending with new samples/contigs.
    pub fn open(
        db_path: &Path,
        new_contigs: &[ContigInfo],
        new_annotations: &[FeatureAnnotation],
        has_bam: bool,
    ) -> Result<Self> {
        let conn = Connection::open(db_path)
            .with_context(|| format!("Failed to open database: {}", db_path.display()))?;

        // Read existing contig name -> id mapping
        let mut contig_name_to_id: HashMap<String, i64> = HashMap::new();
        let mut max_contig_id: i64 = 0;
        {
            let mut stmt = conn.prepare("SELECT Contig_id, Contig_name FROM Contig")?;
            let rows = stmt.query_map([], |row| {
                Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
            })?;
            for row in rows {
                let (id, name) = row?;
                if id > max_contig_id {
                    max_contig_id = id;
                }
                contig_name_to_id.insert(name, id);
            }
        }

        // Insert new contigs with IDs continuing from max
        if !new_contigs.is_empty() {
            let mut appender = conn.appender("Contig")
                .context("Failed to create Contig appender")?;
            for (i, contig) in new_contigs.iter().enumerate() {
                let contig_id = max_contig_id + 1 + i as i64;
                let null_int: Option<i32> = None;
                appender.append_row(params![
                    contig_id,
                    &contig.name,
                    contig.length as i64,
                    null_int, null_int, null_int, null_int, null_int,
                ])?;
                contig_name_to_id.insert(contig.name.clone(), contig_id);
            }
            appender.flush().context("Failed to flush new Contig appender")?;

            // Insert new contig sequences
            insert_contig_sequences_with_offset(&conn, new_contigs, max_contig_id)?;

            // Insert new annotations (remap contig_ids: parser uses 1-based, we offset by max_contig_id)
            if !new_annotations.is_empty() {
                let remapped: Vec<FeatureAnnotation> = new_annotations.iter().map(|a| {
                    let mut ann = a.clone();
                    ann.contig_id += max_contig_id;
                    ann
                }).collect();
                insert_annotations(&conn, &remapped)?;
                // Repopulate Annotated_types from full Contig_annotation table
                conn.execute("DELETE FROM Annotated_types", [])?;
                populate_annotated_types(&conn)?;
            }

            eprintln!("Extended database with {} new contigs", new_contigs.len());
        }

        // Read existing sample name -> id mapping
        let mut sample_name_to_id: HashMap<String, i64> = HashMap::new();
        let mut max_sample_id: i64 = 0;
        {
            let has_sample_table: bool = conn
                .query_row(
                    "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = 'Sample'",
                    [],
                    |row| row.get::<_, i64>(0),
                )? > 0;
            if has_sample_table {
                let mut stmt = conn.prepare("SELECT Sample_id, Sample_name FROM Sample")?;
                let rows = stmt.query_map([], |row| {
                    Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
                })?;
                for row in rows {
                    let (id, name) = row?;
                    if id > max_sample_id {
                        max_sample_id = id;
                    }
                    sample_name_to_id.insert(name, id);
                }
            }
        }

        // Detect existing feature tables
        let mut created_tables: HashSet<String> = HashSet::new();
        {
            let mut stmt = conn.prepare(
                "SELECT table_name FROM information_schema.tables WHERE table_name LIKE 'Feature_%'"
            )?;
            let rows = stmt.query_map([], |row| row.get::<_, String>(0))?;
            for row in rows {
                created_tables.insert(row?);
            }
        }

        // Drop existing materialized tables and views (will be recreated in finalize)
        for table in &[
            "Contig_direct_repeat_count", "Contig_inverted_repeat_count",
            "Contig_direct_repeat_identity", "Contig_inverted_repeat_identity",
        ] {
            conn.execute(&format!("DROP TABLE IF EXISTS {}", table), [])?;
        }
        for view in &[
            "Explicit_coverage", "Explicit_misassembly", "Explicit_microdiversity",
            "Explicit_side_misassembly", "Explicit_topology",
            "Explicit_phage_mechanisms", "Explicit_phage_termini",
        ] {
            conn.execute(&format!("DROP VIEW IF EXISTS {}", view), [])?;
        }

        Ok(Self {
            conn: Mutex::new(conn),
            has_bam,
            is_extend: true,
            contig_name_to_id,
            sample_name_to_id: Mutex::new(sample_name_to_id),
            next_sample_id: Mutex::new(max_sample_id + 1),
            created_tables: Mutex::new(created_tables),
        })
    }

    /// Insert a sample and return its ID.
    pub fn insert_sample(&self, sample_name: &str, sequencing_type: &str, total_reads: u64, mapped_reads: u64, circular_mapping: bool) -> Result<i64> {
        // Get next sample ID
        let sample_id = {
            let mut next_id = self.next_sample_id.lock().unwrap();
            let id = *next_id;
            *next_id += 1;
            id
        };

        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO Sample (Sample_id, Sample_name, Sequencing_type, Number_of_reads, Number_of_mapped_reads, Circular_mapping) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![sample_id, sample_name, sequencing_type, total_reads as i64, mapped_reads as i64, circular_mapping],
        )?;
        drop(conn);

        let mut sample_map = self.sample_name_to_id.lock().unwrap();
        sample_map.insert(sample_name.to_string(), sample_id);

        Ok(sample_id)
    }

    /// Write all data for a sample (coverage, packaging, misassembly/microdiversity/side_misassembly/topology, features).
    pub fn write_sample_data(
        &self,
        sample_name: &str,
        presences: &[PresenceData],
        packaging: &[PackagingData],
        misassembly: &[MisassemblyData],
        microdiversity: &[MicrodiversityData],
        side_misassembly: &[SideMisassemblyData],
        topology: &[TopologyData],
        features: &[FeaturePoint],
        circular: bool,
    ) -> Result<()> {
        let sample_id = {
            let sample_map = self.sample_name_to_id.lock().unwrap();
            *sample_map.get(sample_name).context("Sample not found")?
        };

        let conn = self.conn.lock().unwrap();

        // Appenders handle their own transactions - no explicit BEGIN/COMMIT needed
        // Insert coverage data
        self.write_presences(&conn, sample_id, presences)?;

        // Insert packaging mechanisms
        self.write_packaging(&conn, sample_id, packaging, circular)?;

        // Insert misassembly data
        self.write_misassembly(&conn, sample_id, misassembly)?;

        // Insert microdiversity data
        self.write_microdiversity(&conn, sample_id, microdiversity)?;

        // Insert side misassembly data
        self.write_side_misassembly(&conn, sample_id, side_misassembly)?;

        // Insert topology data
        self.write_topology(&conn, sample_id, topology)?;

        // Insert features
        self.write_features(&conn, sample_id, features)?;

        Ok(())
    }

    fn write_presences(&self, conn: &Connection, sample_id: i64, presences: &[PresenceData]) -> Result<()> {
        let mut appender = conn.appender("Coverage")
            .context("Failed to create Coverage appender")?;

        for p in presences {
            if let Some(&contig_id) = self.contig_name_to_id.get(&p.contig_name) {
                let cov_pct = (p.coverage_pct * 10.0).round() as i32;
                let above_expected = p.above_expected_aligned_fraction;
                let read_count = p.read_count as i64;
                let cov_mean = (p.coverage_mean * 10.0).round() as i32;
                let cov_median = (p.coverage_median * 10.0).round() as i32;
                let cov_trimmed_mean = (p.coverage_trimmed_mean * 10.0).round() as i32;
                let cov_sd = p.coverage_sd.round() as i32;
                let cov_var = p.coverage_variation.round() as i32;
                appender.append_row(params![contig_id, sample_id, cov_pct, above_expected, read_count, cov_mean, cov_median, cov_trimmed_mean, cov_sd, cov_var])?;
            }
        }

        appender.flush().context("Failed to flush Coverage appender")?;
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
                    .filter(|t| t.passed_poisson_test && t.passed_clipping_test)
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

                // Format median clippings as comma-separated integer strings
                let median_left_str: String = pkg.median_left_termini_clippings.iter()
                    .map(|v| format!("{}", v.round() as i32))
                    .collect::<Vec<_>>()
                    .join(",");
                let median_right_str: String = pkg.median_right_termini_clippings.iter()
                    .map(|v| format!("{}", v.round() as i32))
                    .collect::<Vec<_>>()
                    .join(",");

                // Insert into PhageMechanisms
                mechanism_appender.append_row(params![
                    packaging_id,
                    contig_id,
                    sample_id,
                    &pkg.mechanism,
                    pkg.duplication,
                    pkg.repeat_length,
                    terminase_distance,
                    &median_left_str,
                    &median_right_str
                ])?;

                // Insert left termini (status = "start")
                for terminus in &pkg.left_termini {
                    // Store tau as integer (×100)
                    let tau_int = (terminus.tau * 100.0).round() as i32;
                    let passed_poisson_str = if terminus.passed_poisson_test { "yes" } else { "no" };
                    let passed_clipping_str = if terminus.passed_clipping_test { "yes" } else { "no" };
                    let size = terminus.size;
                    termini_appender.append_row(params![
                        terminus_id,
                        packaging_id,
                        terminus.start_pos,
                        terminus.end_pos,
                        size,
                        terminus.center_pos,
                        "start",
                        terminus.total_spc as i32,
                        terminus.median_clippings.round() as i32,
                        terminus.coverage as i32,
                        tau_int,
                        terminus.number_peaks as i32,
                        passed_poisson_str,
                        terminus.expected_spc.round() as i64,
                        compact_pvalue(terminus.pvalue),
                        compact_pvalue(terminus.adjusted_pvalue),
                        passed_clipping_str,
                        terminus.sum_clippings as i64,
                        (terminus.clipped_ratio * 100.0).round() as i32,
                        terminus.expected_clippings.round() as i64
                    ])?;
                    terminus_id += 1;
                }

                // Insert right termini (status = "end")
                for terminus in &pkg.right_termini {
                    let tau_int = (terminus.tau * 100.0).round() as i32;
                    let passed_poisson_str = if terminus.passed_poisson_test { "yes" } else { "no" };
                    let passed_clipping_str = if terminus.passed_clipping_test { "yes" } else { "no" };
                    let size = terminus.size;
                    termini_appender.append_row(params![
                        terminus_id,
                        packaging_id,
                        terminus.start_pos,
                        terminus.end_pos,
                        size,
                        terminus.center_pos,
                        "end",
                        terminus.total_spc as i32,
                        terminus.median_clippings.round() as i32,
                        terminus.coverage as i32,
                        tau_int,
                        terminus.number_peaks as i32,
                        passed_poisson_str,
                        terminus.expected_spc.round() as i64,
                        compact_pvalue(terminus.pvalue),
                        compact_pvalue(terminus.adjusted_pvalue),
                        passed_clipping_str,
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

    fn write_misassembly(&self, conn: &Connection, sample_id: i64, data: &[MisassemblyData]) -> Result<()> {
        let mut appender = conn.appender("Misassembly")
            .context("Failed to create Misassembly appender")?;

        for d in data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&d.contig_name) {
                appender.append_row(params![
                    contig_id, sample_id,
                    d.mismatches_count, d.deletions_count, d.insertions_count, d.clippings_count,
                    d.collapse_bp, d.expansion_bp
                ])?;
            }
        }

        appender.flush().context("Failed to flush Misassembly appender")?;
        Ok(())
    }

    fn write_microdiversity(&self, conn: &Connection, sample_id: i64, data: &[MicrodiversityData]) -> Result<()> {
        let mut appender = conn.appender("Microdiversity")
            .context("Failed to create Microdiversity appender")?;

        for d in data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&d.contig_name) {
                appender.append_row(params![
                    contig_id, sample_id,
                    d.mismatches_count, d.deletions_count, d.insertions_count, d.clippings_count,
                    d.microdiverse_bp_on_reference, d.microdiverse_bp_on_reads
                ])?;
            }
        }

        appender.flush().context("Failed to flush Microdiversity appender")?;
        Ok(())
    }

    fn write_side_misassembly(&self, conn: &Connection, sample_id: i64, data: &[SideMisassemblyData]) -> Result<()> {
        let mut appender = conn.appender("Side_misassembly")
            .context("Failed to create Side_misassembly appender")?;

        for d in data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&d.contig_name) {
                let misjoint = d.contig_end_misjoint_mates.map(|v| v as i64);
                appender.append_row(params![
                    contig_id, sample_id,
                    d.coverage_first_position as i64,
                    d.contig_start_collapse_percentage, d.contig_start_collapse_bp, d.contig_start_expansion_bp,
                    d.coverage_last_position as i64,
                    d.contig_end_collapse_percentage, d.contig_end_collapse_bp, d.contig_end_expansion_bp,
                    misjoint
                ])?;
            }
        }

        appender.flush().context("Failed to flush Side_misassembly appender")?;
        Ok(())
    }

    fn write_topology(&self, conn: &Connection, sample_id: i64, data: &[TopologyData]) -> Result<()> {
        let mut appender = conn.appender("Topology")
            .context("Failed to create Topology appender")?;

        for d in data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&d.contig_name) {
                let circ_reads = d.circularising_reads.map(|v| v as i64);
                let circ_pct = d.circularising_reads_percentage;
                let circ_inserts = d.circularising_inserts.map(|v| v as i64);
                let circ_dev = d.circularising_insert_size_deviation;
                appender.append_row(params![
                    contig_id, sample_id,
                    circ_reads, circ_pct, circ_inserts, circ_dev
                ])?;
            }
        }

        appender.flush().context("Failed to flush Topology appender")?;
        Ok(())
    }

    fn write_features(&self, conn: &Connection, sample_id: i64, features: &[FeaturePoint]) -> Result<()> {
        // Features that have statistics columns
        let features_with_stats = ["left_clippings", "right_clippings", "insertions", "reads_starts", "reads_ends"];

        // Group features by variable name for efficient batch inserts
        let mut by_variable: HashMap<&str, Vec<&FeaturePoint>> = HashMap::new();
        for f in features {
            by_variable.entry(&f.feature).or_default().push(f);
        }

        // Lock created_tables for lazy table creation
        let mut created_tables = self.created_tables.lock().unwrap();

        for v in VARIABLES {
            // Skip contig-level variables stored in Contig_* tables (not per-sample)
            if v.name == "direct_repeat_count" || v.name == "inverted_repeat_count" || v.name == "direct_repeat_identity" || v.name == "inverted_repeat_identity" || v.name == "gc_content" {
                continue;
            }

            let Some(feature_points) = by_variable.get(v.name) else {
                continue;
            };

            let table_name = feature_table_name(v.name);
            let has_stats = features_with_stats.contains(&v.name);

            // Create table lazily if it doesn't exist yet
            create_feature_table_if_needed(conn, &table_name, has_stats, &mut created_tables)?;

            // Determine if this variable stores scaled integers (now only mapq, since tau is removed)
            let is_scaled = v.name == "mapq";



            // Use Appender for bulk loading
            let mut appender = conn.appender(&table_name)
                .with_context(|| format!("Failed to create appender for {}", table_name))?;

            let has_sequences = FEATURES_WITH_SEQUENCES.contains(&table_name.as_str());
            let has_codons = FEATURES_WITH_CODONS.contains(&table_name.as_str());
            let is_single_pos = SINGLE_POSITION_FEATURES.contains(&table_name.as_str());

            if has_stats && has_sequences {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let mean = f.mean.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let median = f.median.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let std = f.std.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let last_pos: Option<i32> = if is_single_pos { None } else { Some(f.end_pos) };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, last_pos, value, mean, median, std, &f.sequence, f.sequence_prevalence])?;

                    }
                }
            } else if has_stats {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let mean = f.mean.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let median = f.median.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let std = f.std.map(|v| if is_scaled { (v * 100.0).round() as i32 } else { v.round() as i32 });
                        let last_pos: Option<i32> = if is_single_pos { None } else { Some(f.end_pos) };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, last_pos, value, mean, median, std])?;

                    }
                }
            } else if has_sequences && has_codons {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let last_pos: Option<i32> = if is_single_pos { None } else { Some(f.end_pos) };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, last_pos, value, &f.sequence, f.sequence_prevalence, &f.codon_category, &f.codon_change, &f.aa_change])?;

                    }
                }
            } else if has_sequences {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let last_pos: Option<i32> = if is_single_pos { None } else { Some(f.end_pos) };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, last_pos, value, &f.sequence, f.sequence_prevalence])?;

                    }
                }
            } else {
                for f in feature_points {
                    if let Some(&contig_id) = self.contig_name_to_id.get(&f.contig_name) {
                        let value = if is_scaled { (f.value * 100.0).round() as i32 } else { f.value.round() as i32 };
                        let last_pos: Option<i32> = if is_single_pos { None } else { Some(f.end_pos) };
                        appender.append_row(params![contig_id, sample_id, f.start_pos, last_pos, value])?;

                    }
                }
            }

            appender.flush()
                .with_context(|| format!("Failed to flush appender for {}", table_name))?;

        }

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
            let percentage = ((total_bp as f64 / contig_length as f64) * 1000.0).round() as i32;

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
                (data.skew_stats.percent_positive * 10.0).round() as i32, // Positive_GC_skew_windows_percentage as int (0-1000, ×10)
                &data.contig_name
            ])?;
        }

        Ok(())
    }

    /// Finalize the database after all samples are processed.
    pub fn finalize(self) -> Result<()> {
        let conn = self.conn.into_inner().unwrap();
        let created_tables = self.created_tables.into_inner().unwrap();

        // Create derived views (in extend mode, views were dropped in open())
        create_views(&conn, self.has_bam)?;

        if !self.is_extend {
            // Delete Variable entries for features without tables (only for fresh DBs)
            cleanup_unused_variables(&conn, &created_tables)?;
        }

        // Drop empty module tables and their views (prevents empty UI sections)
        if self.has_bam {
            drop_empty_tables(&conn)?;
        }

        // Force checkpoint to compress data and write to disk
        conn.execute("CHECKPOINT", [])
            .context("Failed to checkpoint database")?;

        Ok(())
    }

    /// Get the contig ID for a contig name.
    pub fn get_contig_id(&self, name: &str) -> Option<i64> {
        self.contig_name_to_id.get(name).copied()
    }

    /// Get all contig names currently in the database.
    pub fn contig_names(&self) -> Vec<String> {
        self.contig_name_to_id.keys().cloned().collect()
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

/// Compact a p-value into a short text representation with 2 significant digits.
///
/// Examples:
/// - 3.45e-5  → "35e-6"  (3.45 rounded to 35, exponent shifted by -1)
/// - 1.5e-10  → "15e-11"
/// - 1.0      → "10e-1"
/// - 0.0      → "0"
fn compact_pvalue(p: f64) -> String {
    if p == 0.0 {
        return "0".to_string();
    }
    if p.is_nan() || p.is_infinite() {
        return "0".to_string();
    }

    // Get the exponent (floor of log10)
    let log10 = p.log10();
    let exponent = log10.floor() as i32;

    // Get the significand as a 2-digit integer
    // p = significand * 10^exponent, where significand is in [1.0, 10.0)
    // We want a 2-digit integer, so multiply by 10
    let significand = p / 10.0_f64.powi(exponent);
    let two_digits = (significand * 10.0).round() as i32;

    // The effective exponent shifts by -1 since we multiplied significand by 10
    let effective_exponent = exponent - 1;

    format!("{}e{}", two_digits, effective_exponent)
}

/// Create core database tables (no indexes - DuckDB uses zone maps).
fn create_core_tables(conn: &Connection, has_bam: bool) -> Result<()> {
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

    if has_bam {
    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY,
            Sample_name TEXT UNIQUE,
            Sequencing_type TEXT,
            Number_of_reads INTEGER,
            Number_of_mapped_reads INTEGER,
            Circular_mapping BOOLEAN
        )",
        [],
    )
    .context("Failed to create Sample table")?;

    // Coverage table
    conn.execute(
        "CREATE TABLE Coverage (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Aligned_fraction_percentage INTEGER,
            Above_expected_aligned_fraction BOOLEAN,
            Read_count INTEGER,
            Coverage_mean INTEGER,
            Coverage_median INTEGER,
            Coverage_trimmed_mean INTEGER,
            Coverage_sd INTEGER,
            Coverage_variation INTEGER,
            PRIMARY KEY (Contig_id, Sample_id)
        )",
        [],
    )
    .context("Failed to create Coverage table")?;

    conn.execute(
        "CREATE TABLE PhageMechanisms (
            Packaging_id INTEGER PRIMARY KEY,
            Contig_id INTEGER NOT NULL,
            Sample_id INTEGER NOT NULL,
            Packaging_mechanism TEXT NOT NULL,
            Duplication BOOLEAN,
            Repeat_length INTEGER,
            Terminase_distance INTEGER,
            Median_left_termini_clippings TEXT,
            Median_right_termini_clippings TEXT
        )",
        [],
    )
    .context("Failed to create PhageMechanisms table")?;

    // PhageTermini table - stores individual terminus areas with full metadata
    // One row per terminus area, linked to PhageMechanisms via Packaging_id
    // Includes filtering diagnostics for both Poisson and clipping tests
    conn.execute(
        "CREATE TABLE PhageTermini (
            Terminus_id INTEGER PRIMARY KEY,
            Packaging_id INTEGER NOT NULL,
            \"Start\" INTEGER NOT NULL,
            \"End\" INTEGER NOT NULL,
            \"Size\" INTEGER NOT NULL,
            Center INTEGER NOT NULL,
            Status TEXT NOT NULL,
            SPC INTEGER NOT NULL,
            Median_clippings INTEGER NOT NULL,
            Coverage INTEGER NOT NULL,
            Tau INTEGER NOT NULL,
            NumberPeaks INTEGER NOT NULL,
            Passed_PoissonTest TEXT NOT NULL,
            Expected_SPC INTEGER NOT NULL,
            Pvalue TEXT NOT NULL,
            Adjusted_pvalue TEXT NOT NULL,
            Passed_ClippingTest TEXT NOT NULL,
            Clippings INTEGER NOT NULL,
            Clipping_excess INTEGER NOT NULL,
            Expected_clippings INTEGER NOT NULL
        )",
        [],
    )
    .context("Failed to create PhageTermini table")?;
    } // end if has_bam (Sample, Coverage, PhageMechanisms, PhageTermini)

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

    if has_bam {
    // Misassembly table - counts at ≥50% prevalence threshold
    conn.execute(
        "CREATE TABLE Misassembly (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Mismatches_count INTEGER,
            Deletions_count INTEGER,
            Insertions_count INTEGER,
            Clippings_count INTEGER,
            Collapse_bp INTEGER,
            Expansion_bp INTEGER
        )",
        [],
    )
    .context("Failed to create Misassembly table")?;

    // Microdiversity table - counts at ≥10% prevalence threshold
    conn.execute(
        "CREATE TABLE Microdiversity (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Mismatches_count INTEGER,
            Deletions_count INTEGER,
            Insertions_count INTEGER,
            Clippings_count INTEGER,
            Microdiverse_bp_on_reference INTEGER,
            Microdiverse_bp_on_reads INTEGER
        )",
        [],
    )
    .context("Failed to create Microdiversity table")?;

    // Side_misassembly table - left/right clipping events with ≥50% prevalence
    conn.execute(
        "CREATE TABLE Side_misassembly (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Coverage_first_position INTEGER,
            Contig_start_collapse_percentage INTEGER,
            Contig_start_collapse_bp INTEGER,
            Contig_start_expansion_bp INTEGER,
            Coverage_last_position INTEGER,
            Contig_end_collapse_percentage INTEGER,
            Contig_end_collapse_bp INTEGER,
            Contig_end_expansion_bp INTEGER,
            Contig_end_misjoint_mates INTEGER
        )",
        [],
    )
    .context("Failed to create Side_misassembly table")?;

    // Topology table - circularisation metrics
    conn.execute(
        "CREATE TABLE Topology (
            Contig_id INTEGER,
            Sample_id INTEGER,
            Circularising_reads INTEGER,
            Circularising_reads_percentage INTEGER,
            Circularising_inserts INTEGER,
            Circularising_insert_size_deviation INTEGER
        )",
        [],
    )
    .context("Failed to create Topology table")?;
    } // end if has_bam (Misassembly, Microdiversity, Side_misassembly, Topology)

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
            Longest_isoform BOOLEAN,
            Nucleotide_sequence TEXT,
            Protein_sequence TEXT
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

    // Contig_sequence table - stores full contig sequences for visualization
    conn.execute(
        "CREATE TABLE Contig_sequence (
            Contig_id INTEGER PRIMARY KEY,
            Sequence TEXT NOT NULL
        )",
        [],
    )
    .context("Failed to create Contig_sequence table")?;

    // Codon_table - standard genetic code lookup (64 codons)
    conn.execute(
        "CREATE TABLE Codon_table (
            Codon TEXT PRIMARY KEY,
            AminoAcid TEXT NOT NULL,
            AminoAcid_name TEXT NOT NULL,
            Color TEXT NOT NULL
        )",
        [],
    )
    .context("Failed to create Codon_table table")?;

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

/// Insert contig sequences into the database.
fn insert_contig_sequences(conn: &Connection, contigs: &[ContigInfo]) -> Result<()> {
    let mut appender = conn.appender("Contig_sequence")
        .context("Failed to create Contig_sequence appender")?;

    let mut count = 0;
    for (i, contig) in contigs.iter().enumerate() {
        if let Some(ref seq) = contig.sequence {
            let seq_str = String::from_utf8_lossy(seq).into_owned();
            appender.append_row(params![(i + 1) as i64, &seq_str])?;
            count += 1;
        }
    }

    appender.flush().context("Failed to flush Contig_sequence appender")?;
    eprintln!("Stored sequences for {}/{} contigs", count, contigs.len());
    Ok(())
}

/// Insert contig sequences for new contigs with an ID offset (for extend mode).
fn insert_contig_sequences_with_offset(conn: &Connection, contigs: &[ContigInfo], id_offset: i64) -> Result<()> {
    let mut appender = conn.appender("Contig_sequence")
        .context("Failed to create Contig_sequence appender")?;

    let mut count = 0;
    for (i, contig) in contigs.iter().enumerate() {
        if let Some(ref seq) = contig.sequence {
            let seq_str = String::from_utf8_lossy(seq).into_owned();
            appender.append_row(params![id_offset + 1 + i as i64, &seq_str])?;
            count += 1;
        }
    }

    appender.flush().context("Failed to flush Contig_sequence appender")?;
    if count > 0 {
        eprintln!("Stored sequences for {}/{} new contigs", count, contigs.len());
    }
    Ok(())
}

/// Insert the standard genetic code codon table (64 codons).
fn insert_codon_table(conn: &Connection) -> Result<()> {
    let mut appender = conn.appender("Codon_table")
        .context("Failed to create Codon_table appender")?;

    // Standard genetic code: (codon, 1-letter AA, full name, color)
    // Colors by biochemical property:
    //   Hydrophobic aliphatic (G,A,V,L,I,P): #2d6a4f
    //   Aromatic (F,W,Y): #b5838d
    //   Polar uncharged (S,T,N,Q,C,M): #457b9d
    //   Positively charged (K,R,H): #9b2226
    //   Negatively charged (D,E): #e9c46a
    //   Stop (*): #6a3d9a
    let codons: &[(&str, &str, &str, &str)] = &[
        ("TTT", "F", "Phenylalanine", "#b5838d"), ("TTC", "F", "Phenylalanine", "#b5838d"),
        ("TTA", "L", "Leucine", "#2d6a4f"), ("TTG", "L", "Leucine", "#2d6a4f"),
        ("TCT", "S", "Serine", "#457b9d"), ("TCC", "S", "Serine", "#457b9d"),
        ("TCA", "S", "Serine", "#457b9d"), ("TCG", "S", "Serine", "#457b9d"),
        ("TAT", "Y", "Tyrosine", "#b5838d"), ("TAC", "Y", "Tyrosine", "#b5838d"),
        ("TAA", "*", "Stop", "#6a3d9a"), ("TAG", "*", "Stop", "#6a3d9a"),
        ("TGT", "C", "Cysteine", "#457b9d"), ("TGC", "C", "Cysteine", "#457b9d"),
        ("TGA", "*", "Stop", "#6a3d9a"), ("TGG", "W", "Tryptophan", "#b5838d"),
        ("CTT", "L", "Leucine", "#2d6a4f"), ("CTC", "L", "Leucine", "#2d6a4f"),
        ("CTA", "L", "Leucine", "#2d6a4f"), ("CTG", "L", "Leucine", "#2d6a4f"),
        ("CCT", "P", "Proline", "#2d6a4f"), ("CCC", "P", "Proline", "#2d6a4f"),
        ("CCA", "P", "Proline", "#2d6a4f"), ("CCG", "P", "Proline", "#2d6a4f"),
        ("CAT", "H", "Histidine", "#9b2226"), ("CAC", "H", "Histidine", "#9b2226"),
        ("CAA", "Q", "Glutamine", "#457b9d"), ("CAG", "Q", "Glutamine", "#457b9d"),
        ("CGT", "R", "Arginine", "#9b2226"), ("CGC", "R", "Arginine", "#9b2226"),
        ("CGA", "R", "Arginine", "#9b2226"), ("CGG", "R", "Arginine", "#9b2226"),
        ("ATT", "I", "Isoleucine", "#2d6a4f"), ("ATC", "I", "Isoleucine", "#2d6a4f"),
        ("ATA", "I", "Isoleucine", "#2d6a4f"), ("ATG", "M", "Methionine", "#457b9d"),
        ("ACT", "T", "Threonine", "#457b9d"), ("ACC", "T", "Threonine", "#457b9d"),
        ("ACA", "T", "Threonine", "#457b9d"), ("ACG", "T", "Threonine", "#457b9d"),
        ("AAT", "N", "Asparagine", "#457b9d"), ("AAC", "N", "Asparagine", "#457b9d"),
        ("AAA", "K", "Lysine", "#9b2226"), ("AAG", "K", "Lysine", "#9b2226"),
        ("AGT", "S", "Serine", "#457b9d"), ("AGC", "S", "Serine", "#457b9d"),
        ("AGA", "R", "Arginine", "#9b2226"), ("AGG", "R", "Arginine", "#9b2226"),
        ("GTT", "V", "Valine", "#2d6a4f"), ("GTC", "V", "Valine", "#2d6a4f"),
        ("GTA", "V", "Valine", "#2d6a4f"), ("GTG", "V", "Valine", "#2d6a4f"),
        ("GCT", "A", "Alanine", "#2d6a4f"), ("GCC", "A", "Alanine", "#2d6a4f"),
        ("GCA", "A", "Alanine", "#2d6a4f"), ("GCG", "A", "Alanine", "#2d6a4f"),
        ("GAT", "D", "Aspartate", "#e9c46a"), ("GAC", "D", "Aspartate", "#e9c46a"),
        ("GAA", "E", "Glutamate", "#e9c46a"), ("GAG", "E", "Glutamate", "#e9c46a"),
        ("GGT", "G", "Glycine", "#2d6a4f"), ("GGC", "G", "Glycine", "#2d6a4f"),
        ("GGA", "G", "Glycine", "#2d6a4f"), ("GGG", "G", "Glycine", "#2d6a4f"),
    ];

    for &(codon, aa, name, color) in codons {
        appender.append_row(params![codon, aa, name, color])?;
    }

    appender.flush().context("Failed to flush Codon_table appender")?;
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
        if !has_pharokka_function && ann.feature_type == "CDS" {
            if let Some(function) = &ann.function {
                let func_lower = function.to_lowercase();
                for pharokka_key in &pharokka_keys {
                    if func_lower.contains(pharokka_key) {
                        has_pharokka_function = true;
                        break;  // Exit inner pharokka_key loop only
                    }
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
            &ann.nucleotide_sequence,
            &ann.protein_sequence,
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
        // Special case: contig-level data is stored in separate Contig_* tables/views
        let table_name = match v.name {
            "direct_repeat_count" => "Contig_direct_repeat_count".to_string(),
            "direct_repeat_identity" => "Contig_direct_repeat_identity".to_string(),
            "inverted_repeat_count" => "Contig_inverted_repeat_count".to_string(),
            "inverted_repeat_identity" => "Contig_inverted_repeat_identity".to_string(),
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

/// Feature tables that have Sequence and Sequence_prevalence columns.
const FEATURES_WITH_SEQUENCES: &[&str] = &[
    "Feature_mismatches",
    "Feature_insertions",
    "Feature_left_clippings",
    "Feature_right_clippings",
    "Feature_reads_starts",
    "Feature_reads_ends",
];

/// Feature tables that have Codon_category, Codon_change, AA_change columns.
const FEATURES_WITH_CODONS: &[&str] = &[
    "Feature_mismatches",
];

/// Feature tables where Last_position is always equal to First_position (single-position bar spikes).
/// For these tables, Last_position is stored as NULL to save storage space.
/// At read time, COALESCE(Last_position, First_position) recovers the value.
const SINGLE_POSITION_FEATURES: &[&str] = &[
    "Feature_insertions",
    "Feature_mismatches",
    "Feature_left_clippings",
    "Feature_right_clippings",
    "Feature_reads_starts",
    "Feature_reads_ends",
];

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

    let has_sequences = FEATURES_WITH_SEQUENCES.contains(&table_name);
    let has_codons = FEATURES_WITH_CODONS.contains(&table_name);

    let table_sql = if has_stats && has_sequences {
        format!(
            "CREATE TABLE IF NOT EXISTS {} (
                Contig_id INTEGER,
                Sample_id INTEGER,
                First_position INTEGER,
                Last_position INTEGER,
                Value INTEGER,
                Mean INTEGER,
                Median INTEGER,
                Std INTEGER,
                Sequence TEXT,
                Sequence_prevalence INTEGER
            )",
            table_name
        )
    } else if has_stats {
        format!(
            "CREATE TABLE IF NOT EXISTS {} (
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
    } else if has_sequences && has_codons {
        format!(
            "CREATE TABLE IF NOT EXISTS {} (
                Contig_id INTEGER,
                Sample_id INTEGER,
                First_position INTEGER,
                Last_position INTEGER,
                Value INTEGER,
                Sequence TEXT,
                Sequence_prevalence INTEGER,
                Codon_category TEXT,
                Codon_change TEXT,
                AA_change TEXT
            )",
            table_name
        )
    } else if has_sequences {
        format!(
            "CREATE TABLE IF NOT EXISTS {} (
                Contig_id INTEGER,
                Sample_id INTEGER,
                First_position INTEGER,
                Last_position INTEGER,
                Value INTEGER,
                Sequence TEXT,
                Sequence_prevalence INTEGER
            )",
            table_name
        )
    } else {
        format!(
            "CREATE TABLE IF NOT EXISTS {} (
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

/// Create materialized tables and views after all data is inserted.
fn create_views(conn: &Connection, has_bam: bool) -> Result<()> {
    // Feature_primary_reads is now written directly as a table from Rust processing data
    // (no longer reconstructed from plus+minus strands via SQL boundary sweep).

    // Materialized repeat tables (sweep-line algorithm, computed once after all repeat data is inserted).
    // These produce (Contig_id, First_position, Last_position, Value) like other Contig_* tables,
    // allowing the standard binning logic in get_feature_data() to work with repeats.

    conn.execute(
        "CREATE TABLE Contig_direct_repeat_count AS
        WITH boundaries AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_directRepeats
            UNION
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_directRepeats
        ),
        segments AS (
            SELECT Contig_id, pos AS seg_start,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
            FROM boundaries
        )
        SELECT s.Contig_id,
               s.seg_start AS First_position,
               s.seg_end AS Last_position,
               (SELECT COUNT(*) FROM Contig_directRepeats r
                WHERE r.Contig_id = s.Contig_id
                  AND LEAST(r.Position1, r.Position2) <= s.seg_end
                  AND GREATEST(r.Position1, r.Position2) >= s.seg_start) AS Value
        FROM segments s
        WHERE s.seg_end IS NOT NULL",
        [],
    )
    .context("Failed to create Contig_direct_repeat_count table")?;

    conn.execute(
        "CREATE TABLE Contig_inverted_repeat_count AS
        WITH boundaries AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_invertedRepeats
            UNION
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_invertedRepeats
        ),
        segments AS (
            SELECT Contig_id, pos AS seg_start,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
            FROM boundaries
        )
        SELECT s.Contig_id,
               s.seg_start AS First_position,
               s.seg_end AS Last_position,
               (SELECT COUNT(*) FROM Contig_invertedRepeats r
                WHERE r.Contig_id = s.Contig_id
                  AND LEAST(r.Position1, r.Position2) <= s.seg_end
                  AND GREATEST(r.Position1, r.Position2) >= s.seg_start) AS Value
        FROM segments s
        WHERE s.seg_end IS NOT NULL",
        [],
    )
    .context("Failed to create Contig_inverted_repeat_count table")?;

    conn.execute(
        "CREATE TABLE Contig_direct_repeat_identity AS
        WITH boundaries AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_directRepeats
            UNION
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_directRepeats
        ),
        segments AS (
            SELECT Contig_id, pos AS seg_start,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
            FROM boundaries
        )
        SELECT s.Contig_id,
               s.seg_start AS First_position,
               s.seg_end AS Last_position,
               COALESCE((SELECT MAX(r.Pident) / 100.0 FROM Contig_directRepeats r
                WHERE r.Contig_id = s.Contig_id
                  AND LEAST(r.Position1, r.Position2) <= s.seg_end
                  AND GREATEST(r.Position1, r.Position2) >= s.seg_start), 0) AS Value
        FROM segments s
        WHERE s.seg_end IS NOT NULL",
        [],
    )
    .context("Failed to create Contig_direct_repeat_identity table")?;

    conn.execute(
        "CREATE TABLE Contig_inverted_repeat_identity AS
        WITH boundaries AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_invertedRepeats
            UNION
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_invertedRepeats
        ),
        segments AS (
            SELECT Contig_id, pos AS seg_start,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
            FROM boundaries
        )
        SELECT s.Contig_id,
               s.seg_start AS First_position,
               s.seg_end AS Last_position,
               COALESCE((SELECT MAX(r.Pident) / 100.0 FROM Contig_invertedRepeats r
                WHERE r.Contig_id = s.Contig_id
                  AND LEAST(r.Position1, r.Position2) <= s.seg_end
                  AND GREATEST(r.Position1, r.Position2) >= s.seg_start), 0) AS Value
        FROM segments s
        WHERE s.seg_end IS NOT NULL",
        [],
    )
    .context("Failed to create Contig_inverted_repeat_identity table")?;

    if has_bam {
    // Explicit_coverage VIEW with RPKM and TPM
    conn.execute(
        "CREATE VIEW Explicit_coverage AS
         WITH rpkm_base AS (
             SELECT
                 p.Contig_id,
                 p.Sample_id,
                 CASE WHEN c.Contig_length > 0 AND s.Number_of_mapped_reads > 0
                      THEN (CAST(p.Read_count AS DOUBLE) * 1e9) / (CAST(c.Contig_length AS DOUBLE) * CAST(s.Number_of_mapped_reads AS DOUBLE))
                      ELSE 0 END AS RPKM
             FROM Coverage p
             JOIN Contig c ON p.Contig_id = c.Contig_id
             JOIN Sample s ON p.Sample_id = s.Sample_id
         ),
         rpkm_sum AS (
             SELECT Sample_id, SUM(RPKM) AS total_rpkm FROM rpkm_base GROUP BY Sample_id
         )
         SELECT
             c.Contig_name,
             s.Sample_name,
             p.Aligned_fraction_percentage / 10.0 AS Aligned_fraction_percentage,
             p.Above_expected_aligned_fraction,
             p.Read_count,
             p.Coverage_mean / 10.0 AS Coverage_mean,
             p.Coverage_median / 10.0 AS Coverage_median,
             p.Coverage_trimmed_mean / 10.0 AS Coverage_trimmed_mean,
             rb.RPKM,
             CASE WHEN rs.total_rpkm > 0 THEN (rb.RPKM / rs.total_rpkm) * 1e6 ELSE 0 END AS TPM,
             ROUND(p.Coverage_sd / 1000000.0, 2) AS Coverage_sd,
             ROUND(p.Coverage_variation / 1000000.0, 4) AS Coverage_variation
         FROM Coverage p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id
         JOIN rpkm_base rb ON p.Contig_id = rb.Contig_id AND p.Sample_id = rb.Sample_id
         JOIN rpkm_sum rs ON p.Sample_id = rs.Sample_id",
        [],
    )
    .context("Failed to create Explicit_coverage VIEW")?;

    // Explicit_misassembly VIEW - per 100kbp normalization
    conn.execute(
        "CREATE VIEW Explicit_misassembly AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             CASE WHEN c.Contig_length > 0 THEN m.Mismatches_count * 100000.0 / c.Contig_length ELSE 0 END AS Mismatches_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN m.Deletions_count * 100000.0 / c.Contig_length ELSE 0 END AS Deletions_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN m.Insertions_count * 100000.0 / c.Contig_length ELSE 0 END AS Insertions_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN m.Clippings_count * 100000.0 / c.Contig_length ELSE 0 END AS Clippings_per_100kbp,
             m.Collapse_bp,
             CASE WHEN c.Contig_length > 0 THEN m.Collapse_bp * 100000.0 / c.Contig_length ELSE 0 END AS Collapse_per_100kbp,
             m.Expansion_bp,
             CASE WHEN c.Contig_length > 0 THEN m.Expansion_bp * 100000.0 / c.Contig_length ELSE 0 END AS Expansion_per_100kbp
         FROM Misassembly m
         JOIN Contig c ON m.Contig_id = c.Contig_id
         JOIN Sample s ON m.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_misassembly VIEW")?;

    // Explicit_microdiversity VIEW
    conn.execute(
        "CREATE VIEW Explicit_microdiversity AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             CASE WHEN c.Contig_length > 0 THEN md.Mismatches_count * 100000.0 / c.Contig_length ELSE 0 END AS Mismatches_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN md.Deletions_count * 100000.0 / c.Contig_length ELSE 0 END AS Deletions_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN md.Insertions_count * 100000.0 / c.Contig_length ELSE 0 END AS Insertions_per_100kbp,
             CASE WHEN c.Contig_length > 0 THEN md.Clippings_count * 100000.0 / c.Contig_length ELSE 0 END AS Clippings_per_100kbp,
             md.Microdiverse_bp_on_reads,
             CASE WHEN c.Contig_length > 0 THEN md.Microdiverse_bp_on_reads * 100000.0 / c.Contig_length ELSE 0 END AS Microdiverse_bp_per_100kbp_on_reads,
             md.Microdiverse_bp_on_reference,
             CASE WHEN c.Contig_length > 0 THEN md.Microdiverse_bp_on_reference * 100000.0 / c.Contig_length ELSE 0 END AS Microdiverse_bp_per_100kbp_on_reference
         FROM Microdiversity md
         JOIN Contig c ON md.Contig_id = c.Contig_id
         JOIN Sample s ON md.Sample_id = s.Sample_id",
        [],
    )
    .context("Failed to create Explicit_microdiversity VIEW")?;

    // Explicit_side_misassembly VIEW (JOINs with Coverage for normalization)
    conn.execute(
        "CREATE VIEW Explicit_side_misassembly AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             COALESCE(sm.Coverage_first_position, 0) AS Coverage_first_position,
             COALESCE(sm.Contig_start_collapse_percentage, 0) / 10.0 AS Contig_start_collapse_prevalence,
             COALESCE(sm.Contig_start_collapse_bp, 0) AS Contig_start_collapse_bp,
             COALESCE(sm.Contig_start_expansion_bp, 0) AS Contig_start_expansion_bp,
             COALESCE(sm.Coverage_last_position, 0) AS Coverage_last_position,
             COALESCE(sm.Contig_end_collapse_percentage, 0) / 10.0 AS Contig_end_collapse_prevalence,
             COALESCE(sm.Contig_end_collapse_bp, 0) AS Contig_end_collapse_bp,
             COALESCE(sm.Contig_end_expansion_bp, 0) AS Contig_end_expansion_bp,
             CASE WHEN s.Sequencing_type = 'paired-short' THEN COALESCE(sm.Contig_end_misjoint_mates, 0) ELSE NULL END AS Contig_end_misjoint_mates,
             CASE WHEN s.Sequencing_type != 'paired-short' THEN NULL WHEN cov.Coverage_mean > 0 THEN COALESCE(sm.Contig_end_misjoint_mates, 0) * 1000.0 / cov.Coverage_mean ELSE 0 END AS Normalized_contig_end_misjoint_mates
         FROM Side_misassembly sm
         JOIN Contig c ON sm.Contig_id = c.Contig_id
         JOIN Sample s ON sm.Sample_id = s.Sample_id
         LEFT JOIN Coverage cov ON sm.Contig_id = cov.Contig_id AND sm.Sample_id = cov.Sample_id",
        [],
    )
    .context("Failed to create Explicit_side_misassembly VIEW")?;

    // Explicit_topology VIEW (JOINs with Coverage for normalization)
    conn.execute(
        "CREATE VIEW Explicit_topology AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             COALESCE(t.Circularising_reads, 0) AS Circularising_reads,
             COALESCE(t.Circularising_reads_percentage, 0) AS Circularising_reads_prevalence,
             CASE WHEN s.Sequencing_type = 'paired-short' THEN COALESCE(t.Circularising_inserts, 0) ELSE NULL END AS Circularising_inserts,
             CASE WHEN s.Sequencing_type = 'paired-short' THEN COALESCE(t.Circularising_insert_size_deviation, 0) ELSE NULL END AS Circularising_insert_size_deviation,
             CASE WHEN s.Sequencing_type != 'paired-short' THEN NULL WHEN cov.Coverage_mean > 0 THEN COALESCE(t.Circularising_inserts, 0) * 1000.0 / cov.Coverage_mean ELSE 0 END AS Normalized_circularising_inserts
         FROM Topology t
         JOIN Contig c ON t.Contig_id = c.Contig_id
         JOIN Sample s ON t.Sample_id = s.Sample_id
         LEFT JOIN Coverage cov ON t.Contig_id = cov.Contig_id AND t.Sample_id = cov.Sample_id",
        [],
    )
    .context("Failed to create Explicit_topology VIEW")?;

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
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'start' AND pt.Passed_PoissonTest = 'yes' AND pt.Passed_ClippingTest = 'yes'),
                 ''
             ) AS Left_termini,
             m.Median_left_termini_clippings,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM PhageTermini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'end' AND pt.Passed_PoissonTest = 'yes' AND pt.Passed_ClippingTest = 'yes'),
                 ''
             ) AS Right_termini,
             m.Median_right_termini_clippings,
             CASE WHEN m.Duplication = true THEN 'DTR'
                  WHEN m.Duplication = false THEN 'ITR'
                  ELSE NULL END AS Duplication,
             COALESCE(
                 (SELECT COUNT(*)
                  FROM PhageTermini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Passed_PoissonTest = 'yes' AND pt.Passed_ClippingTest = 'yes'),
                 0
             ) AS Total_peaks,
             m.Repeat_length,
             m.Terminase_distance,
             CASE WHEN m.Terminase_distance IS NOT NULL AND c.Contig_length > 0
                  THEN CAST(ROUND(CAST(m.Terminase_distance AS REAL) / c.Contig_length * 100) AS INTEGER)
                  ELSE NULL END AS Terminase_percentage
         FROM Coverage p
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
             t.\"Size\",
             t.Center,
             t.Status,
             t.SPC,
             t.Median_clippings,
             t.Coverage,
             t.Tau,
             t.NumberPeaks,
             t.Passed_PoissonTest,
             t.Expected_SPC,
             t.Pvalue,
             t.Adjusted_pvalue,
             t.Passed_ClippingTest,
             t.Clippings,
             t.Clipping_excess,
             t.Expected_clippings
         FROM PhageTermini t
         JOIN PhageMechanisms m ON t.Packaging_id = m.Packaging_id
         JOIN Contig c ON m.Contig_id = c.Contig_id
         JOIN Sample s ON m.Sample_id = s.Sample_id",
    )
    .context("Failed to create Explicit_phage_termini VIEW")?;
    } // end if has_bam (Explicit_* views)

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
        if var_name == "direct_repeat_count" || var_name == "inverted_repeat_count" || var_name == "direct_repeat_identity" || var_name == "inverted_repeat_identity" || var_name == "gc_content" || var_name == "gc_skew" {
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

/// Drop empty module-dependent tables and their views.
/// Prevents empty sections in the Filtering/Summary UI when modules weren't calculated.
fn drop_empty_tables(conn: &Connection) -> Result<()> {
    let table_view_map: &[(&str, &[&str])] = &[
        ("Misassembly", &["Explicit_misassembly"]),
        ("Microdiversity", &["Explicit_microdiversity"]),
        ("Side_misassembly", &["Explicit_side_misassembly"]),
        ("Topology", &["Explicit_topology"]),
        ("PhageTermini", &["Explicit_phage_termini"]),
        ("PhageMechanisms", &["Explicit_phage_mechanisms"]),
    ];

    for (table, views) in table_view_map {
        let count: i64 = conn
            .query_row(&format!("SELECT COUNT(*) FROM {}", table), [], |row| row.get(0))
            .unwrap_or(0);
        if count == 0 {
            for view in *views {
                conn.execute(&format!("DROP VIEW IF EXISTS {}", view), [])?;
            }
            conn.execute(&format!("DROP TABLE IF EXISTS {}", table), [])?;
        }
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

/// Misassembly data for a contig (≥50% prevalence threshold).
#[derive(Clone, Debug)]
pub struct MisassemblyData {
    pub contig_name: String,
    pub mismatches_count: i64,
    pub deletions_count: i64,
    pub insertions_count: i64,
    pub clippings_count: i64,
    pub collapse_bp: i64,
    pub expansion_bp: i64,
}

/// Microdiversity data for a contig (≥10% prevalence threshold).
#[derive(Clone, Debug)]
pub struct MicrodiversityData {
    pub contig_name: String,
    pub mismatches_count: i64,
    pub deletions_count: i64,
    pub insertions_count: i64,
    pub clippings_count: i64,
    pub microdiverse_bp_on_reference: i64,
    pub microdiverse_bp_on_reads: i64,
}

/// Side misassembly data for a contig (left/right clipping events with ≥50% prevalence).
#[derive(Clone, Debug)]
pub struct SideMisassemblyData {
    pub contig_name: String,
    pub coverage_first_position: u64,
    pub contig_start_collapse_percentage: Option<i32>,
    pub contig_start_collapse_bp: Option<i32>,
    pub contig_start_expansion_bp: Option<i32>,
    pub coverage_last_position: u64,
    pub contig_end_collapse_percentage: Option<i32>,
    pub contig_end_collapse_bp: Option<i32>,
    pub contig_end_expansion_bp: Option<i32>,
    pub contig_end_misjoint_mates: Option<u64>,
}

/// Topology data for a contig (circularisation metrics).
#[derive(Clone, Debug)]
pub struct TopologyData {
    pub contig_name: String,
    pub circularising_reads: Option<u64>,
    pub circularising_reads_percentage: Option<i32>,
    pub circularising_inserts: Option<u64>,
    pub circularising_insert_size_deviation: Option<i32>,
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
