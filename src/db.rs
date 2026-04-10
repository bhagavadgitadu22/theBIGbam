//! DuckDB database operations.
//!
//! Database schema:
//! - `Contig`: Contig metadata (name, length, annotation tool)
//! - `Sample`: Sample names from BAM files
//! - `Coverage`: Coverage metrics per contig per sample
//! - `Contig_annotation`: Gene annotations from GenBank
//! - `Variable`: Feature metadata (name, type, plot configuration)
//! - `Feature_blob`: Per-sample feature data as compressed binary blobs
//! - `Contig_blob`: Contig-level feature data as compressed binary blobs
//!
//! DuckDB advantages over SQLite:
//! - Columnar storage with automatic compression
//! - No need for explicit indexes (uses zone maps)
//! - Better write performance with batch inserts
//! - Simpler concurrency model (single writer, multiple readers)

use anyhow::{Context, Result};
use chrono::Local;
use duckdb::{params, Connection};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;

use crate::features::translate_codon;
use crate::gc_content::{GCSkewStats, GCStats};
use crate::types::{feature_name_to_id, ContigInfo, FeatureAnnotation, PackagingData, PresenceData, VARIABLES};
// Re-export new metric data structs (defined below in this file)

/// Thread-safe database connection wrapper for sequential writes.
pub struct DbWriter {
    conn: Mutex<Connection>,
    has_bam: bool,
    contig_name_to_id: HashMap<String, i64>,
    sample_name_to_id: Mutex<HashMap<String, i64>>,
    next_sample_id: Mutex<i64>,
}

impl DbWriter {
    /// Create a new database and return a writer for sequential sample insertion.
    pub fn create(
        db_path: &Path,
        contigs: &[ContigInfo],
        annotations: &[FeatureAnnotation],
        contig_qualifiers: &[(i64, HashMap<String, String>)],
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

        // Insert annotations and their qualifiers
        insert_annotations(&conn, annotations)?;

        // Insert contig-level qualifiers from source features
        insert_contig_qualifiers(&conn, contig_qualifiers)?;

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
            contig_name_to_id,
            sample_name_to_id: Mutex::new(HashMap::new()),
            next_sample_id: Mutex::new(1),
        })
    }

    /// Open an existing database for extending with new samples/contigs.
    pub fn open(
        db_path: &Path,
        new_contigs: &[ContigInfo],
        new_annotations: &[FeatureAnnotation],
        new_contig_qualifiers: &[(i64, HashMap<String, String>)],
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

            // Insert contig-level qualifiers (remapped to new contig IDs)
            if !new_contig_qualifiers.is_empty() {
                let remapped_quals: Vec<(i64, HashMap<String, String>)> = new_contig_qualifiers
                    .iter()
                    .map(|(cid, quals)| (cid + max_contig_id, quals.clone()))
                    .collect();
                insert_contig_qualifiers(&conn, &remapped_quals)?;
            }

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
            contig_name_to_id,
            sample_name_to_id: Mutex::new(sample_name_to_id),
            next_sample_id: Mutex::new(max_sample_id + 1),
        })
    }

    /// Insert a sample and return its ID.
    pub fn insert_sample(&self, sample_name: &str, sequencing_type: &str, total_reads: u64, mapped_reads: u64, circular_mapping: bool) -> Result<i64> {
        // Get next sample ID
        let sample_id = {
            let mut next_id = self.next_sample_id.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
            let id = *next_id;
            *next_id += 1;
            id
        };

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        conn.execute(
            "INSERT INTO Sample (Sample_id, Sample_name, Sequencing_type, Number_of_reads, Number_of_mapped_reads, Circular_mapping) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![sample_id, sample_name, sequencing_type, total_reads as i64, mapped_reads as i64, circular_mapping],
        )?;
        drop(conn);

        let mut sample_map = self.sample_name_to_id.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        sample_map.insert(sample_name.to_string(), sample_id);

        Ok(sample_id)
    }

    /// Write all data for a sample (coverage, packaging, misassembly/microdiversity/side_misassembly/topology, feature blobs).
    pub fn write_sample_data(
        &self,
        sample_name: &str,
        presences: &[PresenceData],
        packaging: &[PackagingData],
        misassembly: &[MisassemblyData],
        microdiversity: &[MicrodiversityData],
        side_misassembly: &[SideMisassemblyData],
        topology: &[TopologyData],
        feature_blobs: &[(String, String, crate::blob::EncodedBlob)],
        circular: bool,
    ) -> Result<()> {
        let sample_id = {
            let sample_map = self.sample_name_to_id.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
            *sample_map.get(sample_name).context("Sample not found")?
        };

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

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

        // Insert feature BLOBs (compressed format)
        self.write_feature_blobs(&conn, sample_id, feature_blobs)?;

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
                let cov_sd = p.coverage_coefficient_of_variation.round() as i32;
                let cov_var = p.coverage_relative_coverage_roughness.round() as i32;
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
            "SELECT COALESCE(MAX(Packaging_id), 0) FROM Phage_mechanisms",
            [],
            |row| row.get(0),
        ).unwrap_or(0) + 1;

        let mut terminus_id: i64 = conn.query_row(
            "SELECT COALESCE(MAX(Terminus_id), 0) FROM Phage_termini",
            [],
            |row| row.get(0),
        ).unwrap_or(0) + 1;

        let mut mechanism_appender = conn.appender("Phage_mechanisms")
            .context("Failed to create Phage_mechanisms appender")?;
        let mut termini_appender = conn.appender("Phage_termini")
            .context("Failed to create Phage_termini appender")?;

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

                // Insert into Phage_mechanisms
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

        mechanism_appender.flush().context("Failed to flush Phage_mechanisms appender")?;
        termini_appender.flush().context("Failed to flush Phage_termini appender")?;
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
                let median_circ_len = d.median_circularising_len;
                let circ_inserts = d.circularising_inserts.map(|v| v as i64);
                let circ_dev = d.circularising_insert_size_deviation;
                appender.append_row(params![
                    contig_id, sample_id,
                    circ_reads, circ_pct, median_circ_len, circ_inserts, circ_dev
                ])?;
            }
        }

        appender.flush().context("Failed to flush Topology appender")?;
        Ok(())
    }

    /// Write feature BLOBs to the Feature_blob table.
    ///
    /// Each tuple is (feature_name, contig_name, blob_bytes).
    /// Feature names are mapped to feature_ids via VARIABLES index.
    pub fn write_feature_blobs(&self, conn: &Connection, sample_id: i64, blobs: &[(String, String, crate::blob::EncodedBlob)]) -> Result<()> {
        if blobs.is_empty() {
            return Ok(());
        }

        let mut appender = conn.appender("Feature_blob")
            .context("Failed to create Feature_blob appender")?;
        let mut chunk_appender = conn.appender("Feature_blob_chunk")
            .context("Failed to create Feature_blob_chunk appender")?;

        for (feature_name, contig_name, encoded) in blobs {
            // Skip empty blobs (all-zero dense or no-event sparse)
            if encoded.zoom.is_empty() {
                continue;
            }
            if let (Some(&contig_id), Some(fid)) = (
                self.contig_name_to_id.get(contig_name),
                feature_name_to_id(feature_name),
            ) {
                appender.append_row(params![contig_id, sample_id, fid as i32, encoded.zoom.as_slice()])?;

                // Write individual base-resolution chunks for fast zoomed-in queries
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    chunk_appender.append_row(params![
                        contig_id, sample_id, fid as i32, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Feature_blob appender")?;
        chunk_appender.flush().context("Failed to flush Feature_blob_chunk appender")?;
        Ok(())
    }

    /// Write repeats data to the database (both direct and inverted).
    /// This is called once during database creation (not per-sample).
    pub fn write_repeats(&self, repeats: &[RepeatsData]) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

        if repeats.is_empty() {
            // Autoblast ran but found no repeats — set 0% for all contigs
            conn.execute("UPDATE Contig SET Duplication_percentage = 0", [])?;
            return Ok(());
        }
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

    /// Write GC content data to Contig_blob table.
    /// Stores raw window-level GC percentages as dense BLOBs (one value per 500bp window).
    /// This is called once during database creation (not per-sample).
    pub fn write_gc_content(&self, gc_data: &[GCContentData]) -> Result<()> {
        use crate::blob::{encode_contig_dense_blob, ValueScale};

        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        const GC_CONTENT_FEATURE_ID: i16 = 1;
        const GC_WINDOW_SIZE: u32 = 500; // Standard GC content window size

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                if data.gc_values.is_empty() {
                    continue;
                }

                // Convert u8 values to i32 for BLOB encoding
                let values: Vec<i32> = data.gc_values.iter().map(|&v| v as i32).collect();

                // Estimate contig length from number of windows
                let contig_length = (values.len() as u32) * GC_WINDOW_SIZE;

                // Encode as dense BLOB (array index is window number)
                let encoded = encode_contig_dense_blob(
                    &values,
                    ValueScale::Raw,
                    contig_length,
                    GC_WINDOW_SIZE,
                );

                appender.append_row(params![
                    contig_id,
                    GC_CONTENT_FEATURE_ID,
                    encoded.zoom
                ])?;
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    chunk_appender.append_row(params![
                        contig_id, GC_CONTENT_FEATURE_ID, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;

        Ok(())
    }

    /// Write GC skew data to Contig_blob table.
    /// Stores raw window-level GC skew values as dense BLOBs (one value per 1000bp window).
    /// This is called once during database creation (not per-sample).
    pub fn write_gc_skew(&self, gc_data: &[GCContentData]) -> Result<()> {
        use crate::blob::{encode_contig_dense_blob, ValueScale};

        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        const GC_SKEW_FEATURE_ID: i16 = 2;
        const GC_SKEW_WINDOW_SIZE: u32 = 1000; // Standard GC skew window size

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                if data.skew_values.is_empty() {
                    continue;
                }

                // Convert i16 values to i32 for BLOB encoding
                let values: Vec<i32> = data.skew_values.iter().map(|&v| v as i32).collect();

                // Estimate contig length from number of windows
                let contig_length = (values.len() as u32) * GC_SKEW_WINDOW_SIZE;

                // Encode as dense BLOB (array index is window number)
                // Values are stored as integer ×100 (e.g. -40 = -0.40), use Times100 scale
                let encoded = encode_contig_dense_blob(
                    &values,
                    ValueScale::Times100,
                    contig_length,
                    GC_SKEW_WINDOW_SIZE,
                );

                appender.append_row(params![
                    contig_id,
                    GC_SKEW_FEATURE_ID,
                    encoded.zoom
                ])?;
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    chunk_appender.append_row(params![
                        contig_id, GC_SKEW_FEATURE_ID, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;

        Ok(())
    }

    /// Write repeat count and identity features to Contig_blob table from materialized tables.
    /// Converts RLE segments into dense per-bp arrays, then encodes as dense contig BLOBs.
    /// Delta+zstd compression handles constant-value runs very efficiently.
    fn write_repeat_features_to_blob_static(conn: &Connection) -> Result<()> {
        use crate::blob::{encode_dense_blob, ValueScale};

        // Feature IDs for repeat features
        const DIRECT_REPEAT_COUNT_ID: i16 = 3;
        const INVERTED_REPEAT_COUNT_ID: i16 = 4;
        const DIRECT_REPEAT_IDENTITY_ID: i16 = 5;
        const INVERTED_REPEAT_IDENTITY_ID: i16 = 6;

        // List of (table_name, feature_id) to process
        let tables = vec![
            ("Contig_direct_repeat_count", DIRECT_REPEAT_COUNT_ID),
            ("Contig_inverted_repeat_count", INVERTED_REPEAT_COUNT_ID),
            ("Contig_direct_repeat_identity", DIRECT_REPEAT_IDENTITY_ID),
            ("Contig_inverted_repeat_identity", INVERTED_REPEAT_IDENTITY_ID),
        ];

        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        for (table_name, feature_id) in tables {
            // Check if table exists (it might not if there are no repeats)
            let table_exists: bool = conn.query_row(
                "SELECT COUNT(*) > 0 FROM information_schema.tables WHERE table_name = ?",
                params![table_name],
                |row| row.get(0),
            ).unwrap_or(false);

            if !table_exists {
                eprintln!("Skipping {} (table not found)", table_name);
                continue;
            }

            // Identity features store 0.0-1.0 values → scale ×100 for integer storage
            let is_identity = feature_id == DIRECT_REPEAT_IDENTITY_ID
                           || feature_id == INVERTED_REPEAT_IDENTITY_ID;
            let scale = if is_identity { ValueScale::Times100 } else { ValueScale::Raw };

            // Get contig length map
            let mut contig_lengths: HashMap<i64, u32> = HashMap::new();
            {
                let mut stmt = conn.prepare("SELECT Contig_id, Contig_length FROM Contig")?;
                let rows = stmt.query_map([], |row| {
                    Ok((row.get::<_, i64>(0)?, row.get::<_, i32>(1)?))
                })?;
                for row in rows {
                    let (contig_id, length) = row?;
                    contig_lengths.insert(contig_id, length as u32);
                }
            }

            // Read RLE segments from the materialized table
            let query = format!("SELECT Contig_id, First_position, Last_position, Value FROM {} ORDER BY Contig_id, First_position", table_name);
            let mut stmt = conn.prepare(&query)?;
            let rows = stmt.query_map([], |row| {
                Ok((
                    row.get::<_, i64>(0)?,    // Contig_id
                    row.get::<_, i32>(1)?,    // First_position
                    row.get::<_, i32>(2)?,    // Last_position
                    row.get::<_, f64>(3)?,    // Value
                ))
            })?;

            // Collect all segments grouped by contig, then encode each as dense blob
            let mut contig_segments: HashMap<i64, Vec<(i32, i32, f64)>> = HashMap::new();

            for row in rows {
                let (contig_id, first_pos, last_pos, value) = row?;
                contig_segments.entry(contig_id).or_default().push((first_pos, last_pos, value));
            }

            for (contig_id, segments) in &contig_segments {
                let contig_length = *contig_lengths.get(contig_id).unwrap_or(&0);
                if contig_length == 0 || segments.is_empty() {
                    continue;
                }
                // Build dense per-bp array (zeros, then fill RLE runs)
                let mut dense = vec![0i32; contig_length as usize];
                for &(first_pos, last_pos, value) in segments {
                    let v = if is_identity { (value * 100.0) as i32 } else { value as i32 };
                    let start = (first_pos as usize).min(dense.len());
                    let end = ((last_pos + 1) as usize).min(dense.len());
                    for pos in start..end {
                        dense[pos] = v;
                    }
                }
                let encoded = encode_dense_blob(&dense, scale, contig_length);
                appender.append_row(params![*contig_id, feature_id, encoded.zoom])?;
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    chunk_appender.append_row(params![
                        *contig_id, feature_id, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;
        Ok(())
    }

    /// Update Contig table with GC statistics (average, sd) and GC skew stats.
    /// Called after computing GC content and GC skew for all contigs.
    pub fn update_contig_gc_stats(&self, gc_data: &[GCContentData]) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
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
    /// Write initial metadata rows into the Database_metadata table.
    pub fn write_metadata(
        &self,
        modules: &[String],
        min_aligned_fraction: f64,
        min_coverage_depth: f64,
        curve_ratio: f64,
        bar_ratio: f64,
    ) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let now = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

        // Build version string with optional git short hash
        let mut version = env!("CARGO_PKG_VERSION").to_string();
        if let Ok(output) = std::process::Command::new("git")
            .args(["rev-parse", "--short", "HEAD"])
            .output()
        {
            let hash = String::from_utf8_lossy(&output.stdout).trim().to_string();
            if !hash.is_empty() {
                version = format!("{}+{}", version, hash);
            }
        }

        let modules_str = modules.join(",");

        let rows: Vec<(&str, String)> = vec![
            ("Date_of_creation", now.clone()),
            ("Date_of_last_modification", now),
            ("Tool_version_used_for_creation", version.clone()),
            ("Tool_version_used_for_last_modification", version),
            ("Modules", modules_str),
            ("Min_aligned_fraction", min_aligned_fraction.to_string()),
            ("Min_coverage_depth", min_coverage_depth.to_string()),
            ("Variation_percentage", curve_ratio.to_string()),
            ("Coverage_percentage", bar_ratio.to_string()),
        ];

        for (key, value) in &rows {
            conn.execute(
                "INSERT INTO Database_metadata (Key, Value) VALUES (?, ?)",
                params![key, value],
            )
            .with_context(|| format!("Failed to insert metadata key '{}'", key))?;
        }

        Ok(())
    }

    /// Update modification-related metadata (for extend mode).
    pub fn update_metadata_modification(&self) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let now = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

        let mut version = env!("CARGO_PKG_VERSION").to_string();
        if let Ok(output) = std::process::Command::new("git")
            .args(["rev-parse", "--short", "HEAD"])
            .output()
        {
            let hash = String::from_utf8_lossy(&output.stdout).trim().to_string();
            if !hash.is_empty() {
                version = format!("{}+{}", version, hash);
            }
        }

        conn.execute(
            "UPDATE Database_metadata SET Value = ? WHERE Key = 'Date_of_last_modification'",
            params![now],
        )
        .context("Failed to update Date_of_last_modification")?;

        conn.execute(
            "UPDATE Database_metadata SET Value = ? WHERE Key = 'Tool_version_used_for_last_modification'",
            params![version],
        )
        .context("Failed to update Tool_version_used_for_last_modification")?;

        Ok(())
    }

    pub fn finalize(self) -> Result<()> {
        let conn = self.conn.into_inner().unwrap();

        // Create derived views (in extend mode, views were dropped in open())
        create_views(&conn, self.has_bam)?;

        // Delete Variable entries for features without data
        cleanup_unused_variables(&conn)?;

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
        "SELECT ca.\"Start\", ca.\"End\"
         FROM Contig_annotation_core ca
         JOIN Annotation_qualifier aq ON aq.Annotation_id = ca.Annotation_id
         WHERE ca.Contig_id = ? AND aq.\"Key\" = 'product' AND aq.\"Value\" LIKE '%terminase%'"
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
            Contig_id INTEGER REFERENCES Contig(Contig_id),
            Sample_id INTEGER REFERENCES Sample(Sample_id),
            Aligned_fraction_percentage INTEGER,
            Above_expected_aligned_fraction BOOLEAN,
            Read_count INTEGER,
            Coverage_mean INTEGER,
            Coverage_median INTEGER,
            Coverage_trimmed_mean INTEGER,
            Coverage_coefficient_of_variation INTEGER,
            Coverage_relative_coverage_roughness INTEGER,
            PRIMARY KEY (Contig_id, Sample_id)
        )",
        [],
    )
    .context("Failed to create Coverage table")?;

    conn.execute(
        "CREATE TABLE Phage_mechanisms (
            Packaging_id INTEGER PRIMARY KEY,
            Contig_id INTEGER NOT NULL REFERENCES Contig(Contig_id),
            Sample_id INTEGER NOT NULL REFERENCES Sample(Sample_id),
            Packaging_mechanism TEXT NOT NULL,
            Duplication BOOLEAN,
            Repeat_length INTEGER,
            Terminase_distance INTEGER,
            Median_left_termini_clippings TEXT,
            Median_right_termini_clippings TEXT
        )",
        [],
    )
    .context("Failed to create Phage_mechanisms table")?;

    // Phage_termini table - stores individual terminus areas with full metadata
    // One row per terminus area, linked to Phage_mechanisms via Packaging_id
    // Includes filtering diagnostics for both Poisson and clipping tests
    conn.execute(
        "CREATE TABLE Phage_termini (
            Terminus_id INTEGER PRIMARY KEY,
            Packaging_id INTEGER NOT NULL REFERENCES Phage_mechanisms(Packaging_id),
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
    .context("Failed to create Phage_termini table")?;
    // Feature_blob table - zoom-level summaries for per-position feature data.
    // Base-resolution data is stored in Feature_blob_chunk (one row per 65kbp chunk).
    conn.execute(
        "CREATE TABLE Feature_blob (
            Contig_id  INTEGER NOT NULL REFERENCES Contig(Contig_id),
            Sample_id  INTEGER NOT NULL REFERENCES Sample(Sample_id),
            Feature_id SMALLINT NOT NULL,
            Zoom_data  BLOB NOT NULL,
            PRIMARY KEY (Contig_id, Sample_id, Feature_id)
        )",
        [],
    )
    .context("Failed to create Feature_blob table")?;

    // Feature_blob_chunk table - individual base-resolution chunks (65,536 bp each).
    // Python fetches only the 1-2 chunks overlapping a zoomed-in view instead of
    // the entire multi-MB Data BLOB from Feature_blob.
    conn.execute(
        "CREATE TABLE Feature_blob_chunk (
            Contig_id  INTEGER NOT NULL,
            Sample_id  INTEGER NOT NULL,
            Feature_id SMALLINT NOT NULL,
            Chunk_idx  SMALLINT NOT NULL,
            Data       BLOB NOT NULL,
            PRIMARY KEY (Contig_id, Sample_id, Feature_id, Chunk_idx)
        )",
        [],
    )
    .context("Failed to create Feature_blob_chunk table")?;

    } // end if has_bam (Sample, Coverage, Phage_mechanisms, Phage_termini, Feature_blob, Feature_blob_chunk)

    // Contig_blob table - zoom-level summaries for contig-level features (GC content, GC skew, repeats).
    // Base-resolution data is stored in Contig_blob_chunk.
    conn.execute(
        "CREATE TABLE Contig_blob (
            Contig_id   INTEGER NOT NULL REFERENCES Contig,
            Feature_id  SMALLINT NOT NULL,
            Zoom_data   BLOB NOT NULL,
            PRIMARY KEY (Contig_id, Feature_id)
        )",
        [],
    )
    .context("Failed to create Contig_blob table")?;

    // Contig_blob_chunk table - individual base-resolution chunks for contig-level features
    conn.execute(
        "CREATE TABLE Contig_blob_chunk (
            Contig_id   INTEGER NOT NULL,
            Feature_id  SMALLINT NOT NULL,
            Chunk_idx   SMALLINT NOT NULL,
            Data        BLOB NOT NULL,
            PRIMARY KEY (Contig_id, Feature_id, Chunk_idx)
        )",
        [],
    )
    .context("Failed to create Contig_blob_chunk table")?;

    // Column_scales table - documents scaling factors for all features/columns
    // Both Rust (encoding) and Python (decoding) can read from here instead of hardcoding
    conn.execute(
        "CREATE TABLE Column_scales (
            Feature_name TEXT NOT NULL,
            Column_name TEXT NOT NULL,
            Scale INTEGER NOT NULL,
            PRIMARY KEY (Feature_name, Column_name)
        )",
        [],
    )
    .context("Failed to create Column_scales table")?;

    // Populate default scaling factors
    let scale_entries: &[(&str, &str, i32)] = &[
        // Coverage table columns (×100 for integer storage)
        ("Coverage", "Coverage_mean", 100),
        ("Coverage", "Coverage_median", 100),
        ("Coverage", "Coverage_trimmed_mean", 100),
        ("Coverage", "Coverage_coefficient_of_variation", 1_000_000),
        ("Coverage", "Coverage_relative_coverage_roughness", 1_000_000),
        ("Coverage", "Aligned_fraction_percentage", 100),
        // MapQ (×100 for decimal precision)
        ("mapq", "Value", 100),
        // Coverage-relative features (×1000 for per-mille)
        ("mismatches", "Value_relative", 1000),
        ("deletions", "Value_relative", 1000),
        ("insertions", "Value_relative", 1000),
        ("left_clippings", "Value_relative", 1000),
        ("right_clippings", "Value_relative", 1000),
        ("non_inward_pairs", "Value_relative", 1000),
        ("mate_not_mapped", "Value_relative", 1000),
        ("mate_on_another_contig", "Value_relative", 1000),
        ("reads_starts", "Value_relative", 1000),
        ("reads_ends", "Value_relative", 1000),
        // GC content/skew (×1000 for 3 decimal places)
        ("gc_content", "Value", 1000),
        ("gc_skew", "Value", 1000),
        // Contig-level GC stats (×100)
        ("Contig", "GC_mean", 100),
        ("Contig", "GC_sd", 100),
        ("Contig", "GC_skew_amplitude", 100),
        // Repeat identity (×100)
        ("Contig_directRepeats", "Pident", 100),
        ("Contig_invertedRepeats", "Pident", 100),
    ];

    {
        let mut appender = conn.appender("Column_scales")
            .context("Failed to create Column_scales appender")?;
        for (feature, column, scale) in scale_entries {
            appender.append_row(params![*feature, *column, *scale])
                .context("Failed to append Column_scales row")?;
        }
        appender.flush().context("Failed to flush Column_scales appender")?;
    }

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

    // Note: Contig_GCContent and Contig_GCSkew tables have been removed.
    // GC content and GC skew are now stored in Contig_blob table.

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
            Median_circularising_len INTEGER,
            Circularising_inserts INTEGER,
            Circularising_insert_size_deviation INTEGER
        )",
        [],
    )
    .context("Failed to create Topology table")?;
    } // end if has_bam (Misassembly, Microdiversity, Side_misassembly, Topology)

    // Structural annotation data only — qualifiers go in Annotation_qualifier.
    // Main_isoform flags the canonical transcript per (locus_tag, feature_type):
    //   - The first mRNA seen in file order per locus_tag is the main mRNA.
    //   - CDS/exon rows whose Parent_annotation_id points at the main mRNA
    //     inherit Main_isoform = TRUE.
    //   - For features without an mRNA parent (prokaryotic CDS, tRNA, ...),
    //     the first occurrence per (locus_tag, feature_type) is main.
    // Parent_annotation_id links CDS/exon children to their parent mRNA
    // (populated from GFF3 Parent= or GenBank segment containment matching).
    conn.execute(
        "CREATE TABLE Contig_annotation_core (
            Annotation_id INTEGER PRIMARY KEY,
            Contig_id INTEGER REFERENCES Contig(Contig_id),
            \"Start\" INTEGER,
            \"End\" INTEGER,
            Strand INTEGER,
            \"Type\" TEXT,
            Main_isoform BOOLEAN,
            Parent_annotation_id INTEGER
        )",
        [],
    )
    .context("Failed to create Contig_annotation_core table")?;

    // Sub-intervals of spliced CDS / mRNA features. Only populated when a
    // feature has two or more segments (GenBank join(...), GFF3 exon/CDS
    // children collapsed into their parent mRNA). Unspliced features have
    // no rows here — their bounding box in Contig_annotation_core is
    // sufficient.
    conn.execute(
        "CREATE TABLE Annotation_segments (
            Annotation_id INTEGER,
            Segment_index INTEGER,
            Start_segment INTEGER,
            End_segment INTEGER
        )",
        [],
    )
    .context("Failed to create Annotation_segments table")?;

    conn.execute(
        "CREATE INDEX idx_annotation_segments_id ON Annotation_segments(Annotation_id)",
        [],
    )
    .context("Failed to create Annotation_segments index")?;

    // Key-value store for all annotation qualifiers (product, function, locus_tag, phrog, gene, note, db_xref, ...)
    conn.execute(
        "CREATE TABLE Annotation_qualifier (
            Annotation_id INTEGER,
            \"Key\" TEXT,
            \"Value\" TEXT
        )",
        [],
    )
    .context("Failed to create Annotation_qualifier table")?;

    // Key-value store for contig-level qualifiers from GenBank "source" features (organism, host, strain, ...)
    conn.execute(
        "CREATE TABLE Contig_qualifier (
            Contig_id INTEGER,
            \"Key\" TEXT,
            \"Value\" TEXT
        )",
        [],
    )
    .context("Failed to create Contig_qualifier table")?;

    // Heavy sequence data stored separately — only queried for codon analysis and export
    conn.execute(
        "CREATE TABLE Annotation_sequence (
            Annotation_id INTEGER PRIMARY KEY,
            Nucleotide_sequence TEXT,
            Protein_sequence TEXT,
            S_sites INTEGER,
            N_sites INTEGER
        )",
        [],
    )
    .context("Failed to create Annotation_sequence table")?;

    // Convenience view that pivots common qualifiers into columns for backward compatibility.
    // All existing Python queries use this view.
    // `Segments` is a list of {start, end} structs in genomic order when the
    // feature is spliced, NULL otherwise.
    conn.execute(
        "CREATE VIEW Contig_annotation AS
         SELECT ca.Annotation_id, ca.Contig_id, ca.\"Start\", ca.\"End\", ca.Strand, ca.\"Type\",
            ca.Main_isoform, ca.Parent_annotation_id,
            aseq.Nucleotide_sequence,
            aseq.Protein_sequence,
            aseq.S_sites,
            aseq.N_sites,
            seg.Segments
         FROM Contig_annotation_core ca
         LEFT JOIN Annotation_sequence aseq ON ca.Annotation_id = aseq.Annotation_id
         LEFT JOIN (
            SELECT Annotation_id,
                   list({'start': Start_segment, 'end': End_segment} ORDER BY Segment_index) AS Segments
            FROM Annotation_segments
            GROUP BY Annotation_id
         ) seg ON ca.Annotation_id = seg.Annotation_id",
        [],
    )
    .context("Failed to create Contig_annotation view")?;

    conn.execute(
        "CREATE TABLE Variable (
            Variable_id INTEGER PRIMARY KEY,
            Variable_name TEXT UNIQUE,
            Subplot TEXT,
            Module TEXT,
            Module_order INTEGER,
            \"Type\" TEXT,
            Encoding TEXT,
            Value_scale TEXT,
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

    // Constants_for_plotting table - stores metadata flags about database content
    // Pre-computed at database creation time (not per-request)
    conn.execute(
        "CREATE TABLE Constants_for_plotting (
            Constant TEXT PRIMARY KEY,
            Status BOOLEAN
        )",
        [],
    )
    .context("Failed to create Constants_for_plotting table")?;

    // Color_templates table - named color schemes (e.g. pharokka, generic)
    conn.execute(
        "CREATE TABLE Color_templates (
            Template_id INTEGER PRIMARY KEY,
            Template_name TEXT UNIQUE NOT NULL
        )",
        [],
    )
    .context("Failed to create Color_templates table")?;

    // Color_rules table - individual coloring rules belonging to a template
    conn.execute(
        "CREATE TABLE Color_rules (
            Rule_id INTEGER PRIMARY KEY,
            Template_id INTEGER NOT NULL REFERENCES Color_templates(Template_id),
            Qualifier_name TEXT NOT NULL,
            Operator TEXT NOT NULL,
            Qualifier_value TEXT NOT NULL,
            Color TEXT NOT NULL
        )",
        [],
    )
    .context("Failed to create Color_rules table")?;

    // Database_metadata table - tracks how the database was created/modified
    conn.execute(
        "CREATE TABLE Database_metadata (
            Key TEXT PRIMARY KEY,
            Value TEXT
        )",
        [],
    )
    .context("Failed to create Database_metadata table")?;

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
            AminoAcid_label TEXT NOT NULL,
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

    // Standard genetic code: (codon, 1-letter AA, full name, 3-letter label, color)
    // Colors by biochemical property:
    //   Hydrophobic aliphatic (G,A,V,L,I,P): #2d6a4f
    //   Aromatic (F,W,Y): #b5838d
    //   Polar uncharged (S,T,N,Q,C,M): #457b9d
    //   Positively charged (K,R,H): #9b2226
    //   Negatively charged (D,E): #e9c46a
    //   Stop (*): #6a3d9a
    let codons: &[(&str, &str, &str, &str, &str)] = &[
        ("TTT", "F", "Phenylalanine", "Phe", "#b5838d"), ("TTC", "F", "Phenylalanine", "Phe", "#b5838d"),
        ("TTA", "L", "Leucine", "Leu", "#2d6a4f"), ("TTG", "L", "Leucine", "Leu", "#2d6a4f"),
        ("TCT", "S", "Serine", "Ser", "#457b9d"), ("TCC", "S", "Serine", "Ser", "#457b9d"),
        ("TCA", "S", "Serine", "Ser", "#457b9d"), ("TCG", "S", "Serine", "Ser", "#457b9d"),
        ("TAT", "Y", "Tyrosine", "Tyr", "#b5838d"), ("TAC", "Y", "Tyrosine", "Tyr", "#b5838d"),
        ("TAA", "*", "Stop", "*", "#6a3d9a"), ("TAG", "*", "Stop", "*", "#6a3d9a"),
        ("TGT", "C", "Cysteine", "Cys", "#457b9d"), ("TGC", "C", "Cysteine", "Cys", "#457b9d"),
        ("TGA", "*", "Stop", "*", "#6a3d9a"), ("TGG", "W", "Tryptophan", "Trp", "#b5838d"),
        ("CTT", "L", "Leucine", "Leu", "#2d6a4f"), ("CTC", "L", "Leucine", "Leu", "#2d6a4f"),
        ("CTA", "L", "Leucine", "Leu", "#2d6a4f"), ("CTG", "L", "Leucine", "Leu", "#2d6a4f"),
        ("CCT", "P", "Proline", "Pro", "#2d6a4f"), ("CCC", "P", "Proline", "Pro", "#2d6a4f"),
        ("CCA", "P", "Proline", "Pro", "#2d6a4f"), ("CCG", "P", "Proline", "Pro", "#2d6a4f"),
        ("CAT", "H", "Histidine", "His", "#9b2226"), ("CAC", "H", "Histidine", "His", "#9b2226"),
        ("CAA", "Q", "Glutamine", "Gln", "#457b9d"), ("CAG", "Q", "Glutamine", "Gln", "#457b9d"),
        ("CGT", "R", "Arginine", "Arg", "#9b2226"), ("CGC", "R", "Arginine", "Arg", "#9b2226"),
        ("CGA", "R", "Arginine", "Arg", "#9b2226"), ("CGG", "R", "Arginine", "Arg", "#9b2226"),
        ("ATT", "I", "Isoleucine", "Ile", "#2d6a4f"), ("ATC", "I", "Isoleucine", "Ile", "#2d6a4f"),
        ("ATA", "I", "Isoleucine", "Ile", "#2d6a4f"), ("ATG", "M", "Methionine", "Met", "#457b9d"),
        ("ACT", "T", "Threonine", "Thr", "#457b9d"), ("ACC", "T", "Threonine", "Thr", "#457b9d"),
        ("ACA", "T", "Threonine", "Thr", "#457b9d"), ("ACG", "T", "Threonine", "Thr", "#457b9d"),
        ("AAT", "N", "Asparagine", "Asn", "#457b9d"), ("AAC", "N", "Asparagine", "Asn", "#457b9d"),
        ("AAA", "K", "Lysine", "Lys", "#9b2226"), ("AAG", "K", "Lysine", "Lys", "#9b2226"),
        ("AGT", "S", "Serine", "Ser", "#457b9d"), ("AGC", "S", "Serine", "Ser", "#457b9d"),
        ("AGA", "R", "Arginine", "Arg", "#9b2226"), ("AGG", "R", "Arginine", "Arg", "#9b2226"),
        ("GTT", "V", "Valine", "Val", "#2d6a4f"), ("GTC", "V", "Valine", "Val", "#2d6a4f"),
        ("GTA", "V", "Valine", "Val", "#2d6a4f"), ("GTG", "V", "Valine", "Val", "#2d6a4f"),
        ("GCT", "A", "Alanine", "Ala", "#2d6a4f"), ("GCC", "A", "Alanine", "Ala", "#2d6a4f"),
        ("GCA", "A", "Alanine", "Ala", "#2d6a4f"), ("GCG", "A", "Alanine", "Ala", "#2d6a4f"),
        ("GAT", "D", "Aspartate", "Asp", "#e9c46a"), ("GAC", "D", "Aspartate", "Asp", "#e9c46a"),
        ("GAA", "E", "Glutamate", "Glu", "#e9c46a"), ("GAG", "E", "Glutamate", "Glu", "#e9c46a"),
        ("GGT", "G", "Glycine", "Gly", "#2d6a4f"), ("GGC", "G", "Glycine", "Gly", "#2d6a4f"),
        ("GGA", "G", "Glycine", "Gly", "#2d6a4f"), ("GGG", "G", "Glycine", "Gly", "#2d6a4f"),
    ];

    for &(codon, aa, name, label, color) in codons {
        appender.append_row(params![codon, aa, name, label, color])?;
    }

    appender.flush().context("Failed to flush Codon_table appender")?;
    Ok(())
}

/// Compute synonymous and nonsynonymous site counts for a nucleotide sequence.
///
/// For each codon, enumerates all 9 possible single-nucleotide changes and classifies
/// each as synonymous (same amino acid) or nonsynonymous (different amino acid).
/// Codons with non-ACGT bases or stop codons are skipped.
/// Returns (S_sites, N_sites) as integers.
fn compute_sn_sites(nuc_seq: &str) -> (i64, i64) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let seq = nuc_seq.as_bytes();
    let mut syn: i64 = 0;
    let mut nonsyn: i64 = 0;

    for chunk in seq.chunks_exact(3) {
        let codon: [u8; 3] = [
            chunk[0].to_ascii_uppercase(),
            chunk[1].to_ascii_uppercase(),
            chunk[2].to_ascii_uppercase(),
        ];

        // Skip codons with non-ACGT bases
        if !codon.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }

        // Translate the original codon; skip stop codons
        let (orig_aa, _) = match translate_codon(&codon) {
            Some(v) => v,
            None => continue,
        };
        if orig_aa == '*' {
            continue;
        }

        // Enumerate all 9 single-nucleotide changes
        for pos in 0..3 {
            for &alt in &BASES {
                if alt == codon[pos] {
                    continue;
                }
                let mut mutant = codon;
                mutant[pos] = alt;
                if let Some((mut_aa, _)) = translate_codon(&mutant) {
                    if mut_aa == orig_aa {
                        syn += 1;
                    } else {
                        nonsyn += 1;
                    }
                }
            }
        }
    }

    (syn, nonsyn)
}

/// Insert annotations into the database using Appender for bulk loading.
///
/// Implements the three-pass Main_isoform / Parent_annotation_id logic:
///   A. Assign sequential Annotation_ids and build self_key → id map.
///   B. Resolve parent_key → parent_annotation_id.
///   C. Compute Main_isoform via first-occurrence rule:
///      - The first mRNA seen per locus_tag is the main mRNA.
///      - CDS/exon rows whose parent is the main mRNA inherit Main_isoform = TRUE.
///      - For rows without a parent, the first per (locus_tag, feature_type) is main.
///
/// Also populates Annotation_segments for any feature with >= 2 segments.
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

    let mut has_pharokka_function = false;
    for ann in annotations {
        if !has_pharokka_function && ann.feature_type == "CDS" {
            if let Some(function) = ann.qualifiers.get("function") {
                let func_lower = function.to_lowercase();
                for pharokka_key in &pharokka_keys {
                    if func_lower.contains(pharokka_key) {
                        has_pharokka_function = true;
                        break;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Pass A: assign sequential ids and build self_key → id map.
    // ------------------------------------------------------------------
    let ids: Vec<i64> = (1..=annotations.len() as i64).collect();
    let mut self_key_to_id: std::collections::HashMap<String, i64> =
        std::collections::HashMap::new();
    for (ann, &id) in annotations.iter().zip(ids.iter()) {
        if let Some(sk) = &ann.self_key {
            self_key_to_id.entry(sk.clone()).or_insert(id);
        }
    }

    // ------------------------------------------------------------------
    // Pass B: resolve parent_key → parent_annotation_id.
    // ------------------------------------------------------------------
    let parent_ids: Vec<Option<i64>> = annotations
        .iter()
        .map(|ann| {
            ann.parent_key
                .as_ref()
                .and_then(|pk| self_key_to_id.get(pk).copied())
        })
        .collect();

    // ------------------------------------------------------------------
    // Pass C: first-occurrence Main_isoform selection.
    // ------------------------------------------------------------------
    // Step 1: first mRNA per locus_tag → main_mrna_ids.
    let mut main_mrna_ids: std::collections::HashSet<i64> = std::collections::HashSet::new();
    let mut seen_mrna_locus: std::collections::HashSet<String> =
        std::collections::HashSet::new();
    for (ann, &id) in annotations.iter().zip(ids.iter()) {
        if ann.feature_type != "mRNA" {
            continue;
        }
        if let Some(lt) = ann.qualifiers.get("locus_tag") {
            if seen_mrna_locus.insert(lt.clone()) {
                main_mrna_ids.insert(id);
            }
        }
    }

    // Step 2: first flat row per (locus_tag, feature_type) for rows WITHOUT a
    // resolved parent link — used only when parent_annotation_id is None.
    let mut main_flat_ids: std::collections::HashSet<i64> = std::collections::HashSet::new();
    let mut seen_flat: std::collections::HashSet<(String, String)> =
        std::collections::HashSet::new();
    for (i, (ann, &id)) in annotations.iter().zip(ids.iter()).enumerate() {
        if parent_ids[i].is_some() || ann.feature_type == "mRNA" {
            continue;
        }
        if let Some(lt) = ann.qualifiers.get("locus_tag") {
            let key = (lt.clone(), ann.feature_type.clone());
            if seen_flat.insert(key) {
                main_flat_ids.insert(id);
            }
        }
    }

    // Step 3: per-annotation Main_isoform resolution.
    let main_isoforms: Vec<Option<bool>> = annotations
        .iter()
        .enumerate()
        .map(|(i, ann)| {
            if ann.qualifiers.get("locus_tag").is_none() {
                return None;
            }
            let id = ids[i];
            if ann.feature_type == "mRNA" {
                return Some(main_mrna_ids.contains(&id));
            }
            if let Some(pid) = parent_ids[i] {
                return Some(main_mrna_ids.contains(&pid));
            }
            Some(main_flat_ids.contains(&id))
        })
        .collect();

    // ------------------------------------------------------------------
    // Insertion.
    // ------------------------------------------------------------------
    let mut ann_appender = conn
        .appender("Contig_annotation_core")
        .context("Failed to create Contig_annotation_core appender")?;
    let mut seq_appender = conn
        .appender("Annotation_sequence")
        .context("Failed to create Annotation_sequence appender")?;
    let mut qual_appender = conn
        .appender("Annotation_qualifier")
        .context("Failed to create Annotation_qualifier appender")?;
    let mut seg_appender = conn
        .appender("Annotation_segments")
        .context("Failed to create Annotation_segments appender")?;

    for (i, ann) in annotations.iter().enumerate() {
        let annotation_id = ids[i];
        let main_isoform = main_isoforms[i];
        let parent_annotation_id = parent_ids[i];

        let (s_sites, n_sites): (Option<i64>, Option<i64>) = match &ann.nucleotide_sequence {
            Some(seq) if !seq.is_empty() => {
                let (s, n) = compute_sn_sites(seq);
                (Some(s), Some(n))
            }
            _ => (None, None),
        };

        // Structural annotation row
        ann_appender
            .append_row(params![
                annotation_id,
                ann.contig_id,
                ann.start,
                ann.end,
                ann.strand,
                &ann.feature_type,
                main_isoform,
                parent_annotation_id,
            ])
            .with_context(|| {
                format!(
                    "Failed to append annotation: contig_id={}, start={}, end={}, type={}",
                    ann.contig_id, ann.start, ann.end, ann.feature_type
                )
            })?;

        // All qualifiers → KV table (skip "translation" — stored separately in Annotation_sequence)
        for (key, value) in &ann.qualifiers {
            if key == "translation" {
                continue;
            }
            qual_appender
                .append_row(params![annotation_id, key, value])
                .with_context(|| {
                    format!(
                        "Failed to append qualifier: annotation_id={}, key={}",
                        annotation_id, key
                    )
                })?;
        }

        // Heavy sequence row
        seq_appender
            .append_row(params![
                annotation_id,
                &ann.nucleotide_sequence,
                &ann.protein_sequence,
                s_sites,
                n_sites,
            ])
            .with_context(|| {
                format!(
                    "Failed to append annotation sequence: annotation_id={}",
                    annotation_id
                )
            })?;

        // Segment rows — only for features with two or more parts.
        for (seg_idx, &(seg_start, seg_end)) in ann.segments.iter().enumerate() {
            seg_appender
                .append_row(params![annotation_id, seg_idx as i64, seg_start, seg_end])
                .with_context(|| {
                    format!(
                        "Failed to append segment: annotation_id={}, seg_idx={}",
                        annotation_id, seg_idx
                    )
                })?;
        }
    }

    ann_appender
        .flush()
        .context("Failed to flush Contig_annotation_core appender")?;
    qual_appender
        .flush()
        .context("Failed to flush Annotation_qualifier appender")?;
    seq_appender
        .flush()
        .context("Failed to flush Annotation_sequence appender")?;
    seg_appender
        .flush()
        .context("Failed to flush Annotation_segments appender")?;

    // Check if any locus_tag appears more than once (isoforms present)
    let mut locus_tag_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    for ann in annotations {
        if let Some(locus_tag) = ann.qualifiers.get("locus_tag") {
            *locus_tag_counts.entry(locus_tag.clone()).or_insert(0) += 1;
        }
    }
    let has_isoforms = locus_tag_counts.values().any(|&count| count > 1);

    // Insert constants into Constants_for_plotting table
    conn.execute(
        "INSERT INTO Constants_for_plotting (Constant, Status) VALUES ('isoforms', ?)",
        params![has_isoforms],
    )
    .context("Failed to insert constants")?;

    // Always insert generic color template
    let mut template_id: i64 = 1;
    conn.execute(
        "INSERT INTO Color_templates (Template_id, Template_name) VALUES (?, 'generic')",
        params![template_id],
    )
    .context("Failed to insert generic template")?;

    // (Qualifier_name, Operator, Qualifier_value, Color)
    let generic_rules: Vec<(&str, &str, &str, &str)> = vec![
        ("Type", "=", "gene", "#555555"),
        ("Type", "=", "mRNA", "#777777"),
        ("Type", "=", "CDS", "#cccccc"),
        ("Type", "=", "tRNA", "#66c2a5"),
        ("Type", "=", "tmRNA", "#99d8c9"),
        ("Type", "=", "rRNA", "#238b45"),
        ("Type", "=", "ncRNA", "#33a02c"),
        ("Type", "=", "precursor_RNA", "#a1d99b"),
        ("Type", "=", "misc_RNA", "#74c476"),
        ("Type", "=", "exon", "#fdae61"),
        ("Type", "=", "5'UTR", "#fee08b"),
        ("Type", "=", "3'UTR", "#f46d43"),
        ("Type", "=", "repeat_region", "#6a3d9a"),
        ("Type", "=", "mobile_element", "#cab2d6"),
        ("Type", "=", "misc_feature", "#3c5bfe"),
        ("Type", "=", "gap", "#e5049c"),
        ("Type", "=", "pseudogene", "#e31a1c"),
    ];
    for (i, (qname, op, qvalue, color)) in generic_rules.iter().enumerate() {
        conn.execute(
            "INSERT INTO Color_rules (Rule_id, Template_id, Qualifier_name, Operator, Qualifier_value, Color) VALUES (?, ?, ?, ?, ?, ?)",
            params![i as i64 + 1, template_id, qname, op, qvalue, color],
        )
        .context("Failed to insert generic color rule")?;
    }

    // Insert pharokka color template only if pharokka functions detected
    if has_pharokka_function {
        template_id = 2;
        conn.execute(
            "INSERT INTO Color_templates (Template_id, Template_name) VALUES (?, 'pharokka')",
            params![template_id],
        )
        .context("Failed to insert pharokka template")?;

        // (Qualifier_name, Operator, Qualifier_value, Color)
        let pharokka_rules: Vec<(&str, &str, &str, &str)> = vec![
            ("function", "=", "vfdb_card", "#FF0000"),
            ("function", "=", "unknown function", "#AAAAAA"),
            ("function", "=", "other", "#4deeea"),
            ("function", "=", "tail", "#74ee15"),
            ("function", "=", "transcription regulation", "#ffe700"),
            ("function", "=", "dna, rna and nucleotide metabolism", "#f000ff"),
            ("function", "=", "lysis", "#001eff"),
            ("function", "=", "moron, auxiliary metabolic gene and host takeover", "#8900ff"),
            ("function", "=", "integration and excision", "#E0B0FF"),
            ("function", "=", "head and packaging", "#ff008d"),
            ("function", "=", "connector", "#5A5A5A"),
        ];
        let rule_offset = generic_rules.len() as i64;
        for (i, (qname, op, qvalue, color)) in pharokka_rules.iter().enumerate() {
            conn.execute(
                "INSERT INTO Color_rules (Rule_id, Template_id, Qualifier_name, Operator, Qualifier_value, Color) VALUES (?, ?, ?, ?, ?, ?)",
                params![i as i64 + 1 + rule_offset, template_id, qname, op, qvalue, color],
            )
            .context("Failed to insert pharokka color rule")?;
        }
    }

    Ok(())
}

/// Insert contig-level qualifiers (from GenBank "source" features) into the Contig_qualifier table.
fn insert_contig_qualifiers(conn: &Connection, contig_qualifiers: &[(i64, HashMap<String, String>)]) -> Result<()> {
    if contig_qualifiers.is_empty() {
        return Ok(());
    }
    let mut appender = conn.appender("Contig_qualifier")
        .context("Failed to create Contig_qualifier appender")?;
    for (contig_id, quals) in contig_qualifiers {
        for (key, value) in quals {
            appender.append_row(params![contig_id, key, value])
                .with_context(|| format!(
                    "Failed to append contig qualifier: contig_id={}, key={}",
                    contig_id, key
                ))?;
        }
    }
    appender.flush().context("Failed to flush Contig_qualifier appender")?;
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
         FROM Contig_annotation_core
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
        // Special case: contig-level data is now stored in Contig_blob table
        // Per-sample features use Feature_blob table with Feature_id
        let table_name = match v.name {
            "direct_repeat_count" => "Contig_blob".to_string(),
            "direct_repeat_identity" => "Contig_blob".to_string(),
            "inverted_repeat_count" => "Contig_blob".to_string(),
            "inverted_repeat_identity" => "Contig_blob".to_string(),
            "gc_content" => "Contig_blob".to_string(),
            "gc_skew" => "Contig_blob".to_string(),
            _ => "Feature_blob".to_string(),
        };

        conn.execute(
            "INSERT INTO Variable (Variable_id, Variable_name, Subplot, Module, Module_order, \"Type\", Encoding, Value_scale, Color, Alpha, Fill_alpha, \"Size\", Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15)",
            params![
                (i + 1) as i64,
                v.name,
                v.subplot,
                v.module,
                v.module_order,
                v.plot_type.as_str(),
                v.encoding.as_str(),
                v.value_scale.as_str(),
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

/// Create materialized tables and views after all data is inserted.
fn create_views(conn: &Connection, has_bam: bool) -> Result<()> {
    // All per-sample features (including primary_reads) are now stored in Feature_blob.
    // Previously, some features were stored in separate Feature_* tables.

    // Materialized repeat tables (sweep-line algorithm, computed once after all repeat data is inserted).
    // These produce (Contig_id, First_position, Last_position, Value) like other Contig_* tables,
    // allowing the standard binning logic in get_feature_data() to work with repeats.

    // Optimized repeat count tables using event-based running sum (O(n log n) instead of O(segments × repeats))
    conn.execute(
        "CREATE TABLE Contig_direct_repeat_count AS
        WITH events AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos, 1 AS delta
            FROM Contig_directRepeats
            UNION ALL
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos, -1 AS delta
            FROM Contig_directRepeats
        ),
        grouped AS (
            SELECT Contig_id, pos, SUM(delta) AS net FROM events GROUP BY Contig_id, pos
        ),
        running AS (
            SELECT Contig_id, pos,
                   SUM(net) OVER (PARTITION BY Contig_id ORDER BY pos) AS count,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) AS next_pos
            FROM grouped
        )
        SELECT Contig_id, pos AS First_position, next_pos - 1 AS Last_position, count AS Value
        FROM running WHERE next_pos IS NOT NULL AND count > 0",
        [],
    )
    .context("Failed to create Contig_direct_repeat_count table")?;

    conn.execute(
        "CREATE TABLE Contig_inverted_repeat_count AS
        WITH events AS (
            SELECT Contig_id, LEAST(Position1, Position2) AS pos, 1 AS delta
            FROM Contig_invertedRepeats
            UNION ALL
            SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos, -1 AS delta
            FROM Contig_invertedRepeats
        ),
        grouped AS (
            SELECT Contig_id, pos, SUM(delta) AS net FROM events GROUP BY Contig_id, pos
        ),
        running AS (
            SELECT Contig_id, pos,
                   SUM(net) OVER (PARTITION BY Contig_id ORDER BY pos) AS count,
                   LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) AS next_pos
            FROM grouped
        )
        SELECT Contig_id, pos AS First_position, next_pos - 1 AS Last_position, count AS Value
        FROM running WHERE next_pos IS NOT NULL AND count > 0",
        [],
    )
    .context("Failed to create Contig_inverted_repeat_count table")?;

    // Optimized repeat identity tables using explicit JOIN instead of correlated subquery
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
        SELECT s.Contig_id, s.seg_start AS First_position, s.seg_end AS Last_position,
               MAX(r.Pident) / 100.0 AS Value
        FROM segments s
        JOIN Contig_directRepeats r ON r.Contig_id = s.Contig_id
          AND LEAST(r.Position1, r.Position2) <= s.seg_end
          AND GREATEST(r.Position1, r.Position2) >= s.seg_start
        WHERE s.seg_end IS NOT NULL
        GROUP BY s.Contig_id, s.seg_start, s.seg_end",
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
        SELECT s.Contig_id, s.seg_start AS First_position, s.seg_end AS Last_position,
               MAX(r.Pident) / 100.0 AS Value
        FROM segments s
        JOIN Contig_invertedRepeats r ON r.Contig_id = s.Contig_id
          AND LEAST(r.Position1, r.Position2) <= s.seg_end
          AND GREATEST(r.Position1, r.Position2) >= s.seg_start
        WHERE s.seg_end IS NOT NULL
        GROUP BY s.Contig_id, s.seg_start, s.seg_end",
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
             ROUND(p.Coverage_coefficient_of_variation / 1000000.0, 2) AS Coverage_coefficient_of_variation,
             ROUND(p.Coverage_relative_coverage_roughness / 1000000.0, 4) AS Relative_coverage_roughness
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
             t.Circularising_reads AS Circularising_reads,
             t.Circularising_reads_percentage AS Circularising_reads_prevalence,
             t.Median_circularising_len AS Median_circularising_len,
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
    // Includes aggregated diagnostics columns from Phage_termini
    conn.execute(
        "CREATE VIEW Explicit_phage_mechanisms AS
         SELECT
             c.Contig_name,
             s.Sample_name,
             m.Packaging_mechanism,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM Phage_termini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'start' AND pt.Passed_PoissonTest = 'yes' AND pt.Passed_ClippingTest = 'yes'),
                 ''
             ) AS Left_termini,
             m.Median_left_termini_clippings,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM Phage_termini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'end' AND pt.Passed_PoissonTest = 'yes' AND pt.Passed_ClippingTest = 'yes'),
                 ''
             ) AS Right_termini,
             m.Median_right_termini_clippings,
             CASE WHEN m.Duplication = true THEN 'DTR'
                  WHEN m.Duplication = false THEN 'ITR'
                  ELSE NULL END AS Duplication,
             COALESCE(
                 (SELECT COUNT(*)
                  FROM Phage_termini pt
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
         LEFT JOIN Phage_mechanisms m ON p.Contig_id = m.Contig_id AND p.Sample_id = m.Sample_id",
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
         FROM Phage_termini t
         JOIN Phage_mechanisms m ON t.Packaging_id = m.Packaging_id
         JOIN Contig c ON m.Contig_id = c.Contig_id
         JOIN Sample s ON m.Sample_id = s.Sample_id",
    )
    .context("Failed to create Explicit_phage_termini VIEW")?;
    } // end if has_bam (Explicit_* views)

    // Convert materialized repeat tables to Contig_blob
    DbWriter::write_repeat_features_to_blob_static(conn)?;

    // Drop intermediate materialized tables after conversion to Contig_blob
    // Keep base tables (Contig_directRepeats, Contig_invertedRepeats) for reference
    for table in &[
        "Contig_direct_repeat_count", "Contig_inverted_repeat_count",
        "Contig_direct_repeat_identity", "Contig_inverted_repeat_identity",
    ] {
        conn.execute(&format!("DROP TABLE IF EXISTS {}", table), [])
            .context(format!("Failed to drop intermediate table {}", table))?;
    }
    eprintln!("Dropped intermediate materialized tables (data now in Contig_blob)");

    Ok(())
}

/// Delete Variable entries for features that don't have data in Feature_blob.
fn cleanup_unused_variables(conn: &Connection) -> Result<()> {
    // Get all variable names and their table names
    let mut stmt = conn.prepare("SELECT Variable_name, Feature_table_name FROM Variable")?;
    let rows: Vec<(String, String)> = stmt
        .query_map([], |row| Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?)))?
        .filter_map(|r| r.ok())
        .collect();

    let mut deleted_count = 0;
    for (var_name, table_name) in rows {
        if table_name == "Feature_blob" {
            // Check if any row exists for this feature's Feature_id
            if let Some(fid) = feature_name_to_id(&var_name) {
                let has_data: bool = conn
                    .query_row(
                        "SELECT EXISTS(SELECT 1 FROM Feature_blob WHERE Feature_id = ?1)",
                        params![fid as i32],
                        |row| row.get(0),
                    )
                    .unwrap_or(false);
                if has_data {
                    continue;
                }
            }
        } else if table_name == "Contig_blob" {
            // Contig-level variables — always kept
            continue;
        }

        conn.execute("DELETE FROM Variable WHERE Variable_name = ?1", params![var_name])?;
        deleted_count += 1;
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
        ("Phage_termini", &["Explicit_phage_termini"]),
        ("Phage_mechanisms", &["Explicit_phage_mechanisms"]),
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
/// Contains per-window GC content and GC skew values with statistics.
#[derive(Clone, Debug)]
pub struct GCContentData {
    pub contig_name: String,
    pub gc_values: Vec<u8>,    // raw GC percentages, one per window
    pub skew_values: Vec<i16>, // raw GC skew × 100, one per window
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
    pub median_circularising_len: Option<i64>,
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
