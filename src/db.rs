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
//! - Zone maps for column scans; explicit B-tree indexes for point lookups
//! - Better write performance with batch inserts
//! - Simpler concurrency model (single writer, multiple readers)

use anyhow::{Context, Result};
use chrono::Local;
use duckdb::{params, Connection};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;

use crate::features::translate_codon;
use crate::gc_content::{GCSkewStats, GCStats, DEFAULT_GC_CONTENT_WINDOW_SIZE, DEFAULT_GC_SKEW_WINDOW_SIZE};
use crate::types::{feature_name_to_id, ContigInfo, FeatureAnnotation, PackagingData, PresenceData, VARIABLES};
// Re-export new metric data structs (defined below in this file)

/// Wall times and actual encoded payload written by a sample. Collected only
/// with --diagnose; blob bytes exclude DuckDB/WAL overhead and disk compression.
#[derive(Default, Clone)]
pub struct SampleWriteTimings {
    pub coverage_ns: u64,
    pub metrics_ns: u64,
    pub checkpoint_ns: u64,
    pub blob_append_ns: u64,
    pub zoom_flush_ns: u64,
    pub chunk_flush_ns: u64,
    pub zoom_rows: u64,
    pub chunk_rows: u64,
    pub zoom_bytes: u64,
    pub chunk_bytes: u64,
}

impl SampleWriteTimings {
    pub(crate) fn add(&mut self, other: &Self) {
        self.coverage_ns += other.coverage_ns;
        self.metrics_ns += other.metrics_ns;
        self.checkpoint_ns += other.checkpoint_ns;
        self.blob_append_ns += other.blob_append_ns;
        self.zoom_flush_ns += other.zoom_flush_ns;
        self.chunk_flush_ns += other.chunk_flush_ns;
        self.zoom_rows += other.zoom_rows;
        self.chunk_rows += other.chunk_rows;
        self.zoom_bytes += other.zoom_bytes;
        self.chunk_bytes += other.chunk_bytes;
    }
}

fn timed_start(enabled: bool) -> Option<std::time::Instant> {
    enabled.then(std::time::Instant::now)
}
fn elapsed_ns(start: Option<std::time::Instant>) -> u64 {
    start.map(|t| t.elapsed().as_nanos() as u64).unwrap_or(0)
}

pub(crate) fn detect_available_memory_gb() -> (u64, &'static str) {
    let mut limits = Vec::new();
    // Consider scheduler, cgroup ancestors, and physical memory together.
    if let Ok(val) = std::env::var("SLURM_MEM_PER_NODE") {
        if let Ok(mb) = val.trim().parse::<u64>() {
            if mb > 0 {
                limits.push((mb / 1024, "SLURM_MEM_PER_NODE"));
            }
        }
    }
    if let (Ok(per_cpu), Ok(cpus)) = (
        std::env::var("SLURM_MEM_PER_CPU"),
        std::env::var("SLURM_JOB_CPUS_PER_NODE"),
    ) {
        if let (Ok(mb), Ok(n)) = (per_cpu.trim().parse::<u64>(), cpus.trim().split(|c: char| !c.is_ascii_digit()).next().unwrap_or("").parse::<u64>()) {
            if mb > 0 && n > 0 {
                limits.push((mb.saturating_mul(n) / 1024, "SLURM_MEM_PER_CPU"));
            }
        }
    }

    // 2. Process-aware cgroups — read /proc/self/cgroup to find our actual cgroup path
    if let Ok(cgroup_info) = std::fs::read_to_string("/proc/self/cgroup") {
        // cgroups v2: line like "0::/<path>"
        for line in cgroup_info.lines() {
            if line.starts_with("0::/") {
                let path = line.trim_start_matches("0::");
                if let Some(gb) = read_cgroupv2_limit(path) {
                    limits.push((gb, "cgroup-v2"));
                }
            }
        }
        // cgroups v1: line like "N:memory:/<path>"
        for line in cgroup_info.lines() {
            let parts: Vec<&str> = line.splitn(3, ':').collect();
            if parts.len() == 3 && parts[1].split(',').any(|c| c == "memory") {
                let path = parts[2];
                if let Some(gb) = read_cgroupv1_limit(path) {
                    limits.push((gb, "cgroup-v1"));
                }
            }
        }
    }

    // 3. /proc/meminfo (bare metal, no cgroup)
    if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
        for line in s.lines() {
            if line.starts_with("MemTotal:") {
                if let Some(kb_str) = line.split_whitespace().nth(1) {
                    if let Ok(kb) = kb_str.parse::<u64>() {
                        limits.push((kb / (1024 * 1024), "procmeminfo"));
                    }
                }
            }
        }
    }

    limits.into_iter().min_by_key(|(gb, _)| *gb).unwrap_or((16, "default"))
}

fn read_cgroupv2_limit(cgroup_path: &str) -> Option<u64> {
    let mut dir = format!("/sys/fs/cgroup{}", cgroup_path);
    let mut limit: Option<u64> = None;
    loop {
        let file = format!("{}/memory.max", dir);
        if let Ok(s) = std::fs::read_to_string(&file) {
            let val = s.trim();
            if val != "max" {
                if let Ok(bytes) = val.parse::<u64>() {
                    let gb = bytes / (1024 * 1024 * 1024);
                    limit = Some(limit.map_or(gb, |old| old.min(gb)));
                }
            }
        }
        if dir == "/sys/fs/cgroup" {
            break;
        }
        dir = match dir.rsplit_once('/') {
            Some((parent, _)) if !parent.is_empty() => parent.to_string(),
            _ => break,
        };
    }
    limit
}

fn read_cgroupv1_limit(cgroup_path: &str) -> Option<u64> {
    let mut dir = format!("/sys/fs/cgroup/memory{}", cgroup_path);
    let mut limit: Option<u64> = None;
    loop {
        let file = format!("{}/memory.limit_in_bytes", dir);
        if let Ok(s) = std::fs::read_to_string(&file) {
            let bytes = s.trim().parse::<u64>().unwrap_or(0);
            if bytes > 0 && bytes < u64::MAX / 2 {
                let gb = bytes / (1024 * 1024 * 1024);
                limit = Some(limit.map_or(gb, |old| old.min(gb)));
            }
        }
        if dir == "/sys/fs/cgroup/memory" {
            break;
        }
        dir = match dir.rsplit_once('/') {
            Some((parent, _)) if !parent.is_empty() => parent.to_string(),
            _ => break,
        };
    }
    limit
}

/// Keep a reserve for Rust allocations, allocator overhead, and worker processes.
fn memory_budget_mib(total: u64, rss: u64, requested: Option<u64>) -> Result<u64> {
    let reserve = (total / 5).max(256);
    let headroom = total.saturating_sub(rss).saturating_sub(reserve);
    anyhow::ensure!(headroom >= 128, "Insufficient memory headroom for DuckDB: total {total} MiB, RSS {rss} MiB");
    let desired = requested.unwrap_or((total / 4).min(16 * 1024));
    anyhow::ensure!(desired > 0, "DuckDB memory budget must be positive");
    Ok(desired.min(headroom).max(1))
}

pub(crate) fn duckdb_memory_budget() -> Result<(u64, u64, &'static str)> {
    let (total_gb, source) = detect_available_memory_gb();
    let rss_mib = std::fs::read_to_string("/proc/self/status").ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmRSS:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<u64>().ok())).unwrap_or(0) / 1024;
    let requested = std::env::var("THEBIGBAM_DUCKDB_MEMORY_GB").ok()
        .map(|v| v.parse::<u64>().context("THEBIGBAM_DUCKDB_MEMORY_GB must be a positive integer"))
        .transpose()?.map(|gb| gb.saturating_mul(1024));
    Ok((memory_budget_mib(total_gb * 1024, rss_mib, requested)?, total_gb, source))
}

pub(crate) fn configure_for_bulk_writes(conn: &Connection) -> Result<()> {
    conn.execute_batch("SET threads TO 1; SET preserve_insertion_order=false; SET wal_autocheckpoint='1GB'")?;
    let (mib, total_gb, source) = duckdb_memory_budget()?;
    eprintln!("DuckDB memory_limit: {mib} MiB (total {total_gb} GiB from {source}; Rust/worker reserve applied)");
    conn.execute_batch(&format!("SET memory_limit='{mib}MiB'"))
        .context("Failed to set DuckDB memory limit")?;
    Ok(())
}

fn build_contig_length_cache(conn: &Connection) -> HashMap<i64, i32> {
    let mut map = HashMap::new();
    if let Ok(mut stmt) = conn.prepare("SELECT Contig_id, Contig_length FROM Contig") {
        if let Ok(rows) = stmt.query_map([], |row| Ok((row.get::<_, i64>(0)?, row.get::<_, i32>(1)?))) {
            for row in rows.flatten() {
                map.insert(row.0, row.1);
            }
        }
    }
    map
}

fn build_terminase_cache(conn: &Connection) -> HashMap<i64, Vec<(i32, i32)>> {
    let mut map: HashMap<i64, Vec<(i32, i32)>> = HashMap::new();
    let sql = r#"
        SELECT ca."Contig_id", ca."Start", ca."End"
        FROM Contig_annotation_core ca
        JOIN Annotation_qualifier aq ON aq.Annotation_id = ca.Annotation_id
        WHERE aq."Key" = 'product' AND aq."Value" LIKE '%terminase%'
    "#;
    if let Ok(mut stmt) = conn.prepare(sql) {
        if let Ok(rows) = stmt.query_map([], |row| {
            Ok((row.get::<_, i64>(0)?, row.get::<_, i32>(1)?, row.get::<_, i32>(2)?))
        }) {
            for row in rows.flatten() {
                map.entry(row.0).or_default().push((row.1, row.2));
            }
        }
    }
    map
}

/// Thread-safe database connection wrapper for sequential writes.
#[derive(Default)]
struct BlobSize {
    zoom_rows: u64, chunk_rows: u64, zoom_bytes: u64, chunk_bytes: u64, skipped: u64,
}
impl BlobSize {
    fn record(&mut self, encoded: &crate::blob::EncodedBlob) {
        self.zoom_rows += 1;
        self.zoom_bytes += encoded.zoom.len() as u64;
        for chunk in &encoded.chunks {
            if encoded.should_store_chunk(chunk) {
                self.chunk_rows += 1;
                self.chunk_bytes += chunk.len() as u64;
            } else { self.skipped += 1; }
        }
    }
    fn add(&mut self, other: &Self) {
        self.zoom_rows += other.zoom_rows; self.chunk_rows += other.chunk_rows;
        self.zoom_bytes += other.zoom_bytes; self.chunk_bytes += other.chunk_bytes;
        self.skipped += other.skipped;
    }
}
struct CalculationSizeReport {
    database: std::path::PathBuf,
    blobs: BTreeMap<(String, String), BlobSize>,
    imports: BTreeMap<String, (u64, u64)>,
}
impl CalculationSizeReport {
    fn write(&self, complete: bool) -> Result<()> {
        use std::io::Write;
        let mut path = self.database.as_os_str().to_os_string();
        path.push(".size.tsv");
        let path = std::path::PathBuf::from(path);
        let mut temporary = path.as_os_str().to_os_string(); temporary.push(".tmp");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&temporary)?);
        writeln!(f, "# scope=current calculation only; existing rows in extended databases excluded")?;
        writeln!(f, "# status={}", if complete { "complete" } else { "in_progress" })?;
        writeln!(f, "# database_file_bytes={}", std::fs::metadata(&self.database).map(|m| m.len()).unwrap_or(0))?;
        let mut wal = self.database.as_os_str().to_os_string(); wal.push(".wal");
        writeln!(f, "# wal_file_bytes={}", std::fs::metadata(wal).map(|m| m.len()).unwrap_or(0))?;
        writeln!(f, "# blob bytes are encoded payload, not physical table/index sizes; import bytes count sequence/qualifier text only")?;
        writeln!(f, "table\tfeature\tzoom_rows\tchunk_rows\tzoom_payload_bytes\tchunk_payload_bytes\tmean_zoom_bytes\tmean_chunk_bytes\tempty_sparse_chunks_skipped\timport_rows\timport_text_bytes")?;
        for ((table, feature), b) in &self.blobs {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t0\t0", table, feature,
                b.zoom_rows, b.chunk_rows, b.zoom_bytes, b.chunk_bytes,
                b.zoom_bytes as f64 / b.zoom_rows.max(1) as f64,
                b.chunk_bytes as f64 / b.chunk_rows.max(1) as f64, b.skipped)?;
        }
        for (table, (rows, bytes)) in &self.imports {
            writeln!(f, "{}\t-\t0\t0\t0\t0\t0\t0\t0\t{}\t{}", table, rows, bytes)?;
        }
        f.flush()?; drop(f); std::fs::rename(temporary, path)?; Ok(())
    }
}

pub struct DbWriter {
    conn: Mutex<Connection>,
    size_report: Mutex<Option<CalculationSizeReport>>,
    has_bam: bool,
    /// True when the pipeline is running in MAG view; gates creation of MAG
    /// tables, per-MAG explicit views, and the inter-contig blast-hits table.
    is_mag_mode: bool,
    has_annotations: bool,
    contig_name_to_id: HashMap<String, i64>,
    contig_id_to_length: HashMap<i64, i32>,
    terminase_cache: HashMap<i64, Vec<(i32, i32)>>,
    sample_name_to_id: Mutex<HashMap<String, i64>>,
    next_sample_id: Mutex<i64>,
    next_packaging_id: Mutex<i64>,
    next_terminus_id: Mutex<i64>,
}

pub struct FinalizeTimings {
    pub create_views_secs: f64,
    pub cleanup_vars_secs: f64,
    pub drop_empty_secs: f64,
    pub update_counts_secs: f64,
    pub checkpoint_secs: f64,
    /// Contigs with Number_of_samples = 0 after finalize (contig mode only).
    pub zero_coverage_contigs: u64,
    /// MAGs with Number_of_samples = 0 after finalize (MAG mode only).
    pub zero_coverage_mags: u64,
}

impl DbWriter {
    pub(crate) fn enable_size_report(&self, database: &Path) -> Result<bool> {
        let mut report = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?;
        if report.is_some() { return Ok(false); }
        let initial = CalculationSizeReport { database: database.to_path_buf(), blobs: BTreeMap::new(), imports: BTreeMap::new() };
        initial.write(false)?; *report = Some(initial); Ok(true)
    }
    fn flush_size_report(&self, complete: bool) -> Result<()> {
        if let Some(report) = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?.as_ref() { report.write(complete)?; }
        Ok(())
    }
    fn record_blob_sizes(&self, table: &str, sizes: BTreeMap<String, BlobSize>) -> Result<()> {
        if let Some(report) = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?.as_mut() {
            for (feature, size) in sizes { report.blobs.entry((table.to_owned(), feature)).or_default().add(&size); }
        }
        Ok(())
    }
    pub(crate) fn record_import_sizes(&self, contigs: &[ContigInfo], annotations: &[FeatureAnnotation], qualifiers: &[(i64, HashMap<String, String>)]) -> Result<()> {
        let mut guard = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?;
        if let Some(report) = guard.as_mut() {
            let mut add = |table: &str, rows: u64, bytes: u64| { let e = report.imports.entry(table.to_owned()).or_default(); e.0 += rows; e.1 += bytes; };
            add("Contig", contigs.len() as u64, 0);
            add("Contig_sequence", contigs.iter().filter(|c| c.sequence.is_some()).count() as u64, contigs.iter().filter_map(|c| c.sequence.as_ref()).map(|s| s.len() as u64).sum());
            add("Contig_annotation_core", annotations.len() as u64, 0);
            add("Annotation_sequence", annotations.len() as u64, annotations.iter().map(|a| a.nucleotide_sequence.as_ref().map_or(0, |s| s.len() as u64) + a.protein_sequence.as_ref().map_or(0, |s| s.len() as u64)).sum());
            add("Annotation_segments", annotations.iter().map(|a| a.segments.len() as u64).sum(), 0);
            let aq = annotations.iter().flat_map(|a| a.qualifiers.iter()).filter(|(k, _)| k.as_str() != "translation");
            let (mut rows, mut bytes) = (0, 0); for (k,v) in aq { rows += 1; bytes += (k.len()+v.len()) as u64; }
            add("Annotation_qualifier", rows, bytes);
            let (mut rows, mut bytes) = (0, 0); for (_, q) in qualifiers { for (k,v) in q { rows += 1; bytes += (k.len()+v.len()) as u64; } }
            add("Contig_qualifier", rows, bytes);
        }
        Ok(())
    }

    fn initialize_sample_states(conn:&Connection) -> Result<()> {
        conn.execute_batch("CREATE TABLE IF NOT EXISTS Sample_processing (Sample_id INTEGER PRIMARY KEY, State VARCHAR NOT NULL, Error VARCHAR);
            DELETE FROM Sample_processing WHERE Sample_id NOT IN (SELECT Sample_id FROM Sample);
            INSERT INTO Sample_processing SELECT Sample_id, 'complete', NULL FROM Sample ON CONFLICT DO NOTHING")?;
        Ok(())
    }
    pub(crate) fn sample_complete(&self,name:&str)->Result<()> {
        self.lock_conn()?.execute("UPDATE Sample_processing SET State='complete', Error=NULL WHERE State='processing' AND Sample_id=(SELECT Sample_id FROM Sample WHERE Sample_name=?)",params![name])?;
        Ok(())
    }
    pub(crate) fn sample_failed(&self,name:&str,error:&str)->Result<()> {
        let conn=self.lock_conn()?;
        conn.execute("UPDATE Sample_processing SET State='failed', Error=? WHERE State<>'complete' AND Sample_id=(SELECT Sample_id FROM Sample WHERE Sample_name=?)",params![error,name])?;
        let id={let mut stmt=conn.prepare("SELECT s.Sample_id FROM Sample s JOIN Sample_processing p USING(Sample_id) WHERE s.Sample_name=? AND p.State<>'complete'")?;let mut rows=stmt.query(params![name])?;match rows.next()? {Some(row)=>Some(row.get::<_,i64>(0)?),None=>None}};
        if let Some(id)=id {Self::clear_sample_payload(&conn,id)?;}
        Ok(())
    }
    /// Keep the parent Sample row/ID while clearing a retry's partial children.
    /// Each delete commits separately: a kill during cleanup remains retryable.
    fn clear_sample_payload(conn:&Connection,id:i64)->Result<()> {
        let exists=|table:&str| -> Result<bool> { Ok(conn.query_row("SELECT EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name=?)",params![table],|r|r.get(0))?) };
        if exists("Phage_termini")? && exists("Phage_mechanisms")? {
            conn.execute("DELETE FROM Phage_termini WHERE Packaging_id IN (SELECT Packaging_id FROM Phage_mechanisms WHERE Sample_id=?)",params![id])?;
        }
        let mut stmt=conn.prepare("SELECT c.table_name FROM information_schema.columns c JOIN information_schema.tables t USING (table_catalog,table_schema,table_name) WHERE c.column_name='Sample_id' AND t.table_type='BASE TABLE' AND c.table_name NOT IN ('Sample','Sample_processing')")?;
        let tables=stmt.query_map([],|r|r.get::<_,String>(0))?.collect::<duckdb::Result<Vec<_>>>()?;
        drop(stmt);
        for table in tables {
            let quoted=table.replace('"', "\"\"");
            conn.execute(&format!("DELETE FROM \"{}\" WHERE Sample_id=?",quoted),params![id])?;
        }
        Ok(())
    }

    pub(crate) fn reset_reference_blast(&self)->Result<()> {
        let conn=self.lock_conn()?;
        conn.execute_batch("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_state','processing');
            DELETE FROM Database_metadata WHERE Key='Blast_parameters';
            DELETE FROM Contig_directRepeats; DELETE FROM Contig_invertedRepeats; DELETE FROM Contig_blast_hits")?;
        for feature in ["direct_repeat_count","inverted_repeat_count","direct_repeat_identity","inverted_repeat_identity","hit_count_within_mag","hit_identity_within_mag","hit_count","hit_identity"] {
            if let Some(id)=feature_name_to_id(feature) {
                conn.execute("DELETE FROM Contig_blob_chunk WHERE Feature_id=?",params![id])?;
                conn.execute("DELETE FROM Contig_blob WHERE Feature_id=?",params![id])?;
            }
        }
        Ok(())
    }
    pub(crate) fn update_duplication_from_database(&self)->Result<()> {
        let conn=self.lock_conn()?;
        let scale=crate::types::get_column_scale("Contig","Duplication_percentage");
        conn.execute("UPDATE Contig SET Duplication_percentage=0",[])?;
        conn.execute(&format!("WITH intervals AS (
            SELECT Contig_id, LEAST(Position1,Position2) AS s, GREATEST(Position1,Position2) AS e FROM Contig_directRepeats
            UNION ALL SELECT Contig_id, LEAST(Position1prime,Position2prime),GREATEST(Position1prime,Position2prime) FROM Contig_directRepeats
            UNION ALL SELECT Contig_id, LEAST(Position1,Position2),GREATEST(Position1,Position2) FROM Contig_invertedRepeats
            UNION ALL SELECT Contig_id, LEAST(Position1prime,Position2prime),GREATEST(Position1prime,Position2prime) FROM Contig_invertedRepeats
        ), previous AS (
            SELECT *, MAX(e) OVER(PARTITION BY Contig_id ORDER BY s,e ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prior_end FROM intervals
        ), islands AS (
            SELECT *, SUM(CASE WHEN prior_end IS NULL OR s>prior_end+1 THEN 1 ELSE 0 END) OVER(PARTITION BY Contig_id ORDER BY s,e) AS island FROM previous
        ), merged AS (SELECT Contig_id,island,MIN(s) AS s,MAX(e) AS e FROM islands GROUP BY Contig_id,island),
        totals AS (SELECT Contig_id,SUM(e-s+1) AS bases FROM merged GROUP BY Contig_id)
        UPDATE Contig SET Duplication_percentage=CAST(ROUND(CAST(t.bases AS DOUBLE)*100.0*{scale}/Contig.Contig_length) AS INTEGER)
        FROM totals t WHERE Contig.Contig_id=t.Contig_id AND Contig.Contig_length>0"),[])?;
        if self.is_mag_mode {
            conn.execute("WITH summary AS (
                SELECT a.MAG_id,CAST(ROUND(SUM(CAST(c.Duplication_percentage AS DOUBLE)*c.Contig_length)/NULLIF(SUM(c.Contig_length),0)) AS INTEGER) AS duplication
                FROM MAG_contigs_association a JOIN Contig c USING(Contig_id)
                WHERE c.Duplication_percentage IS NOT NULL GROUP BY a.MAG_id)
                UPDATE MAG SET Duplication_percentage=s.duplication FROM summary s WHERE MAG.MAG_id=s.MAG_id",[])?;
        }
        Ok(())
    }

    pub(crate) fn write_blast_mode(&self, mode: &str) -> Result<()> {
        let conn = self.lock_conn()?;
        // Legacy booleans describe prior scope; they never authorize a new search.
        conn.execute_batch("INSERT INTO Database_metadata SELECT 'Blast_mode', CASE WHEN Value='true' THEN 'all' ELSE 'none' END FROM Database_metadata WHERE Key='Blast_enabled' ON CONFLICT DO NOTHING")?;
        drop(conn);
        self.lock_conn()?.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_mode', ?)", params![mode])?;
        Ok(())
    }
    pub(crate) fn write_blast_parameters(&self) -> Result<()> {
        self.lock_conn()?.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_parameters', ?)", params![crate::processing::BLAST_PARAMETERS])?;
        Ok(())
    }

    pub(crate) fn prepare_terminase_cache(&mut self, enabled: bool) -> Result<()> {
        if enabled {
            let conn = self.conn.get_mut().map_err(|e| anyhow::anyhow!("{e}"))?;
            self.terminase_cache = build_terminase_cache(conn);
        } else { self.terminase_cache.clear(); }
        Ok(())
    }

    /// Expose the underlying connection to sibling modules (e.g. `mag_blob`)
    /// that need to both read existing blob data and append MAG-scale blobs
    /// under a single lock.
    pub(crate) fn lock_conn(&self) -> Result<std::sync::MutexGuard<'_, Connection>> {
        self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))
    }

    pub fn checkpoint(&self) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        conn.execute("FORCE CHECKPOINT", [])
            .context("Failed to checkpoint database")?;
        self.flush_size_report(false)?;
        Ok(())
    }

    /// Create a new database and return a writer for sequential sample insertion.
    pub fn create(
        db_path: &Path,
        contigs: &[ContigInfo],
        annotations: &[FeatureAnnotation],
        contig_qualifiers: &[(i64, HashMap<String, String>)],
        has_bam: bool,
        is_mag_mode: bool,
    ) -> Result<Self> {
        Self::create_empty_or_populated(db_path, contigs, annotations, contig_qualifiers,
            has_bam, is_mag_mode, !annotations.is_empty())
    }

    /// Create the schema before reading annotation files. Heavy inputs are appended
    /// one file at a time; only the lightweight lookup caches survive each batch.
    pub(crate) fn create_streaming(db_path: &Path, has_bam: bool, is_mag_mode: bool) -> Result<Self> {
        let writer = Self::create_empty_or_populated(db_path, &[], &[], &[], has_bam, is_mag_mode, true)?;
        writer.lock_conn()?.execute_batch(
            "CREATE TABLE Annotation_import_keys (Annotation_id BIGINT, Self_key VARCHAR, Parent_key VARCHAR, Locus_tag VARCHAR)"
        )?;
        Ok(writer)
    }

    fn create_empty_or_populated(
        db_path: &Path, contigs: &[ContigInfo], annotations: &[FeatureAnnotation],
        contig_qualifiers: &[(i64, HashMap<String, String>)],
        has_bam: bool, is_mag_mode: bool, has_annotations: bool,
    ) -> Result<Self> {
        // Remove existing database if present
        let _ = std::fs::remove_file(db_path);

        let conn = Connection::open(db_path)
            .with_context(|| format!("Failed to create database: {}", db_path.display()))?;
        configure_for_bulk_writes(&conn)?;

        // Create tables
        create_core_tables(&conn, has_bam, is_mag_mode, has_annotations)?;
        conn.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_hit_storage','query-v1')", [])?;
        create_variable_tables(&conn)?;

        // Insert contigs
        insert_contigs(&conn, contigs, is_mag_mode)?;

        // Insert contig sequences (from GenBank or FASTA)
        insert_contig_sequences(&conn, contigs)?;

        // Insert codon table (standard genetic code)
        insert_codon_table(&conn)?;

        if !annotations.is_empty() {
            // Insert annotations and their qualifiers
            insert_annotations(&conn, annotations)?;

            // Insert contig-level qualifiers from source features
            insert_contig_qualifiers(&conn, contig_qualifiers)?;

            // Populate Annotated_types table with distinct Type values from Contig_annotation
            populate_annotated_types(&conn)?;
        }

        // Build contig name -> id mapping
        let contig_name_to_id: HashMap<String, i64> = contigs
            .iter()
            .enumerate()
            .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
            .collect();

        let contig_id_to_length: HashMap<i64, i32> = contigs
            .iter()
            .enumerate()
            .map(|(i, c)| ((i + 1) as i64, c.length as i32))
            .collect();

        if has_bam { Self::initialize_sample_states(&conn)?; }
        let terminase_cache = HashMap::new();

        Ok(Self {
            conn: Mutex::new(conn),
            size_report: Mutex::new(None),
            has_bam,
            is_mag_mode,
            has_annotations,
            contig_name_to_id,
            contig_id_to_length,
            terminase_cache,
            sample_name_to_id: Mutex::new(HashMap::new()),
            next_sample_id: Mutex::new(1),
            next_packaging_id: Mutex::new(1),
            next_terminus_id: Mutex::new(1),
        })
    }

    pub(crate) fn append_annotation_batch(
        &mut self, contigs: &[ContigInfo], annotations: &[FeatureAnnotation],
        qualifiers: &[(i64, HashMap<String, String>)], contig_offset: i64, annotation_offset: i64,
    ) -> Result<()> {
        let conn = self.conn.get_mut().map_err(|e| anyhow::anyhow!("Lock poisoned: {e}"))?;
        insert_contigs_with_offset(conn, contigs, self.is_mag_mode, contig_offset)?;
        insert_contig_sequences_with_offset(conn, contigs, contig_offset)?;
        insert_annotation_batch(conn, annotations, annotation_offset)?;
        insert_contig_qualifiers(conn, qualifiers)?;
        {
            let mut keys = conn.appender("Annotation_import_keys")?;
            for (i, ann) in annotations.iter().enumerate() {
                keys.append_row(params![annotation_offset + i as i64 + 1,
                    &ann.self_key, &ann.parent_key, ann.qualifiers.get("locus_tag")])?;
            }
            keys.flush()?;
        }
        for (i, contig) in contigs.iter().enumerate() {
            let id = contig_offset + i as i64 + 1;
            self.contig_name_to_id.insert(contig.name.clone(), id);
            self.contig_id_to_length.insert(id, contig.length as i32);
        }
        self.record_import_sizes(contigs, annotations, qualifiers)?;
        self.flush_size_report(false)?;
        Ok(())
    }

    /// Resolve across all input files using spillable SQL, preserving the original
    /// first-occurrence semantics even when a parent or locus occurs in another batch.
    pub(crate) fn finish_annotation_import(&mut self) -> Result<()> {
        let conn = self.conn.get_mut().map_err(|e| anyhow::anyhow!("Lock poisoned: {e}"))?;
        eprintln!("Resolving annotation parents and isoforms across input files...");
        conn.execute_batch(r#"
            CREATE TEMP TABLE Import_parents AS
                SELECT Self_key, min(Annotation_id) AS Parent_id
                FROM Annotation_import_keys WHERE Self_key IS NOT NULL GROUP BY Self_key;
            UPDATE Contig_annotation_core AS a SET Parent_annotation_id = p.Parent_id
                FROM Annotation_import_keys k LEFT JOIN Import_parents p ON k.Parent_key = p.Self_key
                WHERE a.Annotation_id = k.Annotation_id
                  AND a.Parent_annotation_id IS DISTINCT FROM p.Parent_id;
            CREATE TEMP TABLE Import_main_mrna AS
                SELECT k.Locus_tag, min(a.Annotation_id) AS Main_id
                FROM Contig_annotation_core a JOIN Annotation_import_keys k USING (Annotation_id)
                WHERE a.Type = 'mRNA' AND k.Locus_tag IS NOT NULL GROUP BY k.Locus_tag;
            CREATE TEMP TABLE Import_main_flat AS
                SELECT k.Locus_tag, a.Type, min(a.Annotation_id) AS Main_id
                FROM Contig_annotation_core a JOIN Annotation_import_keys k USING (Annotation_id)
                WHERE a.Parent_annotation_id IS NULL AND a.Type <> 'mRNA' AND k.Locus_tag IS NOT NULL
                GROUP BY k.Locus_tag, a.Type;
            UPDATE Contig_annotation_core AS a SET Main_isoform =
                CASE WHEN k.Locus_tag IS NULL THEN NULL
                     WHEN a.Type = 'mRNA' THEN a.Annotation_id = m.Main_id
                     WHEN a.Parent_annotation_id IS NOT NULL THEN pm.Main_id IS NOT NULL
                     ELSE a.Annotation_id = f.Main_id END
                FROM Annotation_import_keys k
                LEFT JOIN Import_main_mrna m ON k.Locus_tag = m.Locus_tag
                LEFT JOIN Contig_annotation_core ac ON ac.Annotation_id = k.Annotation_id
                LEFT JOIN Import_main_mrna pm ON ac.Parent_annotation_id = pm.Main_id
                LEFT JOIN Import_main_flat f ON k.Locus_tag = f.Locus_tag AND ac.Type = f.Type
                WHERE a.Annotation_id = k.Annotation_id;
            INSERT INTO Constants_for_plotting VALUES ('isoforms',
                EXISTS (SELECT 1 FROM Annotation_import_keys WHERE Locus_tag IS NOT NULL
                        GROUP BY Locus_tag HAVING count(*) > 1))
                ON CONFLICT (Constant) DO UPDATE SET Status = excluded.Status;
            DROP TABLE Import_main_flat;
            DROP TABLE Import_main_mrna;
            DROP TABLE Import_parents;
            DROP TABLE Annotation_import_keys;
        "#)?;
        populate_annotated_types(conn)?;
        self.terminase_cache.clear();
        conn.execute_batch("CHECKPOINT")?;
        Ok(())
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
        configure_for_bulk_writes(&conn)?;

        // Read view_mode from existing metadata (defaults to 'contig' for legacy DBs).
        let is_mag_mode = {
            let mut stmt = conn.prepare(
                "SELECT Value FROM Database_metadata WHERE Key = 'View_mode'",
            )?;
            let mut rows = stmt.query([])?;
            match rows.next()? {
                Some(r) => r.get::<_, String>(0).unwrap_or_else(|_| "contig".to_string()) == "mag",
                None => false,
            }
        };

        if has_bam {
            create_mapping_tables(&conn, is_mag_mode)
                .context("Failed to restore mapping tables for extend mode")?;
        }

        create_intercontig_hit_table(&conn)?;
        conn.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_hit_storage','query-v1')", [])?;
        create_variable_tables(&conn)?;
        for variable in VARIABLES {
            conn.execute("INSERT OR IGNORE INTO Column_scales VALUES (?, 'Value', ?)",
                params![variable.name, variable.value_scale.multiplier() as i32])?;
        }

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
            "Explicit_coverage_per_MAG", "Explicit_misassembly_per_MAG",
            "Explicit_microdiversity_per_MAG", "Contig_blast_hits_symmetric", "Contig_blast_hits_query",
        ] {
            conn.execute(&format!("DROP VIEW IF EXISTS {}", view), [])?;
        }

        // Initialize next packaging/terminus IDs from existing data
        let next_packaging_id: i64 = conn.query_row(
            "SELECT COALESCE(MAX(Packaging_id), 0) + 1 FROM Phage_mechanisms",
            [], |row| row.get(0),
        ).unwrap_or(1);
        let next_terminus_id: i64 = conn.query_row(
            "SELECT COALESCE(MAX(Terminus_id), 0) + 1 FROM Phage_termini",
            [], |row| row.get(0),
        ).unwrap_or(1);

        let contig_id_to_length = build_contig_length_cache(&conn);
        if has_bam { Self::initialize_sample_states(&conn)?; }
        let terminase_cache = HashMap::new();

        let has_annotations: bool = conn.query_row(
            "SELECT COUNT(*) > 0 FROM information_schema.tables WHERE table_name = 'Contig_annotation_core'",
            [], |r| r.get(0),
        ).unwrap_or(false);

        Ok(Self {
            conn: Mutex::new(conn),
            size_report: Mutex::new(None),
            has_bam,
            is_mag_mode,
            has_annotations,
            contig_name_to_id,
            contig_id_to_length,
            terminase_cache,
            sample_name_to_id: Mutex::new(sample_name_to_id),
            next_sample_id: Mutex::new(max_sample_id + 1),
            next_packaging_id: Mutex::new(next_packaging_id),
            next_terminus_id: Mutex::new(next_terminus_id),
        })
    }

    /// Returns whether this database is in MAG view mode (derived from the
    /// `View_mode` metadata row on open, or from the create-time flag).
    pub fn is_mag_mode(&self) -> bool {
        self.is_mag_mode
    }

    /// Insert a sample and return its ID.
    pub fn insert_sample(&self, sample_name: &str, sequencing_type: &str, total_reads: u64, mapped_reads: u64, circular_mapping: bool) -> Result<i64> {
        let conn=self.lock_conn()?;
        let existing:Option<(i64,String)>= {
            let mut stmt=conn.prepare("SELECT s.Sample_id, COALESCE(p.State,'complete') FROM Sample s LEFT JOIN Sample_processing p USING(Sample_id) WHERE Sample_name=?")?;
            let mut rows=stmt.query(params![sample_name])?;
            match rows.next()? { Some(row)=>Some((row.get(0)?,row.get(1)?)),None=>None }
        };
        let sample_id=if let Some((id,state))=existing {
            anyhow::ensure!(state!="complete","Sample '{}' is already complete",sample_name);
            conn.execute("UPDATE Sample_processing SET State='processing',Error=NULL WHERE Sample_id=?",params![id])?;
            Self::clear_sample_payload(&conn,id)?;
            conn.execute("UPDATE Sample SET Sequencing_type=?,Number_of_reads=?,Number_of_mapped_reads=?,Circular_mapping=? WHERE Sample_id=?",params![sequencing_type,total_reads as i64,mapped_reads as i64,circular_mapping,id])?;
            id
        } else {
            let mut next=self.next_sample_id.lock().map_err(|e|anyhow::anyhow!("{e}"))?;
            let id=*next; *next+=1;
            // State is inserted first: a kill between these statements cannot
            // cause a subsequently inserted Sample row to be mistaken for complete.
            conn.execute("INSERT INTO Sample_processing VALUES (?,'processing',NULL)",params![id])?;
            conn.execute("INSERT INTO Sample (Sample_id,Sample_name,Sequencing_type,Number_of_reads,Number_of_mapped_reads,Circular_mapping) VALUES (?,?,?,?,?,?)",params![id,sample_name,sequencing_type,total_reads as i64,mapped_reads as i64,circular_mapping])?;
            id
        };
        drop(conn);

        let mut sample_map = self.sample_name_to_id.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        sample_map.insert(sample_name.to_string(), sample_id);

        Ok(sample_id)
    }

    pub(crate) fn update_sample_read_counts(&self, sample_id: i64, total_reads: u64, mapped_reads: u64) -> Result<()> {
        let conn = self.lock_conn()?;
        conn.execute("UPDATE Sample SET Number_of_reads = ?, Number_of_mapped_reads = ? WHERE Sample_id = ?",
            params![total_reads as i64, mapped_reads as i64, sample_id])?;
        Ok(())
    }

    /// Write all data for a sample (coverage, packaging, misassembly/microdiversity/side_misassembly/topology, feature blobs).
    pub fn write_sample_data(
        &self, sample_name: &str, presences: &[PresenceData], packaging: &[PackagingData],
        misassembly: &[MisassemblyData], microdiversity: &[MicrodiversityData],
        side_misassembly: &[SideMisassemblyData], topology: &[TopologyData],
        feature_blobs: &[(String, i64, crate::blob::EncodedBlob)], circular: bool,
    ) -> Result<()> {
        self.write_sample_data_timed(sample_name, presences, packaging, misassembly,
            microdiversity, side_misassembly, topology, feature_blobs, circular, false).map(|_| ())
    }

    pub fn write_sample_data_timed(
        &self,
        sample_name: &str,
        presences: &[PresenceData],
        packaging: &[PackagingData],
        misassembly: &[MisassemblyData],
        microdiversity: &[MicrodiversityData],
        side_misassembly: &[SideMisassemblyData],
        topology: &[TopologyData],
        feature_blobs: &[(String, i64, crate::blob::EncodedBlob)],
        circular: bool,
        enable_timing: bool,
    ) -> Result<SampleWriteTimings> {
        self.write_sample_batch_timed(sample_name, presences, packaging, misassembly,
            microdiversity, side_misassembly, topology, feature_blobs, circular, enable_timing, true)
    }

    pub(crate) fn write_sample_batch_timed(
        &self,
        sample_name: &str,
        presences: &[PresenceData],
        packaging: &[PackagingData],
        misassembly: &[MisassemblyData],
        microdiversity: &[MicrodiversityData],
        side_misassembly: &[SideMisassemblyData],
        topology: &[TopologyData],
        feature_blobs: &[(String, i64, crate::blob::EncodedBlob)],
        circular: bool,
        enable_timing: bool,
        checkpoint_before_blobs: bool,
    ) -> Result<SampleWriteTimings> {
        let mut timings = SampleWriteTimings::default();
        let sample_id = {
            let sample_map = self.sample_name_to_id.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
            *sample_map.get(sample_name).context("Sample not found")?
        };

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

        // Appenders handle their own transactions - no explicit BEGIN/COMMIT needed
        // Insert coverage data
        let start = timed_start(enable_timing);
        self.write_presences(&conn, sample_id, presences)?;
        timings.coverage_ns = elapsed_ns(start);
        let start = timed_start(enable_timing);

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
        timings.metrics_ns = elapsed_ns(start);

        // Flush non-blob WAL data before the heavy blob writes
        let start = timed_start(enable_timing && checkpoint_before_blobs);
        if checkpoint_before_blobs { conn.execute("FORCE CHECKPOINT", []).ok(); }
        timings.checkpoint_ns = elapsed_ns(start);

        // Insert feature BLOBs (compressed format)
        self.write_feature_blobs_timed(&conn, sample_id, feature_blobs, enable_timing, &mut timings)?;

        Ok(timings)
    }

    fn write_presences(&self, conn: &Connection, sample_id: i64, presences: &[PresenceData]) -> Result<()> {
        let mut appender = conn.appender("Coverage")
            .context("Failed to create Coverage appender")?;

        use crate::types::get_column_scale;
        let s_af  = get_column_scale("Coverage", "Aligned_fraction_percentage");
        let s_eaf = get_column_scale("Coverage", "Expected_aligned_fraction");
        let s_m   = get_column_scale("Coverage", "Coverage_mean");
        let s_med = get_column_scale("Coverage", "Coverage_median");
        let s_tm  = get_column_scale("Coverage", "Coverage_trimmed_mean");
        for p in presences {
            let cov_pct = (p.coverage_pct as f64 * s_af).round() as i32;
            let expected_af = (p.expected_aligned_fraction as f64 * s_eaf).round() as i32;
            let read_count = p.read_count as i64;
            let cov_mean = (p.coverage_mean as f64 * s_m).round() as i32;
            let cov_median = (p.coverage_median as f64 * s_med).round() as i32;
            let cov_trimmed_mean = (p.coverage_trimmed_mean as f64 * s_tm).round() as i32;
            let cov_sd = p.coverage_coefficient_of_variation.round() as i32;
            let cov_var = p.coverage_relative_coverage_roughness.round() as i32;
            appender.append_row(params![p.contig_id, sample_id, cov_pct, expected_af, read_count, cov_mean, cov_median, cov_trimmed_mean, cov_sd, cov_var])?;
        }

        appender.flush().context("Failed to flush Coverage appender")?;
        Ok(())
    }

    fn write_packaging(&self, conn: &Connection, sample_id: i64, packaging: &[PackagingData], circular: bool) -> Result<()> {
        let mut pkg_guard = self.next_packaging_id.lock()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut term_guard = self.next_terminus_id.lock()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut packaging_id = *pkg_guard;
        let mut terminus_id = *term_guard;

        let mut mechanism_appender = conn.appender("Phage_mechanisms")
            .context("Failed to create Phage_mechanisms appender")?;
        let mut termini_appender = conn.appender("Phage_termini")
            .context("Failed to create Phage_termini appender")?;

        for pkg in packaging {
                let contig_id = pkg.contig_id;
                let kept_terminus_positions: Vec<i32> = pkg.left_termini.iter()
                    .chain(pkg.right_termini.iter())
                    .filter(|t| t.passed_poisson_test && t.passed_clipping_test)
                    .map(|t| t.center_pos)
                    .collect();

                let genome_length = self.contig_id_to_length.get(&contig_id).copied().unwrap_or(0);

                let terminase_distance = calculate_terminase_distance_cached(
                    &self.terminase_cache, contig_id, &kept_terminus_positions, genome_length, circular,
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
                let s_tau = crate::types::get_column_scale("Phage_termini", "Tau");
                let s_cr = crate::types::get_column_scale("Phage_termini", "Clipped_ratio");
                for terminus in &pkg.left_termini {
                    let tau_int = (terminus.tau * s_tau).round() as i32;
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
                        terminus.passed_poisson_test,
                        terminus.expected_spc.round() as i64,
                        compact_pvalue(terminus.pvalue),
                        compact_pvalue(terminus.adjusted_pvalue),
                        terminus.passed_clipping_test,
                        terminus.sum_clippings as i64,
                        (terminus.clipped_ratio * s_cr).round() as i32,
                        terminus.expected_clippings.round() as i64
                    ])?;
                    terminus_id += 1;
                }

                // Insert right termini (status = "end")
                for terminus in &pkg.right_termini {
                    let tau_int = (terminus.tau * s_tau).round() as i32;
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
                        terminus.passed_poisson_test,
                        terminus.expected_spc.round() as i64,
                        compact_pvalue(terminus.pvalue),
                        compact_pvalue(terminus.adjusted_pvalue),
                        terminus.passed_clipping_test,
                        terminus.sum_clippings as i64,
                        (terminus.clipped_ratio * s_cr).round() as i32,
                        terminus.expected_clippings.round() as i64
                    ])?;
                    terminus_id += 1;
                }

                packaging_id += 1;
        }

        mechanism_appender.flush().context("Failed to flush Phage_mechanisms appender")?;
        termini_appender.flush().context("Failed to flush Phage_termini appender")?;

        *pkg_guard = packaging_id;
        *term_guard = terminus_id;
        Ok(())
    }

    fn write_misassembly(&self, conn: &Connection, sample_id: i64, data: &[MisassemblyData]) -> Result<()> {
        let mut appender = conn.appender("Misassembly")
            .context("Failed to create Misassembly appender")?;

        for d in data {
            appender.append_row(params![
                d.contig_id, sample_id,
                d.mismatches_count, d.deletions_count, d.insertions_count, d.clippings_count,
                d.collapse_bp, d.expansion_bp
            ])?;
        }

        appender.flush().context("Failed to flush Misassembly appender")?;
        Ok(())
    }

    fn write_microdiversity(&self, conn: &Connection, sample_id: i64, data: &[MicrodiversityData]) -> Result<()> {
        let mut appender = conn.appender("Microdiversity")
            .context("Failed to create Microdiversity appender")?;

        for d in data {
            appender.append_row(params![
                d.contig_id, sample_id,
                d.mismatches_count, d.deletions_count, d.insertions_count, d.clippings_count,
                d.microdiverse_bp_on_reference, d.microdiverse_bp_on_reads
            ])?;
        }

        appender.flush().context("Failed to flush Microdiversity appender")?;
        Ok(())
    }

    fn write_side_misassembly(&self, conn: &Connection, sample_id: i64, data: &[SideMisassemblyData]) -> Result<()> {
        let mut appender = conn.appender("Side_misassembly")
            .context("Failed to create Side_misassembly appender")?;

        for d in data {
            let misjoint = d.contig_end_misjoint_mates.map(|v| v as i64);
            appender.append_row(params![
                d.contig_id, sample_id,
                d.coverage_first_position as i64,
                d.contig_start_collapse_percentage, d.contig_start_collapse_bp, d.contig_start_expansion_bp,
                d.coverage_last_position as i64,
                d.contig_end_collapse_percentage, d.contig_end_collapse_bp, d.contig_end_expansion_bp,
                misjoint
            ])?;
        }

        appender.flush().context("Failed to flush Side_misassembly appender")?;
        Ok(())
    }

    fn write_topology(&self, conn: &Connection, sample_id: i64, data: &[TopologyData]) -> Result<()> {
        let mut appender = conn.appender("Topology")
            .context("Failed to create Topology appender")?;

        for d in data {
            let circ_reads = d.circularising_reads.map(|v| v as i64);
            let circ_pct = d.circularising_reads_percentage;
            let median_circ_len = d.median_circularising_len;
            let circ_inserts = d.circularising_inserts.map(|v| v as i64);
            let circ_dev = d.circularising_insert_size_deviation;
            appender.append_row(params![
                d.contig_id, sample_id,
                circ_reads, circ_pct, median_circ_len, circ_inserts, circ_dev
            ])?;
        }

        appender.flush().context("Failed to flush Topology appender")?;
        Ok(())
    }

    /// Write feature BLOBs to the Feature_blob table.
    ///
    /// Each tuple is (feature_name, contig_id, encoded_blob).
    /// Feature names are mapped to feature_ids via VARIABLES index.
    pub fn write_feature_blobs(&self, conn: &Connection, sample_id: i64, blobs: &[(String, i64, crate::blob::EncodedBlob)]) -> Result<()> {
        self.write_feature_blobs_timed(conn, sample_id, blobs, false, &mut SampleWriteTimings::default())
    }

    fn write_feature_blobs_timed(&self, conn: &Connection, sample_id: i64,
        blobs: &[(String, i64, crate::blob::EncodedBlob)], enable_timing: bool,
        timings: &mut SampleWriteTimings) -> Result<()> {
        if blobs.is_empty() {
            return Ok(());
        }

        let mut sizes: BTreeMap<String, BlobSize> = BTreeMap::new();
        let start = timed_start(enable_timing);
        let mut appender = conn.appender("Feature_blob")
            .context("Failed to create Feature_blob appender")?;
        let mut chunk_appender = conn.appender("Feature_blob_chunk")
            .context("Failed to create Feature_blob_chunk appender")?;

        for (feature_name, contig_id, encoded) in blobs {
            encoded.ensure_valid()?;
            if encoded.zoom.is_empty() {
                continue;
            }
            if let Some(fid) = feature_name_to_id(feature_name) {
                if enable_timing {
                    if let Some(size) = sizes.get_mut(feature_name) { size.record(encoded); }
                    else { let mut size = BlobSize::default(); size.record(encoded); sizes.insert(feature_name.clone(), size); }
                    timings.zoom_rows += 1;
                    timings.zoom_bytes += encoded.zoom.len() as u64;
                    timings.chunk_rows += encoded.chunks.iter().filter(|c| encoded.should_store_chunk(c)).count() as u64;
                    timings.chunk_bytes += encoded.chunks.iter().filter(|c| encoded.should_store_chunk(c)).map(|c| c.len() as u64).sum::<u64>();
                }
                appender.append_row(params![contig_id, sample_id, fid as i32, encoded.zoom.as_slice()])?;
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    if !encoded.should_store_chunk(chunk_data) { continue; }
                    chunk_appender.append_row(params![
                        contig_id, sample_id, fid as i32, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        timings.blob_append_ns = elapsed_ns(start);
        let start = timed_start(enable_timing);
        appender.flush().context("Failed to flush Feature_blob appender")?;
        timings.zoom_flush_ns = elapsed_ns(start);
        let start = timed_start(enable_timing);
        chunk_appender.flush().context("Failed to flush Feature_blob_chunk appender")?;
        timings.chunk_flush_ns = elapsed_ns(start);
        if enable_timing { self.record_blob_sizes("Feature_blob", sizes)?; self.flush_size_report(false)?; }
        Ok(())
    }

    /// Write one MAG-level coverage summary row, computed from raw
    /// per-position data (see `mag_blob::compute_mag_coverage_stats_from_raw`).
    pub fn write_mag_coverage(
        &self,
        conn: &Connection,
        mag_id: i64,
        sample_id: i64,
        s: &crate::mag_blob::MagCoverageStats,
    ) -> Result<()> {
        let mut appender = conn.appender("MAG_coverage")
            .context("Failed to create MAG_coverage appender")?;
        appender.append_row(params![
            mag_id,
            sample_id,
            s.mag_length as i64,
            s.aligned_fraction_pct_x10,
            s.expected_aligned_fraction_pct_x10,
            s.read_count as i64,
            s.coverage_mean_x10,
            s.coverage_median_x10,
            s.coverage_trimmed_mean_x10,
            s.coverage_cv_x1e6,
            s.coverage_roughness_x1e6,
        ])?;
        appender.flush().context("Failed to flush MAG_coverage appender")?;
        Ok(())
    }

    pub(crate) fn write_mag_coverages(&self, conn: &Connection, sample_id: i64,
        results: &[(String, i64, crate::mag_blob::MagCoverageStats)]) -> Result<()> {
        if results.is_empty() { return Ok(()); }
        let mut appender = conn.appender("MAG_coverage")?;
        for (_, mag_id, s) in results {
            appender.append_row(params![mag_id, sample_id, s.mag_length as i64,
                s.aligned_fraction_pct_x10, s.expected_aligned_fraction_pct_x10, s.read_count as i64,
                s.coverage_mean_x10, s.coverage_median_x10, s.coverage_trimmed_mean_x10,
                s.coverage_cv_x1e6, s.coverage_roughness_x1e6])?;
        }
        appender.flush()?;
        Ok(())
    }

    /// Persist a bounded directed BLAST batch in one transaction. Exact repeated
    /// query/subject rows retain the first identity value via INSERT OR IGNORE.
    /// Returns database-lock wait time and time spent writing under the lock.
    pub(crate) fn write_mag_blast_batch(
        &self,
        repeats: &[RepeatsData],
        hits: &[(String, String, i32, i32, i32, i32, f64)],
    ) -> Result<(std::time::Duration, std::time::Duration)> {
        use std::time::{Duration, Instant};
        if repeats.is_empty() && hits.is_empty() {
            return Ok((Duration::ZERO, Duration::ZERO));
        }
        let scale = crate::types::get_column_scale("Contig_directRepeats", "Pident");
        let mut seen = HashSet::new();
        let mut normalized = Vec::new();
        for (query, subject, p1, p2, p1p, p2p, identity) in hits {
            let (Some(qid), Some(sid)) = (self.get_contig_id(query), self.get_contig_id(subject)) else {
                continue;
            };
            // Preserve BLAST query/subject direction. Reciprocal rows belong to
            // different query tracks, so neither is mirrored or swapped.
            let key = (qid, sid, *p1, *p2, *p1p, *p2p);
            if seen.insert(key) {
                normalized.push((key, (identity * scale).round() as i32));
            }
        }
        let started = Instant::now();
        let mut conn = self.lock_conn()?;
        let waited = started.elapsed();
        let started = Instant::now();
        let tx = conn.transaction()?;
        if !repeats.is_empty() {
            let mut direct = tx.appender("Contig_directRepeats")?;
            let mut inverted = tx.appender("Contig_invertedRepeats")?;
            for repeat in repeats {
                if let Some(id) = self.get_contig_id(&repeat.contig_name) {
                    let appender = if repeat.is_direct { &mut direct } else { &mut inverted };
                    let repeat_scale = crate::types::get_column_scale(
                        if repeat.is_direct { "Contig_directRepeats" } else { "Contig_invertedRepeats" }, "Pident",
                    );
                    appender.append_row(params![id, repeat.position1, repeat.position2,
                        repeat.position1prime, repeat.position2prime,
                        (repeat.pident * repeat_scale).round() as i32])?;
                }
            }
            direct.flush()?;
            inverted.flush()?;
        }
        if !normalized.is_empty() {
            tx.execute_batch("CREATE TEMP TABLE IF NOT EXISTS Mag_blast_batch AS SELECT * FROM Contig_blast_hits WHERE false")?;
            {
                let mut appender = tx.appender("Mag_blast_batch")?;
                for ((c1, c2, p1, p2, p1p, p2p), identity) in normalized {
                    appender.append_row(params![c1, c2, p1, p2, p1p, p2p, identity])?;
                }
                appender.flush()?;
            }
            tx.execute_batch("INSERT OR IGNORE INTO Contig_blast_hits SELECT * FROM Mag_blast_batch;
                              DELETE FROM Mag_blast_batch")?;
        }
        tx.commit()?;
        Ok((waited, started.elapsed()))
    }

    /// Write repeats data to the database (both direct and inverted).
    /// This is called once during database creation (not per-sample).
    pub fn write_repeats(&self, repeats: &[RepeatsData]) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

        if repeats.is_empty() {
            // Repeat detection ran but found no repeats — set 0% for all contigs
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
                let pident_scale = if dup.is_direct {
                    crate::types::get_column_scale("Contig_directRepeats", "Pident")
                } else {
                    crate::types::get_column_scale("Contig_invertedRepeats", "Pident")
                };
                let pident_int = (dup.pident * pident_scale).round() as i32;

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

    pub fn compute_duplication_percentages(&self, repeats: &[RepeatsData]) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        self.update_duplication_percentages(&conn, repeats)
    }

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

            let dup_scale = crate::types::get_column_scale("Contig", "Duplication_percentage");
            let percentage = ((total_bp as f64 / contig_length as f64) * 100.0 * dup_scale).round() as i32;

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
        let reporting = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?.is_some();
        let mut sizes: BTreeMap<String, BlobSize> = BTreeMap::new();
        use crate::blob::encode_contig_dense_blob;

        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        let gc_content_feature_id = feature_name_to_id("gc_content")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: gc_content"))?;
        let gc_window_size = DEFAULT_GC_CONTENT_WINDOW_SIZE as u32;

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                if data.gc_values.is_empty() {
                    continue;
                }

                // Convert u8 values to i32 for BLOB encoding
                let values: Vec<i32> = data.gc_values.iter().map(|&v| v as i32).collect();

                // Estimate contig length from number of windows
                let contig_length = (values.len() as u32) * gc_window_size;

                // Encode as dense BLOB (array index is window number)
                debug_assert_eq!(crate::types::get_encoding("gc_content"), crate::types::Encoding::Dense);
                let encoded = encode_contig_dense_blob(
                    &values,
                    crate::types::get_value_scale("gc_content"),
                    contig_length,
                    gc_window_size,
                );

                appender.append_row(params![
                    contig_id,
                    gc_content_feature_id,
                    { encoded.ensure_valid()?; encoded.zoom.as_slice() }
                ])?;
                if reporting { sizes.entry("gc_content".to_owned()).or_default().record(&encoded); }
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    if !encoded.should_store_chunk(chunk_data) { continue; }
                    chunk_appender.append_row(params![
                        contig_id, gc_content_feature_id, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;

        if reporting { self.record_blob_sizes("Contig_blob", sizes)?; self.flush_size_report(false)?; }
        Ok(())
    }

    /// Write GC skew data to Contig_blob table.
    /// Stores raw window-level GC skew values as dense BLOBs (one value per 1000bp window).
    /// This is called once during database creation (not per-sample).
    pub fn write_gc_skew(&self, gc_data: &[GCContentData]) -> Result<()> {
        let reporting = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?.is_some();
        let mut sizes: BTreeMap<String, BlobSize> = BTreeMap::new();
        use crate::blob::encode_contig_dense_blob;

        if gc_data.is_empty() {
            return Ok(());
        }

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        let gc_skew_feature_id = feature_name_to_id("gc_skew")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: gc_skew"))?;
        let gc_skew_window_size = DEFAULT_GC_SKEW_WINDOW_SIZE as u32;

        for data in gc_data {
            if let Some(&contig_id) = self.contig_name_to_id.get(&data.contig_name) {
                if data.skew_values.is_empty() {
                    continue;
                }

                // Convert i16 values to i32 for BLOB encoding
                let values: Vec<i32> = data.skew_values.iter().map(|&v| v as i32).collect();

                // Estimate contig length from number of windows
                let contig_length = (values.len() as u32) * gc_skew_window_size;

                // Encode as dense BLOB (array index is window number)
                // Values are stored as integer ×100 (e.g. -40 = -0.40), use Times100 scale
                debug_assert_eq!(crate::types::get_encoding("gc_skew"), crate::types::Encoding::Dense);
                let encoded = encode_contig_dense_blob(
                    &values,
                    crate::types::get_value_scale("gc_skew"),
                    contig_length,
                    gc_skew_window_size,
                );

                appender.append_row(params![
                    contig_id,
                    gc_skew_feature_id,
                    { encoded.ensure_valid()?; encoded.zoom.as_slice() }
                ])?;
                if reporting { sizes.entry("gc_skew".to_owned()).or_default().record(&encoded); }
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    if !encoded.should_store_chunk(chunk_data) { continue; }
                    chunk_appender.append_row(params![
                        contig_id, gc_skew_feature_id, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;

        if reporting { self.record_blob_sizes("Contig_blob", sizes)?; self.flush_size_report(false)?; }
        Ok(())
    }

    /// Build the four materialized repeat tables, encode them into Contig_blob,
    /// and drop the materialized intermediates. Selected blobs are replaced
    /// atomically; repeated calls are safe.
    /// First upgrades regenerate every existing repeat blob. Called after
    /// `write_repeats` but before the per-sample loop, so MAG-scale aggregation
    /// can read repeat data from `Contig_blob` like everything else.
    pub fn convert_repeat_blobs(&self, new_contig_ids: &[i64]) -> Result<()> {
        if new_contig_ids.is_empty() { return Ok(()); }
        let mut conn = self.lock_conn()?;
        let feature_ids = ["direct_repeat_count","inverted_repeat_count","direct_repeat_identity","inverted_repeat_identity"]
            .iter().map(|name| feature_name_to_id(name).expect("Repeat feature definition")).collect::<Vec<_>>();
        let has_repeats: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM Contig_directRepeats LIMIT 1) OR EXISTS(SELECT 1 FROM Contig_invertedRepeats LIMIT 1)", [], |r| r.get(0))?;
        let has_blobs: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM Contig_blob WHERE Feature_id IN (?,?,?,?) LIMIT 1)",params![feature_ids[0],feature_ids[1],feature_ids[2],feature_ids[3]],|r|r.get(0))?;
        if !has_repeats && !has_blobs {return Ok(());}
        let current: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM Database_metadata WHERE Key='Repeat_blob_encoding' AND Value='interval-v2')", [], |r|r.get(0))?;
        let tx = conn.transaction()?;
        tx.execute("CREATE TEMP TABLE Repeat_conversion_ids (Contig_id BIGINT)", [])?;
        {
            let mut ids=tx.appender("Repeat_conversion_ids")?;
            for id in new_contig_ids {ids.append_row(params![id])?;}
            ids.flush()?;
        }
        if !current && has_blobs {
            // Upgrade every old blob and every stored repeat, not just new references.
            tx.execute_batch("INSERT INTO Repeat_conversion_ids
                SELECT Contig_id FROM Contig_directRepeats UNION SELECT Contig_id FROM Contig_invertedRepeats")?;
            tx.execute("INSERT INTO Repeat_conversion_ids SELECT DISTINCT Contig_id FROM Contig_blob WHERE Feature_id IN (?,?,?,?)",
                params![feature_ids[0],feature_ids[1],feature_ids[2],feature_ids[3]])?;
        }
        Self::build_materialized_repeat_tables(&tx)?;
        for id in feature_ids {
            tx.execute("DELETE FROM Contig_blob_chunk WHERE Feature_id=? AND Contig_id IN (SELECT Contig_id FROM Repeat_conversion_ids)",params![id])?;
            tx.execute("DELETE FROM Contig_blob WHERE Feature_id=? AND Contig_id IN (SELECT Contig_id FROM Repeat_conversion_ids)",params![id])?;
        }
        let sizes=self.write_repeat_features_to_blob_static(&tx)?;
        for table in ["Contig_direct_repeat_count","Contig_inverted_repeat_count","Contig_direct_repeat_identity","Contig_inverted_repeat_identity","Repeat_conversion_ids"] {
            tx.execute(&format!("DROP TABLE {}",table),[])?;
        }
        // Blobs, temporary state, and the version marker commit together. Any
        // failure rolls back through Transaction::drop, preserving old blobs.
        tx.commit()?;
        self.record_blob_sizes("Contig_blob",sizes)?;
        self.flush_size_report(false)?;
        Ok(())
    }

    /// Build the four CTAS aggregates over `Contig_directRepeats` /
    /// `Contig_invertedRepeats`. Used by `convert_repeat_blobs`.
    /// Selected IDs live in a temporary table; SEMI JOIN preserves IN semantics
    /// even if the caller supplies duplicate IDs.
    fn build_materialized_repeat_tables(conn: &Connection) -> Result<()> {
        let wh = " SEMI JOIN Repeat_conversion_ids USING (Contig_id)";

        conn.execute(
            &format!("CREATE TEMP TABLE Contig_direct_repeat_count AS
            WITH events AS (
                SELECT Contig_id, LEAST(Position1, Position2) AS pos, 1 AS delta
                FROM Contig_directRepeats{wh}
                UNION ALL
                SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos, -1 AS delta
                FROM Contig_directRepeats{wh}
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
            FROM running WHERE next_pos IS NOT NULL AND count > 0"),
            [],
        )
        .context("Failed to create Contig_direct_repeat_count table")?;

        conn.execute(
            &format!("CREATE TEMP TABLE Contig_inverted_repeat_count AS
            WITH events AS (
                SELECT Contig_id, LEAST(Position1, Position2) AS pos, 1 AS delta
                FROM Contig_invertedRepeats{wh}
                UNION ALL
                SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos, -1 AS delta
                FROM Contig_invertedRepeats{wh}
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
            FROM running WHERE next_pos IS NOT NULL AND count > 0"),
            [],
        )
        .context("Failed to create Contig_inverted_repeat_count table")?;

        let pident_dr = crate::types::get_column_scale("Contig_directRepeats", "Pident");
        let pident_ir = crate::types::get_column_scale("Contig_invertedRepeats", "Pident");

        conn.execute(
            &format!("CREATE TEMP TABLE Contig_direct_repeat_identity AS
            WITH boundaries AS (
                SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_directRepeats{wh}
                UNION
                SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_directRepeats{wh}
            ),
            segments AS (
                SELECT Contig_id, pos AS seg_start,
                       LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
                FROM boundaries
            )
            SELECT s.Contig_id, s.seg_start AS First_position, s.seg_end AS Last_position,
                   MAX(r.Pident) / {pident_dr:.1} AS Value,
                   arg_max(LEAST(r.Position1prime, r.Position2prime),    r.Pident) AS Partner_start,
                   arg_max(GREATEST(r.Position1prime, r.Position2prime), r.Pident) AS Partner_end
            FROM segments s
            JOIN Contig_directRepeats r ON r.Contig_id = s.Contig_id
              AND LEAST(r.Position1, r.Position2) <= s.seg_end
              AND GREATEST(r.Position1, r.Position2) >= s.seg_start
            WHERE s.seg_end IS NOT NULL
            GROUP BY s.Contig_id, s.seg_start, s.seg_end"),
            [],
        )
        .context("Failed to create Contig_direct_repeat_identity table")?;

        conn.execute(
            &format!("CREATE TEMP TABLE Contig_inverted_repeat_identity AS
            WITH boundaries AS (
                SELECT Contig_id, LEAST(Position1, Position2) AS pos FROM Contig_invertedRepeats{wh}
                UNION
                SELECT Contig_id, GREATEST(Position1, Position2) + 1 AS pos FROM Contig_invertedRepeats{wh}
            ),
            segments AS (
                SELECT Contig_id, pos AS seg_start,
                       LEAD(pos) OVER (PARTITION BY Contig_id ORDER BY pos) - 1 AS seg_end
                FROM boundaries
            )
            SELECT s.Contig_id, s.seg_start AS First_position, s.seg_end AS Last_position,
                   MAX(r.Pident) / {pident_ir:.1} AS Value,
                   arg_max(LEAST(r.Position1prime, r.Position2prime),    r.Pident) AS Partner_start,
                   arg_max(GREATEST(r.Position1prime, r.Position2prime), r.Pident) AS Partner_end
            FROM segments s
            JOIN Contig_invertedRepeats r ON r.Contig_id = s.Contig_id
              AND LEAST(r.Position1, r.Position2) <= s.seg_end
              AND GREATEST(r.Position1, r.Position2) >= s.seg_start
            WHERE s.seg_end IS NOT NULL
            GROUP BY s.Contig_id, s.seg_start, s.seg_end"),
            [],
        )
        .context("Failed to create Contig_inverted_repeat_identity table")?;

        Ok(())
    }

    /// Write repeat count and identity features to Contig_blob table from materialized tables.
    /// Intervals cover all overlapping zoom bins and carry chunk-local plateau
    /// endpoints plus explicit gap anchors. Identity metadata retains partner bounds.
    fn write_repeat_features_to_blob_static(&self, conn: &Connection) -> Result<BTreeMap<String,BlobSize>> {
        let reporting = self.size_report.lock().map_err(|e| anyhow::anyhow!("{e}"))?.is_some();
        let mut sizes: BTreeMap<String, BlobSize> = BTreeMap::new();
        use crate::blob::MetadataFlags;
        use crate::types::{get_encoding, get_value_scale, Encoding};

        // Resolve feature IDs from VARIABLES so they stay in sync with types.rs
        // and can't collide with other features like hit_count_within_mag (id 5)
        // and hit_identity_within_mag (id 6).
        let direct_repeat_count_id = feature_name_to_id("direct_repeat_count")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: direct_repeat_count"))?;
        let inverted_repeat_count_id = feature_name_to_id("inverted_repeat_count")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: inverted_repeat_count"))?;
        let direct_repeat_identity_id = feature_name_to_id("direct_repeat_identity")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: direct_repeat_identity"))?;
        let inverted_repeat_identity_id = feature_name_to_id("inverted_repeat_identity")
            .ok_or_else(|| anyhow::anyhow!("Missing VARIABLES entry: inverted_repeat_identity"))?;

        // (feature_name, table_name, feature_id, is_identity)
        let tables: Vec<(&str, &str, i16, bool)> = vec![
            ("direct_repeat_count", "Contig_direct_repeat_count", direct_repeat_count_id, false),
            ("inverted_repeat_count", "Contig_inverted_repeat_count", inverted_repeat_count_id, false),
            ("direct_repeat_identity", "Contig_direct_repeat_identity", direct_repeat_identity_id, true),
            ("inverted_repeat_identity", "Contig_inverted_repeat_identity", inverted_repeat_identity_id, true),
        ];

        // The writer already maintains this cache during import/open.
        let mut appender = conn.appender("Contig_blob")
            .context("Failed to create Contig_blob appender")?;
        let mut chunk_appender = conn.appender("Contig_blob_chunk")
            .context("Failed to create Contig_blob_chunk appender")?;

        for (feature_name, table_name, feature_id, is_identity) in tables {
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

            debug_assert_eq!(get_encoding(feature_name), Encoding::Sparse);
            let scale = get_value_scale(feature_name);

            // Read RLE segments from the materialized table.
            // Identity tables carry Partner_start / Partner_end columns (positions within
            // the contig of the partner copy of the max-Pident repeat covering the segment).
            let query = if is_identity {
                format!("SELECT Contig_id, First_position, Last_position, Value, Partner_start, Partner_end FROM {} ORDER BY Contig_id, First_position", table_name)
            } else {
                format!("SELECT Contig_id, First_position, Last_position, Value FROM {} ORDER BY Contig_id, First_position", table_name)
            };
            let mut stmt = conn.prepare(&query)?;
            let rows = stmt.query_map([], |row| {
                Ok((
                    row.get::<_, i64>(0)?,    // Contig_id
                    row.get::<_, i32>(1)?,    // First_position
                    row.get::<_, i32>(2)?,    // Last_position
                    row.get::<_, f64>(3)?,    // Value
                    if is_identity { row.get::<_, Option<i32>>(4)? } else { None }, // Partner_start
                    if is_identity { row.get::<_, Option<i32>>(5)? } else { None }, // Partner_end
                ))
            })?;

            // (first_pos, last_pos, value, partner_start, partner_end)
            let mut contig_segments: HashMap<i64, Vec<(i32, i32, f64, Option<i32>, Option<i32>)>> = HashMap::new();
            for row in rows {
                let (contig_id, first_pos, last_pos, value, partner_start, partner_end) = row?;
                contig_segments.entry(contig_id).or_default().push((first_pos, last_pos, value, partner_start, partner_end));
            }

            let flags = MetadataFlags {
                sparse: true,
                has_partner: is_identity,
                ..MetadataFlags::default()
            };

            for (contig_id, segments) in &contig_segments {
                let contig_length = self.contig_id_to_length.get(contig_id).copied().unwrap_or(0).max(0) as u32;
                if contig_length == 0 || segments.is_empty() {
                    continue;
                }

                let intervals = segments.iter().filter_map(|&(first,last,value,partner_start,partner_end)| {
                    if first < 1 || last < first || first as u32 > contig_length {return None;}
                    Some(crate::blob::SparseInterval {
                        start:first as u32-1,
                        end:(last as u32-1).min(contig_length-1),
                        value:if is_identity {(value*scale.multiplier()).round() as i32} else {value as i32},
                        partner_start:partner_start.map(|position|position-1),
                        partner_end:partner_end.map(|position|position-1),
                    })
                }).collect::<Vec<_>>();
                if intervals.is_empty() {continue;}
                let encoded = crate::blob::encode_sparse_intervals(&intervals,flags,scale,contig_length);
                appender.append_row(params![*contig_id, feature_id, {encoded.ensure_valid()?;encoded.zoom.as_slice()}])?;
                if reporting { if let Some(size) = sizes.get_mut(feature_name) { size.record(&encoded); } else { let mut size = BlobSize::default(); size.record(&encoded); sizes.insert(feature_name.to_string(), size); } }
                for (chunk_idx, chunk_data) in encoded.chunks.iter().enumerate() {
                    if !encoded.should_store_chunk(chunk_data) { continue; }
                    chunk_appender.append_row(params![
                        *contig_id, feature_id, chunk_idx as i16, chunk_data.as_slice()
                    ])?;
                }
            }
        }

        appender.flush().context("Failed to flush Contig_blob appender")?;
        chunk_appender.flush().context("Failed to flush Contig_blob_chunk appender")?;
        conn.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Repeat_blob_encoding','interval-v2')", [])?;
        Ok(sizes)
    }

    /// Update Contig table with GC statistics (average, sd) and GC skew stats.
    /// Called after computing GC content and GC skew for all contigs.
    pub fn update_contig_gc_stats(&self, gc_data: &[GCContentData]) -> Result<()> {
        self.stage_contig_gc_stats(gc_data)?;
        self.flush_contig_gc_stats()
    }

    pub(crate) fn stage_contig_gc_stats(&self, gc_data: &[GCContentData]) -> Result<()> {
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

        let s_gc  = crate::types::get_column_scale("Contig", "GC_mean") as f32;
        let s_sd = crate::types::get_column_scale("Contig", "GC_sd") as f32;
        let s_amp = crate::types::get_column_scale("Contig", "GC_skew_amplitude") as f32;
        let s_pos = crate::types::get_column_scale("Contig", "Positive_GC_skew_windows_percentage") as f32;

        conn.execute(
            "CREATE TEMPORARY TABLE IF NOT EXISTS tmp_gc (
                contig_id BIGINT, gc_mean INTEGER, gc_sd INTEGER, skew_amp INTEGER, pos_pct INTEGER
            )", [],
        ).context("Failed to create tmp_gc table")?;

        let mut app = conn.appender("tmp_gc")
            .context("Failed to create tmp_gc appender")?;
        for data in gc_data {
            app.append_row(params![
                self.contig_name_to_id.get(&data.contig_name).copied(),
                (data.stats.average * s_gc).round() as i32,
                (data.stats.sd * s_sd).round() as i32,
                (data.skew_stats.amplitude * s_amp).round() as i32,
                (data.skew_stats.percent_positive * s_pos).round() as i32,
            ])?;
        }
        app.flush().context("Failed to flush tmp_gc appender")?;

        Ok(())
    }

    pub(crate) fn flush_contig_gc_stats(&self) -> Result<()> {
        let conn = self.lock_conn()?;
        conn.execute(
            "UPDATE Contig SET
                GC_mean = t.gc_mean,
                GC_sd = t.gc_sd,
                GC_skew_amplitude = t.skew_amp,
                Positive_GC_skew_windows_percentage = t.pos_pct
            FROM tmp_gc t WHERE Contig.Contig_id = t.contig_id", [],
        ).context("Failed to bulk-update Contig GC stats")?;

        conn.execute("DROP TABLE tmp_gc", [])
            .context("Failed to clear tmp_gc")?;

        Ok(())
    }

    /// Write MAG rows and MAG_contigs_association rows.
    ///
    /// `mags` is a slice of (mag_name, member_contig_names). Aggregates are computed
    /// by length-weighted mean over the member contigs (means for rates, SUM for
    /// counts). Run AFTER per-contig GC stats and Duplication_percentage have been
    /// written. Associations are persisted in descending contig-length order within
    /// each MAG — downstream visualization relies on this ordering.
    pub fn write_mags(&self, mags: &[(String, Vec<String>)]) -> Result<()> {
        if mags.is_empty() {
            return Ok(());
        }
        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;

        // Support extend mode: new MAG ids start past the existing max (0 if empty).
        let existing_max: i64 = conn
            .query_row("SELECT COALESCE(MAX(MAG_id), 0) FROM MAG", [], |r| r.get(0))
            .unwrap_or(0);
        let existing_names: HashSet<String> = {
            let mut stmt = conn.prepare("SELECT MAG_name FROM MAG")?;
            let rows = stmt.query_map([], |r| r.get::<_, String>(0))?;
            rows.collect::<std::result::Result<HashSet<_>, _>>()?
        };

        let tx = conn.unchecked_transaction()?;

        // Stage numeric IDs once, then resolve metadata with one database join.
        // Member_order preserves the existing stable tie order for equal lengths.
        tx.execute("CREATE TEMP TABLE Mag_member_input (Mag_id BIGINT, Contig_id BIGINT, Member_order BIGINT)", [])?;
        let mut selected = Vec::new();
        {
            let mut input = tx.appender("Mag_member_input")?;
            let mut next_id = existing_max;
            for (name, contig_names) in mags {
                if existing_names.contains(name) { continue; }
                next_id += 1;
                selected.push((name, next_id));
                for (order, contig_name) in contig_names.iter().enumerate() {
                    let contig_id = self.contig_name_to_id.get(contig_name)
                        .with_context(|| format!("MAG '{}' references unknown contig '{}'", name, contig_name))?;
                    input.append_row(params![next_id, contig_id, order as i64])?;
                }
            }
            input.flush()?;
        }
        {
        let mut stmt = tx.prepare(
            "SELECT i.Mag_id, c.Contig_id, c.Contig_length, c.Duplication_percentage,
                    c.GC_mean, c.GC_sd, c.GC_skew_amplitude, c.Positive_GC_skew_windows_percentage
             FROM Mag_member_input i JOIN Contig c USING (Contig_id)
             ORDER BY i.Mag_id, i.Member_order")?;
        let mut rows = stmt.query_map([], |r| Ok((r.get::<_, i64>(0)?, (
            r.get::<_, i64>(1)?, r.get::<_, i64>(2)?, r.get::<_, Option<i32>>(3)?,
            r.get::<_, Option<i32>>(4)?, r.get::<_, Option<i32>>(5)?,
            r.get::<_, Option<i32>>(6)?, r.get::<_, Option<i32>>(7)?))))?.peekable();
        for (mag_name, mag_id) in selected {
            let mut members = Vec::new();
            while let Some(row) = rows.peek() {
                match row {
                    Ok((id, _)) if *id != mag_id => break,
                    _ => { let (_, member) = rows.next().unwrap()?; members.push(member); }
                }
            }
            // Sort by length DESC; stable so ties break by insertion order.
            members.sort_by(|a, b| b.1.cmp(&a.1));

            let total_length: i64 = members.iter().map(|m| m.1).sum();
            let n_contigs = members.len() as i64;

            // N50: length of the shortest contig whose inclusion brings the cumulative
            // sum (sorted DESC) to ≥ 50 % of total assembly length.
            let n50: i64 = {
                let half = (total_length + 1) / 2;
                let mut cumsum: i64 = 0;
                let mut n50_val: i64 = members.last().map(|m| m.1).unwrap_or(0);
                for m in &members {
                    cumsum += m.1;
                    if cumsum >= half {
                        n50_val = m.1;
                        break;
                    }
                }
                n50_val
            };

            // Length-weighted mean helper over (value, length); returns None if all missing or length=0.
            let lw_mean = |values: &[(Option<i32>, i64)]| -> Option<f64> {
                let total_len: i64 = values.iter().filter(|(v, _)| v.is_some()).map(|(_, l)| *l).sum();
                if total_len == 0 {
                    return None;
                }
                let num: f64 = values
                    .iter()
                    .filter_map(|(v, l)| v.map(|vv| vv as f64 * *l as f64))
                    .sum();
                Some(num / total_len as f64)
            };

            let dup: Vec<(Option<i32>, i64)> = members.iter().map(|m| (m.2, m.1)).collect();
            let gc_mean_v: Vec<(Option<i32>, i64)> = members.iter().map(|m| (m.3, m.1)).collect();
            let gc_sd_v: Vec<(Option<i32>, i64)> = members.iter().map(|m| (m.4, m.1)).collect();
            let gc_amp_v: Vec<(Option<i32>, i64)> = members.iter().map(|m| (m.5, m.1)).collect();
            let gc_pos_v: Vec<(Option<i32>, i64)> = members.iter().map(|m| (m.6, m.1)).collect();

            let dup_mag = lw_mean(&dup).map(|v| v.round() as i32);
            let gc_mean_mag = lw_mean(&gc_mean_v).map(|v| v.round() as i32);
            let gc_sd_mag = lw_mean(&gc_sd_v).map(|v| v.round() as i32);
            let gc_amp_mag = lw_mean(&gc_amp_v).map(|v| v.round() as i32);
            let gc_pos_mag = lw_mean(&gc_pos_v).map(|v| v.round() as i32);

            {
                let mut mag = tx.appender("MAG")?;
                mag.append_row(params![mag_id, mag_name, total_length, n_contigs, n50,
                    dup_mag, gc_mean_mag, gc_sd_mag, gc_amp_mag, gc_pos_mag,
                    None::<i32>]) // Number_of_samples is populated during finalization.
                    .with_context(|| format!("Failed to append MAG '{}'", mag_name))?;
                mag.flush()?;
            }

            // Cumulative offset under longest-first ordering.
            let mut cumulative_offset: i64 = 0;
            let mut associations = tx.appender("MAG_contigs_association")?;
            for (contig_id, contig_length, _, _, _, _, _) in &members {
                associations.append_row(params![mag_id, contig_id, cumulative_offset])?;
                cumulative_offset += *contig_length;
            }
            associations.flush()?;
        }
        }
        tx.execute("DROP TABLE Mag_member_input", [])?;

        tx.commit().context("Failed to commit MAG rows")?;
        Ok(())
    }

    /// Derive within-MAG and external hit curves on shared Genome subplots.
    pub fn write_intercontig_hit_features(&self) -> Result<()> {
        if self.is_mag_mode {
            self.write_hit_features(true)?;
        }
        self.write_hit_features(false)?;
        Ok(())
    }

    fn write_hit_features(&self, within_mag: bool) -> Result<()> {
        use crate::blob::{encode_sparse_blob, EncodedBlob, EventMeta, MetadataFlags};
        use rayon::prelude::*;
        type Hits = Vec<(i64,i32,i32,i32)>;
        let count_name = if within_mag {"hit_count_within_mag"} else {"hit_count"};
        let identity_name = if within_mag {"hit_identity_within_mag"} else {"hit_identity"};
        let count_id=feature_name_to_id(count_name).context("Missing hit count feature")?;
        let identity_id=feature_name_to_id(identity_name).context("Missing hit identity feature")?;
        let reporting=self.size_report.lock().map_err(|e|anyhow::anyhow!("{e}"))?.is_some();
        let mut sizes:BTreeMap<String,BlobSize>=BTreeMap::new();
        // Reuse imported lengths; no second whole-reference metadata map.
        let lengths=&self.contig_id_to_length;
        let raw_scale=crate::types::get_column_scale("Contig_directRepeats","Pident");
        let feature_scale=crate::types::get_value_scale(identity_name).multiplier();
        let encode = |(cid,hits):(i64,Hits)| -> Result<Option<(i64,EncodedBlob,EncodedBlob)>> {
            let clen=lengths.get(&cid).copied().unwrap_or(0).max(0) as usize;
            if clen==0 {return Ok(None);}
                let mut count_diff = vec![0i32; clen + 1];
                let mut max_ident = vec![(0i32, 0i64); clen];

                for &(partner, a, b, pident) in &hits {
                    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
                    let lo = (lo as i64).max(1);
                    let hi = (hi as i64).min(clen as i64);
                    if lo > hi {continue;}
                    let (lo,hi)=((lo-1) as usize,(hi-1) as usize);

                    count_diff[lo] += 1;
                    if hi + 1 <= clen { count_diff[hi + 1] -= 1; }

                    for pos in lo..=hi.min(clen - 1) {
                        if pident > max_ident[pos].0 || (pident == max_ident[pos].0 && partner < max_ident[pos].1) {
                            max_ident[pos] = (pident, partner);
                        }
                    }
                }

                let mut running = 0i32;
                let mut pos_vec = Vec::new();
                let mut count_vec = Vec::new();
                let mut ident_vec = Vec::new();
                let mut meta_vec = Vec::new();
                for i in 0..clen {
                    running += count_diff[i];
                    if running > 0 {
                        pos_vec.push(i as u32);
                        count_vec.push(running);
                        ident_vec.push((max_ident[i].0 as f64 / raw_scale * feature_scale).round() as i32);
                        meta_vec.push(EventMeta { partner: Some(max_ident[i].1 as i32), ..Default::default() });
                    }
                }


            if pos_vec.is_empty() {return Ok(None);}
            let count=encode_sparse_blob(&pos_vec,&count_vec,None,MetadataFlags::default(),crate::types::get_value_scale(count_name),clen as u32);
            count.ensure_valid()?;
            let identity=encode_sparse_blob(&pos_vec,&ident_vec,Some(&meta_vec),MetadataFlags {has_partner:true,..Default::default()},crate::types::get_value_scale(identity_name),clen as u32);
            identity.ensure_valid()?;
            Ok(Some((cid,count,identity)))
        };
        let conn=self.lock_conn()?;
        for id in [count_id,identity_id] {
            conn.execute("DELETE FROM Contig_blob_chunk WHERE Feature_id=?",params![id])?;
            conn.execute("DELETE FROM Contig_blob WHERE Feature_id=?",params![id])?;
        }
        let same_mag="EXISTS (SELECT 1 FROM MAG_contigs_association a JOIN MAG_contigs_association b ON a.MAG_id=b.MAG_id WHERE a.Contig_id=h.Contig_id_1 AND b.Contig_id=h.Contig_id_2)";
        let predicate=if !self.is_mag_mode {"TRUE".to_owned()} else if within_mag {same_mag.to_owned()} else {format!("NOT ({same_mag})")};
        // DuckDB manages/spills the sort; Rust consumes rows one query contig at a time.
        let sql=format!("SELECT Contig_id_1,Contig_id_2,Position1,Position2,Pident FROM Contig_blast_hits h WHERE Contig_id_1<>Contig_id_2 AND {predicate} ORDER BY Contig_id_1");
        let mut stmt=conn.prepare(&sql)?;
        let rows=stmt.query_map([],|r|Ok((r.get::<_,i64>(0)?,r.get::<_,i64>(1)?,r.get::<_,i32>(2)?,r.get::<_,i32>(3)?,r.get::<_,i32>(4)?)))?;
        let mut appender=conn.appender("Contig_blob")?;
        let mut chunk_appender=conn.appender("Contig_blob_chunk")?;
        let mut flush = |batch:&mut Vec<(i64,Hits)>| -> Result<()> {
            let encoded=std::mem::take(batch).into_par_iter().map(&encode).collect::<Result<Vec<_>>>()?;
            for (cid,count,identity) in encoded.into_iter().flatten() {
                for (name,id,blob) in [(count_name,count_id,count),(identity_name,identity_id,identity)] {
                    appender.append_row(params![cid,id,blob.zoom.as_slice()])?;
                    if reporting {sizes.entry(name.to_owned()).or_default().record(&blob);}
                    for (index,chunk) in blob.chunks.iter().enumerate() {
                        if blob.should_store_chunk(chunk) {chunk_appender.append_row(params![cid,id,index as i16,chunk.as_slice()])?;}
                    }
                }
            }
            appender.flush()?;
            chunk_appender.flush()?;
            Ok(())
        };
        let estimate = |cid:i64,hits:&Hits| -> usize {
            (lengths.get(&cid).copied().unwrap_or(0).max(0) as usize)
                .saturating_mul(std::mem::size_of::<EventMeta>()+64)
                .saturating_add(hits.capacity().saturating_mul(std::mem::size_of::<(i64,i32,i32,i32)>()))
        };
        let budget=64*1024*1024;
        let workers=rayon::current_num_threads().clamp(1,32);
        let mut batch=Vec::new(); let mut bytes=0usize;
        let mut current=None; let mut hits=Hits::new();
        for row in rows {
            let (cid,partner,start,end,identity)=row?;
            if current!=Some(cid) {
                if let Some(previous)=current {
                    bytes=bytes.saturating_add(estimate(previous,&hits));
                    batch.push((previous,std::mem::take(&mut hits)));
                    if bytes>=budget || batch.len()>=workers {flush(&mut batch)?;bytes=0;}
                }
                current=Some(cid);
            }
            hits.push((partner,start,end,identity));
        }
        if let Some(cid)=current {batch.push((cid,hits));}
        flush(&mut batch)?;
        drop(flush);
        if reporting {self.record_blob_sizes("Contig_blob",sizes)?;self.flush_size_report(false)?;}
        Ok(())
    }

    /// Finalize the database after all samples are processed.
    /// Write initial metadata rows into the Database_metadata table.
    pub fn write_metadata(
        &self,
        modules: &[String],
        min_aligned_fraction: f64,
        min_coverage_depth: f64,
        bar_ratio: f64,
        variation_percentage: f64,
        min_occurrences: u32,
        view_mode: &str,
        blast: bool,
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

        let zoom_bins: String = crate::blob::ZOOM_BIN_SIZES.iter()
            .map(|b| b.to_string()).collect::<Vec<_>>().join(",");
        let contig_zoom_bins: String = crate::blob::CONTIG_ZOOM_BIN_SIZES.iter()
            .map(|b| b.to_string()).collect::<Vec<_>>().join(",");

        let rows: Vec<(&str, String)> = vec![
            ("Date_of_creation", now.clone()),
            ("Date_of_last_modification", now),
            ("Tool_version_used_for_creation", version.clone()),
            ("Tool_version_used_for_last_modification", version),
            ("Modules", modules_str),
            ("Min_aligned_fraction", min_aligned_fraction.to_string()),
            ("Min_coverage_depth", min_coverage_depth.to_string()),
            ("Coverage_percentage", bar_ratio.to_string()),
            ("Variation_percentage", variation_percentage.to_string()),
            ("Min_occurrences", min_occurrences.to_string()),
            ("View_mode", view_mode.to_string()),
            ("Blast_enabled", blast.to_string()),
            ("Chunk_size", crate::blob::CHUNK_SIZE.to_string()),
            ("Zoom_bin_sizes", zoom_bins),
            ("Contig_zoom_bin_sizes", contig_zoom_bins),
            ("GC_content_window_size", crate::gc_content::DEFAULT_GC_CONTENT_WINDOW_SIZE.to_string()),
            ("GC_skew_window_size", crate::gc_content::DEFAULT_GC_SKEW_WINDOW_SIZE.to_string()),
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

    pub fn write_compression_metadata(
        &self,
        min_aligned_fraction: f64,
        min_coverage_depth: f64,
        bar_ratio: f64,
        variation_percentage: f64,
        min_occurrences: u32,
    ) -> Result<()> {
        let conn = self.lock_conn()?;
        let keys = [
            ("Min_aligned_fraction", min_aligned_fraction.to_string()),
            ("Min_coverage_depth", min_coverage_depth.to_string()),
            ("Coverage_percentage", bar_ratio.to_string()),
            ("Variation_percentage", variation_percentage.to_string()),
            ("Min_occurrences", min_occurrences.to_string()),
        ];
        for (key, value) in &keys {
            let exists: bool = conn.query_row(
                "SELECT COUNT(*) > 0 FROM Database_metadata WHERE Key = ?",
                params![key], |row| row.get(0),
            ).unwrap_or(false);
            if !exists {
                conn.execute(
                    "INSERT INTO Database_metadata (Key, Value) VALUES (?, ?)",
                    params![key, value],
                )?;
            }
        }
        Ok(())
    }

    pub fn finalize(self) -> Result<FinalizeTimings> {
        let has_bam = self.has_bam;
        let is_mag_mode = self.is_mag_mode;
        let has_annotations = self.has_annotations;
        let report = self.size_report.into_inner().map_err(|e| anyhow::anyhow!("{e}"))?;
        let conn = self.conn.into_inner().unwrap();

        let t = std::time::Instant::now();
        create_views(&conn, has_bam, is_mag_mode, has_annotations)?;
        let create_views_secs = t.elapsed().as_secs_f64();

        let t = std::time::Instant::now();
        cleanup_unused_variables(&conn)?;
        let cleanup_vars_secs = t.elapsed().as_secs_f64();

        let drop_empty_secs = if has_bam {
            let t = std::time::Instant::now();
            drop_empty_tables(&conn)?;
            t.elapsed().as_secs_f64()
        } else {
            0.0
        };

        let (update_counts_secs, zero_coverage_contigs, zero_coverage_mags) = if has_bam {
            let t = std::time::Instant::now();
            update_sample_counts(&conn, is_mag_mode)?;
            let elapsed = t.elapsed().as_secs_f64();
            let zero_contigs: u64 = if !is_mag_mode {
                conn.query_row(
                    "SELECT COUNT(*) FROM Contig WHERE Number_of_samples = 0",
                    [],
                    |r| r.get(0),
                ).unwrap_or(0)
            } else { 0 };
            let zero_mags: u64 = if is_mag_mode {
                conn.query_row(
                    "SELECT COUNT(*) FROM MAG WHERE Number_of_samples = 0",
                    [],
                    |r| r.get(0),
                ).unwrap_or(0)
            } else { 0 };
            (elapsed, zero_contigs, zero_mags)
        } else {
            (0.0, 0, 0)
        };

        let t = std::time::Instant::now();
        conn.execute("FORCE CHECKPOINT", [])
            .context("Failed to checkpoint database")?;
        let checkpoint_secs = t.elapsed().as_secs_f64();
        if let Some(report) = report { report.write(true)?; }

        Ok(FinalizeTimings {
            create_views_secs,
            cleanup_vars_secs,
            drop_empty_secs,
            update_counts_secs,
            checkpoint_secs,
            zero_coverage_contigs,
            zero_coverage_mags,
        })
    }

    /// Get the contig ID for a contig name.
    pub fn get_contig_id(&self, name: &str) -> Option<i64> {
        self.contig_name_to_id.get(name).copied()
    }

    pub fn contig_name_to_id_map(&self) -> &HashMap<String, i64> {
        &self.contig_name_to_id
    }

    /// Get all contig names currently in the database.
    pub fn contig_names(&self) -> Vec<String> {
        self.contig_name_to_id.keys().cloned().collect()
    }

    pub fn create_serve_indexes(&self) -> Result<f64> {
        const SERVE_INDEXES: &[(&str, &str, &str)] = &[
            ("idx_annot_core",           "Contig_annotation_core",  "(Contig_id, \"Type\")"),
            ("idx_mag_contigs_contig",   "MAG_contigs_association", "(Contig_id)"),
            ("idx_annot_qualifier_aid",  "Annotation_qualifier",    "(Annotation_id)"),
            ("idx_annot_qualifier_key",  "Annotation_qualifier",    "(\"Key\")"),
            ("idx_contig_qualifier_key", "Contig_qualifier",        "(\"Key\")"),
            ("idx_phage_termini_pkg",    "Phage_termini",           "(Packaging_id)"),
        ];

        eprintln!();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("### Building query indexes...");
        eprintln!();

        let conn = self.conn.lock().map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let total_start = std::time::Instant::now();

        for &(idx_name, table_name, columns) in SERVE_INDEXES {
            let exists: bool = conn.query_row(
                "SELECT COUNT(*) > 0 FROM information_schema.tables WHERE table_name = ?",
                params![table_name],
                |row| row.get(0),
            ).unwrap_or(false);
            if !exists {
                eprintln!("[index] {} on {} (skipped, table absent)", idx_name, table_name);
                continue;
            }

            let t = std::time::Instant::now();
            conn.execute(
                &format!("CREATE INDEX IF NOT EXISTS {} ON {} {}", idx_name, table_name, columns),
                [],
            ).with_context(|| format!("Failed to create index {} on {}", idx_name, table_name))?;

            eprintln!("[index] {} on {}: {:.1}s", idx_name, table_name, t.elapsed().as_secs_f64());
        }

        eprintln!();
        let total_secs = total_start.elapsed().as_secs_f64();
        Ok(total_secs)
    }
}

/// Populate Number_of_samples: count distinct samples with any coverage entry.
/// In MAG mode: stored on MAG (counts across all member contigs).
/// In contig mode: stored on Contig (counts per contig directly).
fn update_sample_counts(conn: &Connection, is_mag_mode: bool) -> Result<()> {
    if is_mag_mode {
        conn.execute(
            "UPDATE MAG SET Number_of_samples = (
                SELECT COUNT(DISTINCT c.Sample_id)
                FROM Coverage c
                JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id
                WHERE mca.MAG_id = MAG.MAG_id
            )",
            [],
        )
        .context("Failed to update MAG.Number_of_samples")?;
    } else {
        conn.execute(
            "UPDATE Contig SET Number_of_samples = (
                SELECT COUNT(DISTINCT Sample_id) FROM Coverage
                WHERE Coverage.Contig_id = Contig.Contig_id
            )",
            [],
        )
        .context("Failed to update Contig.Number_of_samples")?;
    }
    Ok(())
}

/// Calculate the minimal distance between any kept terminus center position
/// and any terminase gene annotation for a given contig.
/// Returns None if no terminase genes are found or no terminus positions are provided.
/// When circular=true, uses circular distance (min of linear or wrap-around).
fn calculate_terminase_distance_cached(
    cache: &HashMap<i64, Vec<(i32, i32)>>,
    contig_id: i64,
    terminus_positions: &[i32],
    genome_length: i32,
    circular: bool,
) -> Option<i32> {
    if terminus_positions.is_empty() {
        return None;
    }

    let terminase_genes = cache.get(&contig_id)?;
    if terminase_genes.is_empty() {
        return None;
    }

    let mut min_dist = i32::MAX;
    for &term_pos in terminus_positions {
        for &(gene_start, gene_end) in terminase_genes {
            let linear_dist = if term_pos < gene_start {
                gene_start - term_pos
            } else if term_pos > gene_end {
                term_pos - gene_end
            } else {
                0
            };
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

/// Restore mapping tables omitted by genome-only creation or empty-table cleanup.
/// Shared with initial creation so extend mode uses the same schemas and preserves data.
fn create_mapping_tables(conn: &Connection, is_mag_mode: bool) -> Result<()> {
    conn.execute(
        "CREATE TABLE IF NOT EXISTS Sample (
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
        "CREATE TABLE IF NOT EXISTS Coverage (
            Contig_id INTEGER REFERENCES Contig(Contig_id),
            Sample_id INTEGER REFERENCES Sample(Sample_id),
            Aligned_fraction_percentage INTEGER,
            Expected_aligned_fraction INTEGER,
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
        "CREATE TABLE IF NOT EXISTS Phage_mechanisms (
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
        "CREATE TABLE IF NOT EXISTS Phage_termini (
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
            Passed_PoissonTest BOOLEAN NOT NULL,
            Expected_SPC INTEGER NOT NULL,
            Pvalue TEXT NOT NULL,
            Adjusted_pvalue TEXT NOT NULL,
            Passed_ClippingTest BOOLEAN NOT NULL,
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
        "CREATE TABLE IF NOT EXISTS Feature_blob (
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
        "CREATE TABLE IF NOT EXISTS Feature_blob_chunk (
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

    if is_mag_mode {
        // MAG-level coverage metrics recomputed from raw per-position coverage
        // during MAG-by-MAG processing. Written directly by the writer thread.
        conn.execute(
            "CREATE TABLE IF NOT EXISTS MAG_coverage (
                MAG_id                               INTEGER NOT NULL REFERENCES MAG(MAG_id),
                Sample_id                            INTEGER NOT NULL REFERENCES Sample(Sample_id),
                MAG_length                           INTEGER NOT NULL,
                Aligned_fraction_percentage          INTEGER NOT NULL,
                Expected_aligned_fraction             INTEGER NOT NULL,
                Read_count                           INTEGER NOT NULL,
                Coverage_mean                        INTEGER NOT NULL,
                Coverage_median                      INTEGER NOT NULL,
                Coverage_trimmed_mean                INTEGER NOT NULL,
                Coverage_coefficient_of_variation    INTEGER NOT NULL,
                Coverage_relative_coverage_roughness INTEGER NOT NULL,
                PRIMARY KEY (MAG_id, Sample_id)
            )",
            [],
        )
        .context("Failed to create MAG_coverage table")?;
    }

    // Misassembly table - counts at ≥50% prevalence threshold
    conn.execute(
        "CREATE TABLE IF NOT EXISTS Misassembly (
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
        "CREATE TABLE IF NOT EXISTS Microdiversity (
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
        "CREATE TABLE IF NOT EXISTS Side_misassembly (
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
        "CREATE TABLE IF NOT EXISTS Topology (
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
    Ok(())
}

/// Ensure inter-contig hits are available in both views and legacy extensions.
fn create_intercontig_hit_table(conn: &Connection) -> Result<()> {
    let legacy: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM duckdb_constraints() WHERE table_name='Contig_blast_hits' AND constraint_type='CHECK' AND contains(constraint_text, '<='))", [], |r| r.get(0))?;
    anyhow::ensure!(!legacy, "This database stores undirected BLAST hits. Rebuild it to use query-only hit tracks; legacy databases cannot be extended with this hit format.");

        // Directed inter-contig BLAST hits: ID 1 is the query, ID 2 the subject.
        // Do not mirror rows: each alignment contributes only to its query track.
        // Pident uses the Contig_directRepeats Pident scale from COLUMN_SCALES.
        conn.execute(
            "CREATE TABLE IF NOT EXISTS Contig_blast_hits (
                Contig_id_1 INTEGER NOT NULL REFERENCES Contig(Contig_id),
                Contig_id_2 INTEGER NOT NULL REFERENCES Contig(Contig_id),
                Position1 INTEGER NOT NULL,
                Position2 INTEGER NOT NULL,
                Position1prime INTEGER NOT NULL,
                Position2prime INTEGER NOT NULL,
                Pident INTEGER NOT NULL,
                PRIMARY KEY (Contig_id_1, Contig_id_2, Position1, Position2, Position1prime, Position2prime)
            )",
            [],
        )
        .context("Failed to create Contig_blast_hits table")?;
    Ok(())
}

fn create_core_tables(conn: &Connection, has_bam: bool, is_mag_mode: bool, has_annotations: bool) -> Result<()> {
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
    // Number_of_samples is only relevant in contig-mode (in MAG-mode it lives on the MAG table)
    if !is_mag_mode {
        conn.execute("ALTER TABLE Contig ADD COLUMN Number_of_samples INTEGER", [])
            .context("Failed to add Number_of_samples to Contig table")?;
    }

    if is_mag_mode {
        // MAG bin table. Metric columns reuse the same integer encoding as the
        // Contig table so Python decoders can share the scale constants
        // (all metric columns ×10_000; see COLUMN_SCALES in types.rs).
        conn.execute(
            "CREATE TABLE MAG (
                MAG_id INTEGER PRIMARY KEY,
                MAG_name TEXT UNIQUE NOT NULL,
                MAG_length INTEGER NOT NULL,
                Number_of_contigs INTEGER NOT NULL,
                N50 INTEGER NOT NULL,
                Duplication_percentage INTEGER,
                GC_mean INTEGER,
                GC_sd INTEGER,
                GC_skew_amplitude INTEGER,
                Positive_GC_skew_windows_percentage INTEGER,
                Number_of_samples INTEGER
            )",
            [],
        )
        .context("Failed to create MAG table")?;

        // Associations are written longest-contig-first within each MAG. That
        // ordering is the contract the visualization uses for the MAG track
        // (longest → shortest) and for picking the "best hit partner".
        // Offset_in_MAG is the cumulative bp offset of this contig within the
        // MAG under the longest-first ordering — so per-contig coordinates can
        // be mapped to MAG-space with a single read (no window function needed).
        conn.execute(
            "CREATE TABLE MAG_contigs_association (
                MAG_id INTEGER NOT NULL REFERENCES MAG(MAG_id),
                Contig_id INTEGER NOT NULL REFERENCES Contig(Contig_id),
                Offset_in_MAG INTEGER NOT NULL,
                PRIMARY KEY (MAG_id, Contig_id)
            )",
            [],
        )
        .context("Failed to create MAG_contigs_association table")?;

    }

    create_intercontig_hit_table(conn)?;

    if has_bam {
        create_mapping_tables(conn, is_mag_mode)?;
    }


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
    {
        let mut appender = conn.appender("Column_scales")
            .context("Failed to create Column_scales appender")?;
        // SQL INTEGER columns — from COLUMN_SCALES in types.rs
        for (feature, column, scale) in crate::types::COLUMN_SCALES {
            appender.append_row(params![*feature, *column, *scale])
                .context("Failed to append Column_scales row")?;
        }
        // BLOB feature values — from VARIABLES in types.rs
        for v in crate::types::VARIABLES {
            appender.append_row(params![v.name, "Value", v.value_scale.multiplier() as i32])
                .context("Failed to append Column_scales row")?;
        }
        // BLOB sparse metadata scales
        appender.append_row(params!["Sparse_metadata", "Metadata_statistics", crate::types::METADATA_STATS_SCALE as i32])?;
        appender.append_row(params!["Sparse_metadata", "Sequence_prevalence", crate::types::METADATA_PREVALENCE_SCALE as i32])?;
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



    if has_annotations {
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

    conn.execute(
        "CREATE TABLE Annotation_qualifier (
            Annotation_id INTEGER,
            \"Key\" TEXT,
            \"Value\" TEXT
        )",
        [],
    )
    .context("Failed to create Annotation_qualifier table")?;

    conn.execute(
        "CREATE TABLE Contig_qualifier (
            Contig_id INTEGER,
            \"Key\" TEXT,
            \"Value\" TEXT
        )",
        [],
    )
    .context("Failed to create Contig_qualifier table")?;

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

    let has_mag_assoc: bool = conn.query_row(
        "SELECT COUNT(*) > 0 FROM information_schema.tables
         WHERE table_name = 'MAG_contigs_association'",
        [],
        |r| r.get(0),
    ).unwrap_or(false);
    if has_mag_assoc {
        conn.execute(
            "CREATE VIEW MAG_annotation_core AS
             SELECT mca.MAG_id,
                    ca.Annotation_id,
                    ca.Contig_id,
                    ca.\"Start\" + mca.Offset_in_MAG AS \"Start\",
                    ca.\"End\"   + mca.Offset_in_MAG AS \"End\",
                    ca.Strand,
                    ca.\"Type\",
                    ca.Main_isoform,
                    ca.Parent_annotation_id
             FROM Contig_annotation_core ca
             JOIN MAG_contigs_association mca ON mca.Contig_id = ca.Contig_id",
            [],
        )
        .context("Failed to create MAG_annotation_core view")?;
    }
    } // end if has_annotations

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

    if has_annotations {
    conn.execute(
        "CREATE TABLE Annotated_types (
            Type_id INTEGER PRIMARY KEY,
            Type_name TEXT UNIQUE,
            Frequency INTEGER
        )",
        [],
    )
    .context("Failed to create Annotated_types table")?;
    }

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
fn insert_contigs(conn: &Connection, contigs: &[ContigInfo], is_mag_mode: bool) -> Result<()> {
    insert_contigs_with_offset(conn, contigs, is_mag_mode, 0)
}

fn insert_contigs_with_offset(
    conn: &Connection,
    contigs: &[ContigInfo],
    is_mag_mode: bool,
    offset: i64,
) -> Result<()> {
    let mut appender = conn.appender("Contig")
        .context("Failed to create Contig appender")?;

    for (i, contig) in contigs.iter().enumerate() {
        let null_int: Option<i32> = None;
        if is_mag_mode {
            appender.append_row(params![
                offset + (i + 1) as i64,
                &contig.name,
                contig.length as i64,
                null_int,  // Duplication_percentage - set later from autoblast
                null_int,  // GC_mean - set later from GC content computation (×10_000)
                null_int,  // GC_sd - set later from GC content computation (×10_000)
                null_int,  // GC_skew_amplitude - set later from GC skew computation (×10_000)
                null_int,  // Positive_GC_skew_windows_percentage - set later from GC skew computation (×10_000)
            ])
            .with_context(|| format!("Failed to append contig: {}", contig.name))?;
        } else {
            appender.append_row(params![
                offset + (i + 1) as i64,
                &contig.name,
                contig.length as i64,
                null_int,  // Duplication_percentage - set later from autoblast
                null_int,  // GC_mean - set later from GC content computation (×10_000)
                null_int,  // GC_sd - set later from GC content computation (×10_000)
                null_int,  // GC_skew_amplitude - set later from GC skew computation (×10_000)
                null_int,  // Positive_GC_skew_windows_percentage - set later from GC skew computation (×10_000)
                null_int,  // Number_of_samples - set later
            ])
            .with_context(|| format!("Failed to append contig: {}", contig.name))?;
        }
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
            let seq_str = String::from_utf8_lossy(seq);
            appender.append_row(params![(i + 1) as i64, seq_str.as_ref()])?;
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
            let seq_str = String::from_utf8_lossy(seq);
            appender.append_row(params![id_offset + 1 + i as i64, seq_str.as_ref()])?;
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
    let offset: i64 = conn.query_row("SELECT coalesce(max(Annotation_id), 0) FROM Contig_annotation_core", [], |r| r.get(0))?;
    insert_annotation_batch(conn, annotations, offset)
}

fn insert_annotation_batch(conn: &Connection, annotations: &[FeatureAnnotation], offset: i64) -> Result<()> {
    if annotations.is_empty() { return Ok(()); }
    // PHAROKKA function categories (must match plotting_data_per_sample.py)
    let pharokka_keys = vec![
        "vfdb_card",
        "unknown function",
        "other",
        "tail",
        "transcription regulation",
        "DNA, RNA and nucleotide metabolism",
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
    let ids: Vec<i64> = (offset + 1..=offset + annotations.len() as i64).collect();
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
        "INSERT INTO Constants_for_plotting (Constant, Status) VALUES ('isoforms', ?) ON CONFLICT (Constant) DO UPDATE SET Status = Constants_for_plotting.Status OR excluded.Status",
        params![has_isoforms],
    )
    .context("Failed to insert constants")?;

    // Always insert generic color template
    let mut template_id: i64 = 1;
    let has_generic: bool = conn.query_row("SELECT EXISTS (SELECT 1 FROM Color_templates WHERE Template_id = 1)", [], |r| r.get(0))?;
    if !has_generic {
    conn.execute(
        "INSERT OR IGNORE INTO Color_templates (Template_id, Template_name) VALUES (?, 'generic')",
        params![template_id],
    )
    .context("Failed to insert generic template")?;

    }
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
    if !has_generic {
    for (i, (qname, op, qvalue, color)) in generic_rules.iter().enumerate() {
        conn.execute(
            "INSERT OR IGNORE INTO Color_rules (Rule_id, Template_id, Qualifier_name, Operator, Qualifier_value, Color) VALUES (?, ?, ?, ?, ?, ?)",
            params![i as i64 + 1, template_id, qname, op, qvalue, color],
        )
        .context("Failed to insert generic color rule")?;
    }

    }
    // Insert pharokka color template only if pharokka functions detected
    let has_pharokka: bool = conn.query_row("SELECT EXISTS (SELECT 1 FROM Color_templates WHERE Template_id = 2)", [], |r| r.get(0))?;
    if has_pharokka_function && !has_pharokka {
        template_id = 2;
        conn.execute(
            "INSERT OR IGNORE INTO Color_templates (Template_id, Template_name) VALUES (?, 'pharokka')",
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
            ("function", "=", "DNA, RNA and nucleotide metabolism", "#f000ff"),
            ("function", "=", "lysis", "#001eff"),
            ("function", "=", "moron, auxiliary metabolic gene and host takeover", "#8900ff"),
            ("function", "=", "integration and excision", "#E0B0FF"),
            ("function", "=", "head and packaging", "#ff008d"),
            ("function", "=", "connector", "#5A5A5A"),
        ];
        let rule_offset = generic_rules.len() as i64;
        for (i, (qname, op, qvalue, color)) in pharokka_rules.iter().enumerate() {
            conn.execute(
                "INSERT OR IGNORE INTO Color_rules (Rule_id, Template_id, Qualifier_name, Operator, Qualifier_value, Color) VALUES (?, ?, ?, ?, ?, ?)",
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
            "hit_count_within_mag" => "Contig_blob".to_string(),
            "hit_identity_within_mag" => "Contig_blob".to_string(),
            "hit_count" | "hit_identity" => "Contig_blob".to_string(),
            "gc_content" => "Contig_blob".to_string(),
            "gc_skew" => "Contig_blob".to_string(),
            _ => "Feature_blob".to_string(),
        };

        conn.execute(
            "INSERT OR IGNORE INTO Variable (Variable_id, Variable_name, Subplot, Module, Module_order, \"Type\", Encoding, Value_scale, Color, Alpha, Fill_alpha, \"Size\", Title, Help, Feature_table_name)
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
fn create_views(conn: &Connection, has_bam: bool, is_mag_mode: bool, _has_annotations: bool) -> Result<()> {
    // All per-sample features (including primary_reads) are now stored in Feature_blob.
    // Repeat features are converted to Contig_blob via DbWriter::convert_repeat_blobs,
    // called in pre-sample setup so MAG-scale aggregation can read them.
    use crate::types::get_column_scale;
    let s_af  = get_column_scale("Coverage", "Aligned_fraction_percentage");
    let s_eaf = get_column_scale("Coverage", "Expected_aligned_fraction");
    let s_m   = get_column_scale("Coverage", "Coverage_mean");
    let s_med = get_column_scale("Coverage", "Coverage_median");
    let s_tm  = get_column_scale("Coverage", "Coverage_trimmed_mean");
    let s_cv  = get_column_scale("Coverage", "Coverage_coefficient_of_variation");
    let s_rcr = get_column_scale("Coverage", "Coverage_relative_coverage_roughness");
    let s_sc  = get_column_scale("Side_misassembly", "Contig_start_collapse_percentage");
    let s_ec  = get_column_scale("Side_misassembly", "Contig_end_collapse_percentage");

    if has_bam {
    // Explicit_coverage VIEW with RPKM and TPM
    conn.execute(
        &format!("CREATE VIEW Explicit_coverage AS
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
             p.Contig_id,
             p.Sample_id,
             c.Contig_name,
             s.Sample_name,
             p.Aligned_fraction_percentage / {s_af:.1} AS Aligned_fraction_percentage,
             p.Expected_aligned_fraction / {s_eaf:.1} AS Expected_aligned_fraction,
             p.Read_count,
             p.Coverage_mean / {s_m:.1} AS Coverage_mean,
             p.Coverage_median / {s_med:.1} AS Coverage_median,
             p.Coverage_trimmed_mean / {s_tm:.1} AS Coverage_trimmed_mean,
             rb.RPKM,
             CASE WHEN rs.total_rpkm > 0 THEN (rb.RPKM / rs.total_rpkm) * 1e6 ELSE 0 END AS TPM,
             ROUND(p.Coverage_coefficient_of_variation / {s_cv:.1}, 2) AS Coverage_coefficient_of_variation,
             ROUND(p.Coverage_relative_coverage_roughness / {s_rcr:.1}, 4) AS Relative_coverage_roughness
         FROM Coverage p
         JOIN Contig c ON p.Contig_id = c.Contig_id
         JOIN Sample s ON p.Sample_id = s.Sample_id
         JOIN rpkm_base rb ON p.Contig_id = rb.Contig_id AND p.Sample_id = rb.Sample_id
         JOIN rpkm_sum rs ON p.Sample_id = rs.Sample_id"),
        [],
    )
    .context("Failed to create Explicit_coverage VIEW")?;

    // Explicit_misassembly VIEW - per 100kbp normalization
    conn.execute(
        "CREATE VIEW Explicit_misassembly AS
         SELECT
             m.Contig_id,
             m.Sample_id,
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
             md.Contig_id,
             md.Sample_id,
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
        &format!("CREATE VIEW Explicit_side_misassembly AS
         SELECT
             sm.Contig_id,
             sm.Sample_id,
             c.Contig_name,
             s.Sample_name,
             COALESCE(sm.Coverage_first_position, 0) AS Coverage_first_position,
             COALESCE(sm.Contig_start_collapse_percentage, 0) / {s_sc:.1} AS Contig_start_collapse_prevalence,
             COALESCE(sm.Contig_start_collapse_bp, 0) AS Contig_start_collapse_bp,
             COALESCE(sm.Contig_start_expansion_bp, 0) AS Contig_start_expansion_bp,
             COALESCE(sm.Coverage_last_position, 0) AS Coverage_last_position,
             COALESCE(sm.Contig_end_collapse_percentage, 0) / {s_ec:.1} AS Contig_end_collapse_prevalence,
             COALESCE(sm.Contig_end_collapse_bp, 0) AS Contig_end_collapse_bp,
             COALESCE(sm.Contig_end_expansion_bp, 0) AS Contig_end_expansion_bp,
             CASE WHEN s.Sequencing_type = 'paired-short' THEN COALESCE(sm.Contig_end_misjoint_mates, 0) ELSE NULL END AS Contig_end_misjoint_mates,
             CASE WHEN s.Sequencing_type != 'paired-short' THEN NULL WHEN cov.Coverage_mean > 0 THEN COALESCE(sm.Contig_end_misjoint_mates, 0) * 1000.0 / cov.Coverage_mean ELSE 0 END AS Normalized_contig_end_misjoint_mates
         FROM Side_misassembly sm
         JOIN Contig c ON sm.Contig_id = c.Contig_id
         JOIN Sample s ON sm.Sample_id = s.Sample_id
         LEFT JOIN Coverage cov ON sm.Contig_id = cov.Contig_id AND sm.Sample_id = cov.Sample_id"),
        [],
    )
    .context("Failed to create Explicit_side_misassembly VIEW")?;

    // Explicit_topology VIEW (JOINs with Coverage for normalization)
    conn.execute(
        "CREATE VIEW Explicit_topology AS
         SELECT
             t.Contig_id,
             t.Sample_id,
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
             p.Contig_id,
             p.Sample_id,
             c.Contig_name,
             s.Sample_name,
             m.Packaging_mechanism,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM Phage_termini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'start' AND pt.Passed_PoissonTest = true AND pt.Passed_ClippingTest = true),
                 ''
             ) AS Left_termini,
             m.Median_left_termini_clippings,
             COALESCE(
                 (SELECT STRING_AGG(CAST(pt.Center AS VARCHAR), ',' ORDER BY pt.Center)
                  FROM Phage_termini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Status = 'end' AND pt.Passed_PoissonTest = true AND pt.Passed_ClippingTest = true),
                 ''
             ) AS Right_termini,
             m.Median_right_termini_clippings,
             CASE WHEN m.Duplication = true THEN 'DTR'
                  WHEN m.Duplication = false THEN 'ITR'
                  ELSE NULL END AS Duplication,
             COALESCE(
                 (SELECT COUNT(*)
                  FROM Phage_termini pt
                  WHERE pt.Packaging_id = m.Packaging_id AND pt.Passed_PoissonTest = true AND pt.Passed_ClippingTest = true),
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

    if has_bam && is_mag_mode {
        // Per-MAG aggregates. Counts: SUM. Per-100kbp: recomputed from
        // SUM(count) × 100000 / SUM(Contig_length). Coverage metrics come
        // from the MAG_coverage table, populated during MAG-by-MAG processing.

        conn.execute(
            &format!("CREATE VIEW Explicit_coverage_per_MAG AS
             WITH rpkm_base_mag AS (
                 SELECT
                     mca.MAG_id,
                     p.Sample_id,
                     SUM(CASE WHEN c.Contig_length > 0 AND s.Number_of_mapped_reads > 0
                              THEN (CAST(p.Read_count AS DOUBLE) * 1e9) / (CAST(c.Contig_length AS DOUBLE) * CAST(s.Number_of_mapped_reads AS DOUBLE))
                              ELSE 0 END) AS RPKM
                 FROM Coverage p
                 JOIN Contig c ON p.Contig_id = c.Contig_id
                 JOIN Sample s ON p.Sample_id = s.Sample_id
                 JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id
                 GROUP BY mca.MAG_id, p.Sample_id
             ),
             rpkm_sum_mag AS (
                 SELECT Sample_id, SUM(RPKM) AS total_rpkm FROM rpkm_base_mag GROUP BY Sample_id
             )
             SELECT
                 mc.MAG_id,
                 mc.Sample_id,
                 mg.MAG_name,
                 s.Sample_name,
                 mc.Aligned_fraction_percentage / {s_af:.1}                                                AS Aligned_fraction_percentage,
                 mc.Expected_aligned_fraction / {s_eaf:.1}                                                AS Expected_aligned_fraction,
                 mc.Read_count                                                                           AS Read_count,
                 mc.Coverage_mean / {s_m:.1}                                                              AS Coverage_mean,
                 mc.Coverage_median / {s_med:.1}                                                          AS Coverage_median,
                 mc.Coverage_trimmed_mean / {s_tm:.1}                                                     AS Coverage_trimmed_mean,
                 rbm.RPKM                                                                                AS RPKM,
                 CASE WHEN rsm.total_rpkm > 0 THEN (rbm.RPKM / rsm.total_rpkm) * 1e6 ELSE 0 END          AS TPM,
                 ROUND(mc.Coverage_coefficient_of_variation / {s_cv:.1}, 2)                              AS Coverage_coefficient_of_variation,
                 ROUND(mc.Coverage_relative_coverage_roughness / {s_rcr:.1}, 4)                           AS Relative_coverage_roughness
             FROM MAG_coverage mc
             JOIN MAG mg    ON mg.MAG_id   = mc.MAG_id
             JOIN Sample s  ON s.Sample_id = mc.Sample_id
             JOIN rpkm_base_mag rbm ON rbm.MAG_id = mc.MAG_id AND rbm.Sample_id = mc.Sample_id
             JOIN rpkm_sum_mag  rsm ON rsm.Sample_id = mc.Sample_id"),
            [],
        )
        .context("Failed to create Explicit_coverage_per_MAG VIEW")?;

        conn.execute(
            "CREATE VIEW Explicit_misassembly_per_MAG AS
             SELECT
                 mg.MAG_id,
                 s.Sample_id,
                 mg.MAG_name,
                 s.Sample_name,
                 SUM(m.Mismatches_count) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Mismatches_per_100kbp,
                 SUM(m.Deletions_count)  * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Deletions_per_100kbp,
                 SUM(m.Insertions_count) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Insertions_per_100kbp,
                 SUM(m.Clippings_count)  * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Clippings_per_100kbp,
                 SUM(m.Collapse_bp)  AS Collapse_bp,
                 SUM(m.Collapse_bp)  * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Collapse_per_100kbp,
                 SUM(m.Expansion_bp) AS Expansion_bp,
                 SUM(m.Expansion_bp) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Expansion_per_100kbp
             FROM Misassembly m
             JOIN Contig c ON m.Contig_id = c.Contig_id
             JOIN Sample s ON m.Sample_id = s.Sample_id
             JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id
             JOIN MAG mg ON mg.MAG_id = mca.MAG_id
             GROUP BY mg.MAG_id, s.Sample_id, mg.MAG_name, s.Sample_name",
            [],
        )
        .context("Failed to create Explicit_misassembly_per_MAG VIEW")?;

        conn.execute(
            "CREATE VIEW Explicit_microdiversity_per_MAG AS
             SELECT
                 mg.MAG_id,
                 s.Sample_id,
                 mg.MAG_name,
                 s.Sample_name,
                 SUM(md.Mismatches_count) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Mismatches_per_100kbp,
                 SUM(md.Deletions_count)  * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Deletions_per_100kbp,
                 SUM(md.Insertions_count) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Insertions_per_100kbp,
                 SUM(md.Clippings_count)  * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Clippings_per_100kbp,
                 SUM(md.Microdiverse_bp_on_reads)     AS Microdiverse_bp_on_reads,
                 SUM(md.Microdiverse_bp_on_reads)     * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Microdiverse_bp_per_100kbp_on_reads,
                 SUM(md.Microdiverse_bp_on_reference) AS Microdiverse_bp_on_reference,
                 SUM(md.Microdiverse_bp_on_reference) * 100000.0 / NULLIF(SUM(c.Contig_length), 0) AS Microdiverse_bp_per_100kbp_on_reference
             FROM Microdiversity md
             JOIN Contig c ON md.Contig_id = c.Contig_id
             JOIN Sample s ON md.Sample_id = s.Sample_id
             JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id
             JOIN MAG mg ON mg.MAG_id = mca.MAG_id
             GROUP BY mg.MAG_id, s.Sample_id, mg.MAG_name, s.Sample_name",
            [],
        )
        .context("Failed to create Explicit_microdiversity_per_MAG VIEW")?;

    } // end if has_bam && is_mag_mode

    {
        // Query-side access; reciprocal BLAST rows remain separate directed hits.
        conn.execute(
            "CREATE VIEW Contig_blast_hits_query AS
                SELECT Contig_id_1 AS Contig_id, Position1, Position2,
                       Contig_id_2 AS Contig_id_other, Position1prime, Position2prime, Pident
                FROM Contig_blast_hits",
            [],
        ).context("Failed to create Contig_blast_hits_query view")?;
    }

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

    for (var_name, table_name) in rows {
        if table_name == "Feature_blob" {
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
            continue;
        }

        conn.execute("DELETE FROM Variable WHERE Variable_name = ?1", params![var_name])?;
    }

    Ok(())
}

/// Drop empty module-dependent tables and their views.
/// Prevents empty sections in the Filtering/Summary UI when modules weren't calculated.
fn drop_empty_tables(conn: &Connection) -> Result<()> {
    let table_view_map: &[(&str, &[&str])] = &[
        ("Misassembly", &["Explicit_misassembly", "Explicit_misassembly_per_MAG"]),
        ("Microdiversity", &["Explicit_microdiversity", "Explicit_microdiversity_per_MAG"]),
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

/// Repeats data from BLAST self-alignment.
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
    pub contig_id: i64,
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
    pub contig_id: i64,
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
    pub contig_id: i64,
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
    pub contig_id: i64,
    pub contig_name: String,
    pub circularising_reads: Option<u64>,
    pub circularising_reads_percentage: Option<i32>,
    pub median_circularising_len: Option<i64>,
    pub circularising_inserts: Option<u64>,
    pub circularising_insert_size_deviation: Option<i32>,
}

/// Parse BLAST outfmt 6 file and return repeats data.
/// Used by the standalone `thebigbam autoblast` CLI command (Python/BLAST+).
/// The main `calculate` pipeline uses minimap2 instead (see processing.rs).
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

#[cfg(test)]
mod annotation_import_tests {
    use super::*;

    #[test]
    fn memory_budget_reserves_rust_and_worker_headroom() {
        assert_eq!(memory_budget_mib(273 * 1024, 172_033, None).unwrap(), 16 * 1024);
        assert_eq!(memory_budget_mib(1024, 0, None).unwrap(), 256);
        assert_eq!(memory_budget_mib(10_000, 7_000, Some(8_000)).unwrap(), 1_000);
        assert!(memory_budget_mib(1024, 900, None).is_err());
        assert!(memory_budget_mib(1024, 0, Some(0)).is_err());
    }

    fn annotation(contig: i64, kind: &str, locus: Option<&str>, own: Option<&str>, parent: Option<&str>) -> FeatureAnnotation {
        FeatureAnnotation {
            contig_id: contig, start: 1, end: 9, strand: -1, feature_type: kind.into(),
            qualifiers: locus.map(|s| HashMap::from([("locus_tag".into(), s.into())])).unwrap_or_default(),
            nucleotide_sequence: (kind == "CDS").then(|| "ATGAAATAA".into()),
            protein_sequence: (kind == "CDS").then(|| "MK*".into()),
            segments: vec![(1, 3), (7, 9)],
            self_key: own.map(str::to_owned), parent_key: parent.map(str::to_owned),
        }
    }

    fn rows(writer: &DbWriter, table: &str) -> Vec<String> {
        let conn = writer.lock_conn().unwrap();
        let mut stmt = conn.prepare(&format!("SELECT CAST(to_json(t) AS VARCHAR) AS row FROM {table} t ORDER BY row")).unwrap();
        stmt.query_map([], |r| r.get(0)).unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap()
    }

    #[test]
    fn mag_blast_batches_preserve_query_direction_and_repeat_orientation() {
        let directory = tempfile::tempdir().unwrap();
        let contigs: Vec<_> = ["a", "b"].into_iter().map(|name| ContigInfo {
            name: name.into(), aliases: vec![], length: 100, sequence: None,
        }).collect();
        let writer = DbWriter::create(&directory.path().join("blast.db"), &contigs, &[], &[], false, true).unwrap();
        let hits = vec![
            ("a".into(), "b".into(), 1, 20, 80, 61, 99.0),
            ("b".into(), "a".into(), 80, 61, 1, 20, 88.0),
        ];
        let repeats = vec![
            RepeatsData { contig_name: "a".into(), position1: 1, position2: 20,
                position1prime: 50, position2prime: 69, pident: 97.0, is_direct: true },
            RepeatsData { contig_name: "b".into(), position1: 1, position2: 20,
                position1prime: 80, position2prime: 61, pident: 98.0, is_direct: false },
        ];
        writer.write_mag_blast_batch(&repeats, &hits).unwrap();
        // An exact repeat of the directed row is ignored; reciprocal rows are independent.
        writer.write_mag_blast_batch(&[], &hits[1..]).unwrap();
        writer.write_mag_blast_batch(&[], &[("b".into(), "a".into(), 30, 40, 10, 20, 96.0)]).unwrap();
        let conn = writer.lock_conn().unwrap();
        let mut stmt = conn.prepare("SELECT Contig_id_1, Contig_id_2, Position1, Position2, Position1prime, Position2prime, Pident FROM Contig_blast_hits ORDER BY Position1").unwrap();
        let stored = stmt.query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?, r.get::<_, i32>(2)?, r.get::<_, i32>(3)?, r.get::<_, i32>(4)?, r.get::<_, i32>(5)?, r.get::<_, i32>(6)?)))
            .unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
        let scale = crate::types::get_column_scale("Contig_directRepeats", "Pident");
        assert_eq!(stored, vec![(1, 2, 1, 20, 80, 61, (99.0 * scale).round() as i32),
                              (2, 1, 30, 40, 10, 20, (96.0 * scale).round() as i32),
                              (2, 1, 80, 61, 1, 20, (88.0 * scale).round() as i32)]);
        for table in ["Contig_directRepeats", "Contig_invertedRepeats"] {
            assert_eq!(conn.query_row(&format!("SELECT count(*) FROM {table}"), [], |r| r.get::<_, i64>(0)).unwrap(), 1);
        }
        assert_eq!(conn.query_row("SELECT count(*) FROM Mag_blast_batch", [], |r| r.get::<_, i64>(0)).unwrap(), 0);
        drop(stmt);
        drop(conn);
        writer.write_intercontig_hit_features().unwrap();
        let conn = writer.lock_conn().unwrap();
        for id in [1,2] {
            let mut stmt = conn.prepare("SELECT CAST(Chunk_idx AS INTEGER),Data FROM Contig_blob_chunk WHERE Contig_id=? AND Feature_id=? ORDER BY Chunk_idx").unwrap();
            let chunks = stmt.query_map(params![id,feature_name_to_id("hit_count").unwrap()], |r| Ok((r.get::<_,i32>(0)?,r.get::<_,Vec<u8>>(1)?))).unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
            let (_, values, _, _) = crate::blob::decode_sparse_from_chunks(&chunks).unwrap();
            assert!(!values.is_empty());
            assert!(values.iter().all(|&count| count==1), "Reciprocal row was counted twice for query {id}");
        }

    }

    #[test]
    fn batched_annotations_match_monolithic_including_cross_batch_relationships() {
        let directory = tempfile::tempdir().unwrap();
        let contigs = vec![
            ContigInfo { name: "a".into(), aliases: vec![], length: 9, sequence: Some(b"ATGAAATAA".to_vec()) },
            ContigInfo { name: "b".into(), aliases: vec![], length: 9, sequence: Some(b"TTATTTCAT".to_vec()) },
        ];
        let mut annotations = vec![
            annotation(1, "CDS", Some("late-parent"), None, Some("later")),
            annotation(1, "mRNA", Some("shared"), Some("first"), None),
            annotation(1, "CDS", Some("flat"), None, None),
            annotation(2, "mRNA", Some("late-parent"), Some("later"), None),
            annotation(2, "mRNA", Some("shared"), Some("second"), None),
            annotation(2, "CDS", Some("shared"), None, Some("first")),
            annotation(2, "CDS", Some("shared"), None, Some("second")),
            annotation(2, "CDS", Some("flat"), None, Some("missing")),
            annotation(2, "CDS", None, None, None),
        ];
        annotations[8].qualifiers.insert("function".into(), "tail".into());
        annotations[0].qualifiers.insert("translation".into(), "MK*".into());
        let quals = vec![(1, HashMap::from([("organism".into(), "first organism".into())])),
                         (2, HashMap::from([("organism".into(), "second organism".into())]))];
        let reference = DbWriter::create(&directory.path().join("reference.db"), &contigs, &annotations, &quals, false, true).unwrap();
        let mut streamed = DbWriter::create_streaming(&directory.path().join("streamed.db"), false, true).unwrap();
        streamed.append_annotation_batch(&contigs[..1], &annotations[..3], &quals[..1], 0, 0).unwrap();
        streamed.append_annotation_batch(&contigs[1..], &annotations[3..], &quals[1..], 1, 3).unwrap();
        streamed.finish_annotation_import().unwrap();
        for table in ["Contig", "Contig_sequence", "Contig_annotation_core", "Annotation_sequence", "Annotation_qualifier", "Annotation_segments", "Contig_qualifier", "Annotated_types", "Constants_for_plotting", "Color_templates", "Color_rules"] {
            assert_eq!(rows(&reference, table), rows(&streamed, table), "table {table}");
        }
        assert_eq!(streamed.get_contig_id("b"), Some(2));
    }
}

#[cfg(test)]
mod sample_write_timing_tests {
    use super::*;

    #[test]
    fn timed_blob_write_counts_only_stored_rows_and_payload() {
        let directory = tempfile::tempdir().unwrap();
        let writer = DbWriter::create(&directory.path().join("timed.db"),
            &[ContigInfo { name: "a".into(), aliases: vec![], length: 100, sequence: None }],
            &[], &[], true, false).unwrap();
        let sample = writer.insert_sample("sample", "single-short", 1, 1, false).unwrap();
        let encoded = crate::blob::encode_dense_blob(&vec![1; 100], crate::types::get_value_scale("primary_reads"), 100);
        let expected_zoom = encoded.zoom.len() as u64;
        let expected_chunks = encoded.chunks.len() as u64;
        let expected_bytes = encoded.chunks.iter().filter(|c| encoded.should_store_chunk(c)).map(|c| c.len() as u64).sum::<u64>();
        let empty = crate::blob::encode_dense_blob(&vec![0; 100], crate::types::get_value_scale("primary_reads"), 100);
        let blobs = vec![("primary_reads".into(), 1, encoded), ("primary_reads".into(), 1, empty)];
        let timings = writer.write_sample_data_timed("sample", &[], &[], &[], &[], &[], &[], &blobs, false, true).unwrap();
        assert_eq!((timings.zoom_rows, timings.chunk_rows, timings.zoom_bytes, timings.chunk_bytes),
            (1, expected_chunks, expected_zoom, expected_bytes));
        assert!(timings.blob_append_ns > 0);
        assert!(timings.checkpoint_ns > 0);
        let conn = writer.lock_conn().unwrap();
        let zoom: (u64, u64) = conn.query_row("SELECT count(*), sum(octet_length(Zoom_data)) FROM Feature_blob WHERE Sample_id = ?",
            params![sample], |row| Ok((row.get(0)?, row.get(1)?))).unwrap();
        let chunks: (u64, u64) = conn.query_row("SELECT count(*), sum(octet_length(Data)) FROM Feature_blob_chunk WHERE Sample_id = ?",
            params![sample], |row| Ok((row.get(0)?, row.get(1)?))).unwrap();
        assert_eq!(zoom, (timings.zoom_rows, timings.zoom_bytes));
        assert_eq!(chunks, (timings.chunk_rows, timings.chunk_bytes));
        drop(conn);
        writer.insert_sample("untimed", "single-short", 1, 1, false).unwrap();
        let untimed = writer.write_sample_data_timed("untimed", &[], &[], &[], &[], &[], &[], &blobs, false, false).unwrap();
        assert_eq!((untimed.blob_append_ns, untimed.zoom_rows, untimed.chunk_bytes), (0, 0, 0));
    }
}

#[cfg(test)]
mod size_report_tests {
    use super::*;
    #[test]
    fn sparse_write_preserves_indices_and_reports_only_stored_payload() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("sizes.db");
        let length = 3 * crate::blob::CHUNK_SIZE;
        let writer = DbWriter::create(&path, &[ContigInfo { name: "a".into(), aliases: vec![], length: length as usize, sequence: None }], &[], &[], true, false).unwrap();
        assert!(!directory.path().join("sizes.db.size.tsv").exists());
        assert!(writer.enable_size_report(&path).unwrap());
        assert!(!writer.enable_size_report(&path).unwrap());
        let sample = writer.insert_sample("a", "single-short", 1, 1, false).unwrap();
        let positions = vec![crate::blob::CHUNK_SIZE + 7];
        let values = vec![3];
        let encoded = crate::blob::encode_sparse_blob(&positions, &values, None, crate::blob::MetadataFlags::default(), crate::types::get_value_scale("mismatches"), length);
        let expected_zoom = encoded.zoom.len();
        let expected_chunk = encoded.chunks[1].len();
        let conn = writer.lock_conn().unwrap();
        let mut timing = SampleWriteTimings::default();
        writer.write_feature_blobs_timed(&conn, sample, &[("mismatches".into(), 1, encoded)], true, &mut timing).unwrap();
        let mut stmt = conn.prepare("SELECT Chunk_idx, Data FROM Feature_blob_chunk ORDER BY Chunk_idx").unwrap();
        let chunks = stmt.query_map([], |r| Ok((r.get::<_, i32>(0)?, r.get::<_, Vec<u8>>(1)?))).unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].0, 1);
        let (decoded_pos, decoded_values, _, _) = crate::blob::decode_sparse_from_chunks(&chunks).unwrap();
        assert_eq!(decoded_pos, positions); assert_eq!(decoded_values, values);
        assert_eq!(timing.chunk_rows, 1);
        let report = std::fs::read_to_string(directory.path().join("sizes.db.size.tsv")).unwrap();
        assert!(report.contains(&format!("Feature_blob\tmismatches\t1\t1\t{}\t{}\t", expected_zoom, expected_chunk)));
        assert!(report.lines().any(|l| l.starts_with("Feature_blob\tmismatches\t") && l.split('\t').nth(8) == Some("2")));
        drop(stmt); drop(conn);
        writer.finalize().unwrap();
        assert!(std::fs::read_to_string(directory.path().join("sizes.db.size.tsv")).unwrap().contains("# status=complete"));
    }
}

#[cfg(test)]
mod repeat_conversion_memory_tests {
    use super::*;

    fn writer(path: &Path) -> DbWriter {
        let contigs = ["a", "b"].iter().map(|name| ContigInfo {
            name: name.to_string(), aliases: vec![], length: 100, sequence: None,
        }).collect::<Vec<_>>();
        DbWriter::create(path, &contigs, &[], &[], false, false).unwrap()
    }
    fn assert_no_conversion_tables(writer: &DbWriter) {
        let conn = writer.lock_conn().unwrap();
        let count: u64 = conn.query_row(
            "SELECT count(*) FROM information_schema.tables WHERE table_name IN
             ('Repeat_conversion_ids', 'Contig_direct_repeat_count', 'Contig_inverted_repeat_count',
              'Contig_direct_repeat_identity', 'Contig_inverted_repeat_identity')",
            [], |r| r.get(0)).unwrap();
        assert_eq!(count, 0);
    }
    #[test]
    fn empty_repeats_skip_all_conversion_state() {
        let directory = tempfile::tempdir().unwrap();
        let writer = writer(&directory.path().join("empty.db"));
        writer.convert_repeat_blobs(&[]).unwrap();
        writer.convert_repeat_blobs(&[1, 2]).unwrap();
        assert_no_conversion_tables(&writer);
        let conn = writer.lock_conn().unwrap();
        assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob", [], |r| r.get::<_, u64>(0)).unwrap(), 0);
    }
    #[test]
    fn conversion_selects_only_requested_contigs_without_duplicate_events() {
        let directory = tempfile::tempdir().unwrap();
        // Exercise either repeat table as the only populated source.
        for direct in [true, false] {
            let writer = writer(&directory.path().join(format!("{direct}.db")));
            let repeats = ["a", "b"].iter().map(|name| RepeatsData {
                contig_name: name.to_string(), position1: 1, position2: 20,
                position1prime: 50, position2prime: 69, pident: 97.0, is_direct: direct,
            }).collect::<Vec<_>>();
            writer.write_repeats(&repeats).unwrap();
            writer.convert_repeat_blobs(&[1, 1]).unwrap();
            assert_no_conversion_tables(&writer);
            {
                let conn = writer.lock_conn().unwrap();
                assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Contig_id = 2", [], |r| r.get::<_, u64>(0)).unwrap(), 0);
                let feature = if direct { "direct_repeat_count" } else { "inverted_repeat_count" };
                let fid = feature_name_to_id(feature).unwrap();
                let mut stmt = conn.prepare("SELECT Chunk_idx, Data FROM Contig_blob_chunk WHERE Contig_id = 1 AND Feature_id = ? ORDER BY Chunk_idx").unwrap();
                let chunks = stmt.query_map(params![fid], |r| Ok((r.get::<_, i32>(0)?, r.get::<_, Vec<u8>>(1)?)))
                    .unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
                let (positions, values, _, _) = crate::blob::decode_sparse_from_chunks(&chunks).unwrap();
                assert_eq!(positions, vec![0, 19, 20]);
                assert_eq!(values, vec![1, 1, 0]);
                assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Contig_id = 1", [], |r| r.get::<_, u64>(0)).unwrap(), 2);
            }
            // The temporary selection table must be gone so a later conversion can run.
            writer.convert_repeat_blobs(&[2]).unwrap();
            assert_no_conversion_tables(&writer);
            let conn = writer.lock_conn().unwrap();
            assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Contig_id = 2", [], |r| r.get::<_, u64>(0)).unwrap(), 2);
        }
    }
    #[test]
    fn repeat_regeneration_replaces_selected_rows_and_clears_removed_repeats() {
        let directory=tempfile::tempdir().unwrap();
        let writer=writer(&directory.path().join("regenerate.db"));
        writer.write_repeats(&[RepeatsData {contig_name:"a".into(),position1:1,position2:20,position1prime:50,position2prime:69,pident:97.,is_direct:true}]).unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        assert_no_conversion_tables(&writer);
        assert_eq!(writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Contig_blob",[],|r|r.get::<_,u64>(0)).unwrap(),2);
        writer.lock_conn().unwrap().execute("DELETE FROM Contig_directRepeats",[]).unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        assert_eq!(writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Contig_blob",[],|r|r.get::<_,u64>(0)).unwrap(),0);
    }

    #[test]
    fn partial_repeat_upgrade_rebuilds_every_existing_blob() {
        let directory=tempfile::tempdir().unwrap();
        let writer=writer(&directory.path().join("upgrade.db"));
        let repeats=["a","b"].iter().map(|name| RepeatsData {contig_name:name.to_string(),position1:1,position2:20,position1prime:50,position2prime:69,pident:97.,is_direct:true}).collect::<Vec<_>>();
        writer.write_repeats(&repeats).unwrap();
        writer.convert_repeat_blobs(&[1,2]).unwrap();
        writer.lock_conn().unwrap().execute_batch("UPDATE Database_metadata SET Value='interval-v1' WHERE Key='Repeat_blob_encoding'; UPDATE Contig_blob SET Zoom_data='old'::BLOB").unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        assert_no_conversion_tables(&writer);
        let conn=writer.lock_conn().unwrap();
        assert_eq!(conn.query_row("SELECT Value FROM Database_metadata WHERE Key='Repeat_blob_encoding'",[],|r|r.get::<_,String>(0)).unwrap(),"interval-v2");
        assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Zoom_data='old'::BLOB",[],|r|r.get::<_,u64>(0)).unwrap(),0);
        assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Contig_id=2",[],|r|r.get::<_,u64>(0)).unwrap(),2);
    }

    #[test]
    fn failed_repeat_regeneration_rolls_back_blob_deletes_and_marker() {
        let directory=tempfile::tempdir().unwrap();
        let writer=writer(&directory.path().join("rollback.db"));
        writer.write_repeats(&[RepeatsData {contig_name:"a".into(),position1:1,position2:20,position1prime:50,position2prime:69,pident:97.,is_direct:true}]).unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        {
            let conn=writer.lock_conn().unwrap();
            conn.execute_batch("UPDATE Contig_blob SET Zoom_data='old'::BLOB; UPDATE Database_metadata SET Value='interval-v1' WHERE Key='Repeat_blob_encoding'; DROP TABLE Contig_blob_chunk; CREATE TABLE Contig_blob_chunk (Contig_id BIGINT,Feature_id SMALLINT,Chunk_idx SMALLINT,Data INTEGER)").unwrap();
        }
        assert!(writer.convert_repeat_blobs(&[1]).is_err());
        assert_no_conversion_tables(&writer);
        let conn=writer.lock_conn().unwrap();
        assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Zoom_data='old'::BLOB",[],|r|r.get::<_,u64>(0)).unwrap(),2);
        assert_eq!(conn.query_row("SELECT Value FROM Database_metadata WHERE Key='Repeat_blob_encoding'",[],|r|r.get::<_,String>(0)).unwrap(),"interval-v1");
    }

    #[test]
    fn out_of_bounds_repeats_do_not_create_a_last_base_peak() {
        let directory=tempfile::tempdir().unwrap();
        let writer=writer(&directory.path().join("bounds.db"));
        writer.write_repeats(&[RepeatsData {contig_name:"a".into(),position1:200,position2:220,position1prime:300,position2prime:320,pident:97.,is_direct:true}]).unwrap();
        writer.convert_repeat_blobs(&[1]).unwrap();
        assert_eq!(writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Contig_blob",[],|r|r.get::<_,u64>(0)).unwrap(),0);
    }

}

#[cfg(test)]
mod extend_mapping_table_tests {
    use super::*;

    #[test]
    fn extend_restores_mapping_tables_and_preserves_existing_samples() {
        for is_mag_mode in [false, true] {
            for has_existing_mapping in [false, true] {
                let directory = tempfile::tempdir().unwrap();
                let path = directory.path().join("extend.db");
                let contigs = [ContigInfo {
                    name: "a".into(), aliases: vec![], length: 100, sequence: None,
                }];
                let writer = DbWriter::create(&path, &contigs, &[], &[], has_existing_mapping, is_mag_mode).unwrap();
                writer.lock_conn().unwrap().execute(
                    "INSERT INTO Database_metadata (Key, Value) VALUES ('View_mode', ?)",
                    params![if is_mag_mode { "mag" } else { "contig" }],
                ).unwrap();
                if has_existing_mapping {
                    writer.insert_sample("existing", "single-short", 1, 1, false).unwrap();
                    // Finalization removes unused module tables; purge preserves their absence.
                    drop_empty_tables(&writer.lock_conn().unwrap()).unwrap();
                }
                drop(writer);

                let writer = DbWriter::open(&path, &[], &[], &[], true).unwrap();
                let sample_id = writer.insert_sample("new", "single-short", 1, 1, false).unwrap();
                assert_eq!(sample_id, if has_existing_mapping { 2 } else { 1 });
                // Exercise all appenders, including modules without results.
                writer.write_sample_data("new", &[], &[], &[], &[], &[], &[], &[], false).unwrap();
                let conn = writer.lock_conn().unwrap();
                for table in ["Phage_mechanisms", "Phage_termini", "Misassembly", "Microdiversity", "Side_misassembly", "Topology"] {
                    assert_eq!(conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [],
                        |r| r.get::<_, i64>(0)).unwrap(), 0);
                }
                if is_mag_mode {
                    assert_eq!(conn.query_row("SELECT COUNT(*) FROM MAG_coverage", [],
                        |r| r.get::<_, i64>(0)).unwrap(), 0);
                }
                assert_eq!(conn.query_row("SELECT COUNT(*) FROM Contig", [],
                    |r| r.get::<_, i64>(0)).unwrap(), 1);
            }
        }
    }
}

#[cfg(test)]
mod bulk_mag_tests {
    use super::*;
    #[test]
    fn bulk_mag_members_preserve_ties_offsets_weighted_stats_and_extension() {
        let directory = tempfile::tempdir().unwrap();
        let contigs = [("a",100), ("b",100), ("c",50)].iter().map(|(name,length)| ContigInfo {
            name: name.to_string(), aliases: vec![], length: *length, sequence: None }).collect::<Vec<_>>();
        let writer = DbWriter::create(&directory.path().join("mags.db"), &contigs, &[], &[], false, true).unwrap();
        writer.lock_conn().unwrap().execute_batch("UPDATE Contig SET GC_mean = CASE Contig_id WHEN 1 THEN 20 WHEN 2 THEN 40 ELSE NULL END").unwrap();
        let mag = ("first".to_owned(), vec!["b".into(), "c".into(), "a".into()]);
        writer.write_mags(&[mag.clone()]).unwrap();
        writer.write_mags(&[mag, ("second".into(), vec!["c".into()])]).unwrap();
        let conn = writer.lock_conn().unwrap();
        let summary = conn.query_row("SELECT MAG_length, Number_of_contigs, N50, GC_mean FROM MAG WHERE MAG_id=1", [],
            |r| Ok((r.get::<_,i64>(0)?, r.get::<_,i64>(1)?, r.get::<_,i64>(2)?, r.get::<_,i32>(3)?))).unwrap();
        assert_eq!(summary, (250,3,100,30));
        let mut stmt = conn.prepare("SELECT Contig_id, Offset_in_MAG FROM MAG_contigs_association WHERE MAG_id=1 ORDER BY Offset_in_MAG").unwrap();
        let rows = stmt.query_map([], |r| Ok((r.get::<_,i64>(0)?,r.get::<_,i64>(1)?)))
            .unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
        assert_eq!(rows, vec![(2,0),(1,100),(3,200)]);
        assert_eq!(conn.query_row("SELECT count(*) FROM MAG", [], |r| r.get::<_,u64>(0)).unwrap(), 2);
        assert_eq!(conn.query_row("SELECT count(*) FROM MAG WHERE Number_of_samples IS NULL", [],
            |r| r.get::<_,u64>(0)).unwrap(), 2);
        assert_eq!(conn.query_row("SELECT count(*) FROM information_schema.tables WHERE table_name='Mag_member_input'", [], |r| r.get::<_,u64>(0)).unwrap(), 0);
    }
}

#[cfg(test)]
mod staged_gc_tests {
    use super::*;
    #[test]
    fn staged_batches_match_immediate_updates_and_leave_other_contigs_unchanged() {
        let directory = tempfile::tempdir().unwrap();
        let contigs = ["a", "b", "untouched"].iter().map(|name| ContigInfo {
            name: name.to_string(), aliases: vec![], length: 100, sequence: None }).collect::<Vec<_>>();
        let baseline = DbWriter::create(&directory.path().join("baseline.db"), &contigs, &[], &[], false, false).unwrap();
        let staged = DbWriter::create(&directory.path().join("staged.db"), &contigs, &[], &[], false, false).unwrap();
        for (name, mean) in [("a", 42.25), ("b", 60.5)] {
            let data = GCContentData { contig_name: name.into(), gc_values: vec![], skew_values: vec![],
                stats: GCStats { average: mean, sd: 3.25 },
                skew_stats: GCSkewStats { amplitude: 0.25, percent_positive: 51.5 } };
            baseline.update_contig_gc_stats(&[data.clone()]).unwrap();
            staged.stage_contig_gc_stats(&[data]).unwrap();
        }
        staged.flush_contig_gc_stats().unwrap();
        let rows = |writer: &DbWriter| {
            let conn = writer.lock_conn().unwrap();
            let mut stmt = conn.prepare("SELECT GC_mean, GC_sd, GC_skew_amplitude, Positive_GC_skew_windows_percentage FROM Contig ORDER BY Contig_id").unwrap();
            stmt.query_map([], |r| Ok((r.get::<_,Option<i32>>(0)?, r.get::<_,Option<i32>>(1)?,
                r.get::<_,Option<i32>>(2)?, r.get::<_,Option<i32>>(3)?))).unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap()
        };
        assert_eq!(rows(&staged), rows(&baseline));
        assert_eq!(rows(&staged)[2], (None,None,None,None));
    }
}

#[cfg(test)]
mod sample_recovery_tests {
    use super::*;
    #[test]
    fn partial_samples_reuse_ids_clear_payload_and_reject_completed_retries() {
        let directory=tempfile::tempdir().unwrap();
        let path=directory.path().join("recovery.db");
        let contigs=vec![ContigInfo{name:"a".into(),aliases:vec![],length:100,sequence:None}];
        let writer=DbWriter::create(&path,&contigs,&[],&[],true,false).unwrap();
        let id=writer.insert_sample("sample","single-short",10,10,false).unwrap();
        {
            let conn=writer.lock_conn().unwrap();
            writer.write_feature_blobs(&conn,id,&[("primary_reads".into(),1,
                crate::blob::encode_dense_blob(&vec![1;100],crate::types::ValueScale::Raw,100))]).unwrap();
        }
        writer.sample_failed("sample","injected failure").unwrap();
        drop(writer);
        let writer=DbWriter::open(&path,&[],&[],&[],true).unwrap();
        assert_eq!(writer.insert_sample("sample","single-short",20,20,false).unwrap(),id);
        {
            let conn=writer.lock_conn().unwrap();
            for table in ["Feature_blob","Feature_blob_chunk"] {
                assert_eq!(conn.query_row(&format!("SELECT count(*) FROM {table} WHERE Sample_id=?"),params![id],|r|r.get::<_,u64>(0)).unwrap(),0);
            }
        }
        // Simulate a kill while the retried sample is processing, without a
        // failure callback to clear its newly appended rows.
        {
            let conn=writer.lock_conn().unwrap();
            writer.write_feature_blobs(&conn,id,&[("primary_reads".into(),1,
                crate::blob::encode_dense_blob(&vec![2;100],crate::types::ValueScale::Raw,100))]).unwrap();
        }
        drop(writer);
        let writer=DbWriter::open(&path,&[],&[],&[],true).unwrap();
        assert_eq!(writer.insert_sample("sample","single-short",20,20,false).unwrap(),id);
        assert_eq!(writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Feature_blob",[],|r|r.get::<_,u64>(0)).unwrap(),0);
        writer.sample_complete("sample").unwrap();
        assert!(writer.insert_sample("sample","single-short",30,30,false).is_err());
    }
    #[test]
    fn encoding_failure_is_rejected_before_any_blob_row_is_written() {
        let directory=tempfile::tempdir().unwrap();
        let writer=DbWriter::create(&directory.path().join("bad.db"),&[ContigInfo{name:"a".into(),aliases:vec![],length:1,sequence:None}],&[],&[],true,false).unwrap();
        let id=writer.insert_sample("sample","single-short",1,1,false).unwrap();
        let blob=crate::blob::EncodedBlob{error:Some("injected compressor failure".into()),zoom:vec![],chunks:vec![]};
        let conn=writer.lock_conn().unwrap();
        assert!(writer.write_feature_blobs(&conn,id,&[("primary_reads".into(),1,blob)]).is_err());
        assert_eq!(conn.query_row("SELECT count(*) FROM Feature_blob",[],|r|r.get::<_,u64>(0)).unwrap(),0);
    }
    #[test]
    fn intercontig_curves_split_mag_hits_and_clear_stale_rows() {
        for mag_mode in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let contigs = ["a", "b", "c"].iter().map(|name| ContigInfo {
                name: name.to_string(), aliases: vec![], length: 20, sequence: None,
            }).collect::<Vec<_>>();
            let writer = DbWriter::create(&directory.path().join("hits.db"), &contigs, &[], &[], false, mag_mode).unwrap();
            if mag_mode {
                writer.write_mags(&[("one".into(), vec!["a".into(), "b".into()]), ("two".into(), vec!["c".into()])]).unwrap();
            }
            writer.lock_conn().unwrap().execute_batch("INSERT INTO Contig_blast_hits VALUES (1,2,1,10,1,10,990000),(1,3,1,10,1,10,980000),(1,1,1,5,6,10,1000000)").unwrap();
            writer.write_intercontig_hit_features().unwrap();
            for (name, expected) in [("hit_count", 1), ("hit_identity", 1), ("hit_count_within_mag", if mag_mode { 1 } else { 0 }), ("hit_identity_within_mag", if mag_mode { 1 } else { 0 })] {
                let count = writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Contig_blob WHERE Feature_id=?", params![feature_name_to_id(name).unwrap()], |r| r.get::<_, i64>(0)).unwrap();
                assert_eq!(count, expected, "{name} in MAG mode {mag_mode}");
            }
            {
                let conn = writer.lock_conn().unwrap();
                // Subject-only contigs must never acquire derived hit curves.
                assert_eq!(conn.query_row("SELECT count(*) FROM Contig_blob WHERE Contig_id IN (2,3)", [], |r| r.get::<_,i64>(0)).unwrap(), 0);
                let mut stmt = conn.prepare("SELECT CAST(Chunk_idx AS INTEGER),Data FROM Contig_blob_chunk WHERE Contig_id=1 AND Feature_id=? ORDER BY Chunk_idx").unwrap();
                let chunks = stmt.query_map(params![feature_name_to_id("hit_identity").unwrap()], |r| Ok((r.get::<_,i32>(0)?,r.get::<_,Vec<u8>>(1)?))).unwrap().collect::<duckdb::Result<Vec<_>>>().unwrap();
                let (_, values, metadata, _) = crate::blob::decode_sparse_from_chunks(&chunks).unwrap();
                assert_eq!(values.iter().copied().max(), Some(if mag_mode { 9800 } else { 9900 }));
                assert!(metadata.iter().all(|event| event.partner == Some(if mag_mode { 3 } else { 2 })));
            }
            writer.lock_conn().unwrap().execute("DELETE FROM Contig_blast_hits", []).unwrap();
            writer.write_intercontig_hit_features().unwrap();
            let remaining = writer.lock_conn().unwrap().query_row("SELECT count(*) FROM Contig_blob", [], |r| r.get::<_, i64>(0)).unwrap();
            assert_eq!(remaining, 0);
        }
    }

}
