//! On-demand terminal-repeat evidence, independent of general homology tracks.
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use crate::annotation_import::SequenceStore;
use crate::db::{DbWriter, RepeatsData};
use crate::processing_phage_packaging::{is_valid_terminal_repeat, PhageTerminiConfig};

type SearchResult = std::result::Result<Arc<Vec<RepeatsData>>, String>;
type CachedSearch = Arc<OnceLock<SearchResult>>;

/// Shared across samples. OnceLock avoids duplicate subprocesses for a contig
/// when several samples reach it concurrently. Empty successful results are cached.
#[derive(Default, Clone)]
pub struct TerminalTimings { pub search_ns: u64, pub dna_ns: u64, pub wait_ns: u64, pub parse_ns: u64 }
impl TerminalTimings {
    fn add(&mut self, other: &Self) { self.search_ns += other.search_ns; self.dna_ns += other.dna_ns; self.wait_ns += other.wait_ns; self.parse_ns += other.parse_ns; }
}

pub struct TerminalRepeatSearch {
    sequences: SequenceStore,
    indices: HashMap<i64, usize>,
    cache: Mutex<HashMap<i64, CachedSearch>>,
    unavailable: Mutex<Vec<i64>>,
    error: Mutex<Option<String>>,
    config: PhageTerminiConfig,
    completed: AtomicUsize,
    mode: Mutex<crate::processing::BlastMode>,
    general_compatible: Mutex<bool>,
    general_assessed: Mutex<std::collections::HashSet<i64>>,
    unassessed: Mutex<std::collections::HashSet<i64>>,
    timings: Mutex<HashMap<String, TerminalTimings>>,
    invalidate_persisted: Mutex<bool>,
}

pub(crate) fn has_candidate(
    starts: &[u64],
    ends: &[u64],
    coverage: &[u64],
    config: &PhageTerminiConfig,
) -> bool {
    if coverage.is_empty() {
        return false;
    }
    let aligned = coverage.iter().filter(|&&v| v > 0).count();
    if aligned as f64 * 100.0 / (coverage.len() as f64) < config.min_aligned_fraction as f64 {
        return false;
    }
    // Same local thresholds as the first phase of peak detection, BEFORE
    // clipping, Poisson, DTR, or packaging filters. Peaks need not be at an end.
    starts
        .iter()
        .zip(coverage)
        .chain(ends.iter().zip(coverage))
        .any(|(&count, &cov)| {
            count > 0
                && count >= config.min_events as u64
                && count >= (cov as f64 * config.min_frequency) as u64
        })
}

impl TerminalRepeatSearch {
    pub(crate) fn new(sequences: SequenceStore, ids: Vec<i64>, config: PhageTerminiConfig) -> Self {
        Self {
            sequences,
            indices: ids.into_iter().enumerate().map(|(i, id)| (id, i)).collect(),
            cache: Mutex::new(HashMap::new()),
            unavailable: Mutex::new(Vec::new()),
            error: Mutex::new(None),
            config,
            completed: AtomicUsize::new(0),
            mode: Mutex::new(crate::processing::BlastMode::Termini),
            general_compatible: Mutex::new(true),
            general_assessed: Mutex::new(std::collections::HashSet::new()),
            unassessed: Mutex::new(std::collections::HashSet::new()),
            timings: Mutex::new(HashMap::new()), invalidate_persisted: Mutex::new(false),
        }
    }

    fn parameters(&self) -> String { format!("{};initial_window=10000;min_identity={};max_distance={}", crate::processing::BLAST_PARAMETERS, self.config.min_identity_dtr, self.config.max_distance_duplication) }
    pub(crate) fn configure_mode(&self, mode: crate::processing::BlastMode, writer: &DbWriter) -> Result<()> {
        *self.mode.lock().unwrap() = mode;
        let conn = writer.lock_conn()?;
        let metadata = |key: &str| -> String { conn.query_row("SELECT Value FROM Database_metadata WHERE Key=?", duckdb::params![key], |r| r.get(0)).unwrap_or_default() };
        *self.general_compatible.lock().unwrap() = metadata("Blast_parameters") == crate::processing::BLAST_PARAMETERS;
        let exists: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name='Terminal_repeat_search')", [], |r| r.get(0))?;
        if !exists { return Ok(()); }
        if metadata("Terminal_repeat_parameters") != self.parameters() { *self.invalidate_persisted.lock().unwrap() = true; return Ok(()); }
        let mut hits_by_id: HashMap<i64, Vec<RepeatsData>> = HashMap::new();
        let mut stmt = conn.prepare("SELECT e.Contig_id,c.Contig_name,e.Position1,e.Position2,e.Position1prime,e.Position2prime,e.Pident,e.Is_direct FROM Terminal_repeat_evidence e JOIN Contig c USING (Contig_id)")?;
        for row in stmt.query_map([], |r| Ok((r.get::<_,i64>(0)?, RepeatsData {
            contig_name:r.get(1)?, position1:r.get(2)?,position2:r.get(3)?,position1prime:r.get(4)?,position2prime:r.get(5)?,pident:r.get(6)?,is_direct:r.get(7)? })))? {
            let (id,hit)=row?; hits_by_id.entry(id).or_default().push(hit);
        }
        let mut stmt = conn.prepare("SELECT Contig_id,Status FROM Terminal_repeat_search WHERE Status IN ('found','none')")?;
        let mut cache = self.cache.lock().unwrap();
        for row in stmt.query_map([],|r| Ok((r.get::<_,i64>(0)?,r.get::<_,String>(1)?)))? {
            let (id,_)=row?; let cell = Arc::new(OnceLock::new());
            let _ = cell.set(Ok(Arc::new(hits_by_id.remove(&id).unwrap_or_default()))); cache.insert(id,cell);
        }
        Ok(())
    }
    pub(crate) fn mark_general_assessed(&self, ids: &[i64]) { self.general_assessed.lock().unwrap().extend(ids.iter().copied()); let mut cache=self.cache.lock().unwrap(); for id in ids {cache.remove(id);} }
    pub(crate) fn finish_sample(&self, sample: &str) -> TerminalTimings {
        self.timings.lock().unwrap().remove(sample).unwrap_or_default()
    }
    #[cfg(test)]
    pub(crate) fn lookup(&self, id:i64,name:&str,existing:&[RepeatsData]) -> Arc<Vec<RepeatsData>> { self.lookup_for_sample(id,name,existing,"test") }

    pub(crate) fn lookup_for_sample(
        &self,
        id: i64,
        name: &str,
        existing: &[RepeatsData],
        sample: &str,
    ) -> Arc<Vec<RepeatsData>> {
        if self.error.lock().unwrap().is_some() {
            return Arc::new(Vec::new());
        }
        let cache_start = std::time::Instant::now();
        let cell = self
            .cache
            .lock()
            .unwrap()
            .entry(id)
            .or_insert_with(|| Arc::new(OnceLock::new()))
            .clone();
        let wait_start = std::time::Instant::now();
        let mut phases = TerminalTimings { wait_ns: cache_start.elapsed().as_nanos() as u64, ..Default::default() };
        let mut owns_search = false;
        let result = cell.get_or_init(|| {
            owns_search = true;
            let assessed = self.general_assessed.lock().unwrap().contains(&id);
            if !existing.is_empty() && (assessed || *self.general_compatible.lock().unwrap()) { return Ok(Arc::new(existing.to_vec())); }
            if assessed || *self.general_compatible.lock().unwrap() {
                let started=std::time::Instant::now();
                let loaded=(||->Result<Vec<RepeatsData>> {
                    let index=*self.indices.get(&id).context("Missing reference index")?;
                    let length=self.sequences.reference_length(index)?;
                    Ok(self.sequences.general_repeats(index,name)?.into_iter().filter(|hit|hit.pident>=self.config.min_identity_dtr as f64 && is_valid_terminal_repeat(hit,length,self.config.max_distance_duplication)).collect())
                })();
                phases.parse_ns+=started.elapsed().as_nanos() as u64;
                match loaded {Ok(hits) if assessed || !hits.is_empty()=>return Ok(Arc::new(hits)),Ok(_)=>{},Err(error)=>return Err(error.to_string())}
            }
            if *self.mode.lock().unwrap() != crate::processing::BlastMode::Termini {
                if !assessed { self.unassessed.lock().unwrap().insert(id); }
                return Ok(Arc::new(Vec::new()));
            }
            let result = (|| -> Result<Vec<RepeatsData>> {
                let index = *self.indices.get(&id).context("Missing contig sequence index")?;
                let dna_start = std::time::Instant::now();
                let sequence = self.sequences.load_sequence(index)?;
                phases.dna_ns += dna_start.elapsed().as_nanos() as u64;
                let Some(sequence) = sequence.filter(|s| !s.is_empty()) else {
                    let mut missing = self.unavailable.lock().unwrap();
                    missing.push(id);
                    if missing.len() == 1 {
                        eprintln!("WARNING: Sample '{sample}', contig '{name}': terminal-repeat evidence unavailable without reference DNA; read-based termini analysis will continue. Supply annotations with sequences or --assembly.");
                    }
                    return Ok(Vec::new());
                };
                let hits = search_ends_timed(&sequence, name, &self.config, sample, &mut phases)?;
                self.completed.fetch_add(1, Ordering::Relaxed);
                Ok(hits)
            })();
            result.map(Arc::new).map_err(|e| format!("Sample '{sample}', contig '{name}': terminal-repeat search failed: {e:#}"))
        });
        if !owns_search { phases.wait_ns += wait_start.elapsed().as_nanos() as u64; }
        self.timings.lock().unwrap().entry(sample.to_owned()).or_default().add(&phases);
        match result {
            Ok(hits) => hits.clone(),
            Err(message) => {
                self.error
                    .lock()
                    .unwrap()
                    .get_or_insert_with(|| message.clone());
                Arc::new(Vec::new())
            }
        }
    }

    /// Call before returning a sample so a failed search never produces a saved
    /// classification based on silently missing evidence.
    pub(crate) fn check_error(&self) -> Result<()> {
        if let Some(message) = self.error.lock().unwrap().as_ref() {
            anyhow::bail!("{message}");
        }
        Ok(())
    }

    pub(crate) fn persist(&self, writer: &DbWriter) -> Result<()> {
        self.check_error()?;
        let mut conn = writer.lock_conn()?;
        let tx = conn.transaction()?;
        tx.execute_batch("CREATE TABLE IF NOT EXISTS Terminal_repeat_search (Contig_id BIGINT PRIMARY KEY, Status VARCHAR);
            CREATE TABLE IF NOT EXISTS Terminal_repeat_evidence (Contig_id BIGINT, Position1 INTEGER, Position2 INTEGER, Position1prime INTEGER, Position2prime INTEGER, Pident DOUBLE, Is_direct BOOLEAN);
            CREATE TEMP TABLE Terminal_search_batch (Contig_id BIGINT, Status VARCHAR);
            INSERT OR REPLACE INTO Database_metadata VALUES ('Terminal_repeat_mode', 'on-demand')")?;
        if *self.invalidate_persisted.lock().unwrap() {
            tx.execute_batch("DELETE FROM Terminal_repeat_evidence; DELETE FROM Terminal_repeat_search")?;
        }
        tx.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Terminal_repeat_parameters', ?)", duckdb::params![self.parameters()])?;
        tx.execute("INSERT OR REPLACE INTO Database_metadata VALUES ('Terminal_repeat_mode', ?)", duckdb::params![self.mode.lock().unwrap().as_str()])?;
        let unassessed = self.unassessed.lock().unwrap();
        let missing: std::collections::HashSet<i64> =
            self.unavailable.lock().unwrap().iter().copied().collect();
        let cache = self.cache.lock().unwrap();
        {
            let mut appender = tx.appender("Terminal_search_batch")?;
            for (&id, cell) in cache.iter() {
                if let Some(Ok(hits)) = cell.get() {
                    let status = if unassessed.contains(&id) { "not_assessed" } else if missing.contains(&id) {
                        "unavailable"
                    } else if hits.is_empty() {
                        "none"
                    } else {
                        "found"
                    };
                    appender.append_row(duckdb::params![id, status])?;
                }
            }
            appender.flush()?;
        }
        tx.execute_batch(
            "DELETE FROM Terminal_repeat_evidence USING Terminal_search_batch b
                WHERE Terminal_repeat_evidence.Contig_id = b.Contig_id;
            INSERT OR REPLACE INTO Terminal_repeat_search SELECT * FROM Terminal_search_batch;
            DROP TABLE Terminal_search_batch",
        )?;
        {
            let mut appender = tx.appender("Terminal_repeat_evidence")?;
            for (&id, cell) in cache.iter() {
                if let Some(Ok(hits)) = cell.get() {
                    for hit in hits.iter() {
                        appender.append_row(duckdb::params![
                            id,
                            hit.position1,
                            hit.position2,
                            hit.position1prime,
                            hit.position2prime,
                            hit.pident,
                            hit.is_direct
                        ])?;
                    }
                }
            }
            appender.flush()?;
        }
        tx.commit()?;
        Ok(())
    }
}

fn write_fasta(path: &Path, sequence: &[u8]) -> Result<()> {
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    file.write_all(b">terminal_window\n")?;
    for chunk in sequence.chunks(80) {
        file.write_all(chunk)?;
        file.write_all(b"\n")?;
    }
    file.flush()?;
    Ok(())
}

fn parse_hit(
    line: &str,
    offset: usize,
    length: usize,
    name: &str,
    config: &PhageTerminiConfig,
) -> Result<Option<RepeatsData>> {
    let fields: Vec<_> = line.split_whitespace().collect();
    anyhow::ensure!(fields.len() == 5, "Invalid terminal BLAST row: {line}");
    let hit = RepeatsData {
        contig_name: name.into(),
        position1: fields[0].parse()?,
        position2: fields[1].parse()?,
        position1prime: fields[2].parse::<i32>()? + offset as i32,
        position2prime: fields[3].parse::<i32>()? + offset as i32,
        pident: fields[4].parse()?,
        is_direct: (fields[0].parse::<i32>()? < fields[1].parse::<i32>()?)
            == (fields[2].parse::<i32>()? < fields[3].parse::<i32>()?),
    };
    if (hit.position1 == hit.position1prime && hit.position2 == hit.position2prime)
        || hit.pident < config.min_identity_dtr as f64
        || !is_valid_terminal_repeat(&hit, length, config.max_distance_duplication)
    {
        return Ok(None);
    }
    Ok(Some(hit))
}

#[cfg(test)]
fn search_ends(
    sequence: &[u8],
    name: &str,
    config: &PhageTerminiConfig,
) -> Result<Vec<RepeatsData>> {
    search_ends_timed(sequence,name,config,"test",&mut TerminalTimings::default())
}

fn search_ends_timed(sequence:&[u8],name:&str,config:&PhageTerminiConfig,_sample:&str,phases:&mut TerminalTimings) -> Result<Vec<RepeatsData>> {
    if sequence.is_empty() {
        return Ok(Vec::new());
    }
    anyhow::ensure!(
        sequence.len() <= i32::MAX as usize,
        "Reference exceeds supported coordinate range"
    );
    let directory = tempfile::Builder::new()
        .prefix("thebigbam-termini-")
        .tempdir()?;
    let query = directory.path().join("left.fa");
    let subject = directory.path().join("right.fa");
    let output = directory.path().join("hits.tsv");
    let mut window = 10_000.min(sequence.len());
    loop {
        let offset = sequence.len() - window;
        write_fasta(&query, &sequence[..window])?;
        write_fasta(&subject, &sequence[offset..])?;
        let search_start = std::time::Instant::now();
        let stderr_path = directory.path().join("stderr.log");
        let mut child = std::process::Command::new("blastn")
            .arg("-query").arg(&query).arg("-subject").arg(&subject)
            .args(crate::processing::BLAST_ARGS)
            .args(["-outfmt", "6 qstart qend sstart send pident"])
            .arg("-out").arg(&output)
            .stdout(std::process::Stdio::null())
            .stderr(std::fs::File::create(&stderr_path)?)
            .spawn().context("Starting blastn for Phagetermini terminal search")?;
        let status = child.wait().context("Waiting for Phagetermini terminal blastn")?;
        phases.search_ns += search_start.elapsed().as_nanos() as u64;
        anyhow::ensure!(status.success(), "blastn failed: {}", std::fs::read_to_string(stderr_path).unwrap_or_default());
        let parse_start = std::time::Instant::now();
        let mut hits = Vec::new();
        for line in BufReader::new(std::fs::File::open(&output)?).lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            if let Some(hit) = parse_hit(&line, offset, sequence.len(), name, config)? {
                hits.push(hit);
            }
        }
        phases.parse_ns += parse_start.elapsed().as_nanos() as u64;
        let touches_boundary = hits.iter().any(|h| {
            h.position1.max(h.position2) as usize >= window.saturating_sub(32)
                || h.position1prime.min(h.position2prime) as usize <= offset + 32
        });
        if window == sequence.len() || (!hits.is_empty() && !touches_boundary) {
            return Ok(hits);
        }
        // Expand even when no match is found: two short end windows can sample
        // different portions of a long repeat. Fall back to the whole contig.
        window = window.saturating_mul(2).min(sequence.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ContigInfo;

    #[test]
    fn candidates_use_local_evidence_before_repeat_dependent_filters() {
        let config = PhageTerminiConfig::default();
        let coverage = vec![100; 1000];
        let mut starts = vec![0; 1000];
        assert!(!has_candidate(&starts, &starts, &coverage, &config));
        starts[500] = config.min_events as u64;
        assert!(has_candidate(&starts, &[], &coverage, &config));
        assert!(has_candidate(&[], &starts, &coverage, &config));
        starts[500] -= 1;
        assert!(!has_candidate(&starts, &[], &coverage, &config));
        starts[500] = 100;
        assert!(!has_candidate(&starts, &[], &vec![0; 1000], &config));
    }

    #[test]
    fn terminal_coordinates_restore_offset_and_both_orientations() {
        let config = PhageTerminiConfig::default();
        let direct = parse_hit("1 500 1501 2000 99", 8000, 10000, "c", &config)
            .unwrap()
            .unwrap();
        assert_eq!(
            (direct.position1prime, direct.position2prime),
            (9501, 10000)
        );
        assert!(direct.is_direct);
        let inverted = parse_hit("1 500 2000 1501 99", 8000, 10000, "c", &config)
            .unwrap()
            .unwrap();
        assert!(!inverted.is_direct);
        assert_eq!(
            (inverted.position1prime, inverted.position2prime),
            (10000, 9501)
        );
        assert!(parse_hit("1 10000 1 10000 100", 0, 10000, "c", &config)
            .unwrap()
            .is_none());
        assert!(parse_hit("1 500 1501 2000 80", 8000, 10000, "c", &config)
            .unwrap()
            .is_none());
        assert!(
            parse_hit("501 1000 1501 2000 99", 8000, 10000, "c", &config)
                .unwrap()
                .is_none()
        );
        assert!(parse_hit("bad row", 0, 10000, "c", &config).is_err());
    }

    #[test]
    fn missing_sequence_and_existing_evidence_are_cached_across_workers() {
        use rayon::prelude::*;
        let directory = tempfile::tempdir().unwrap();
        let mut store = SequenceStore::new(directory.path()).unwrap();
        store
            .push_batch(&mut [ContigInfo {
                name: "c".into(),
                aliases: vec![],
                length: 10000,
                sequence: None,
            }])
            .unwrap();
        let search = TerminalRepeatSearch::new(store, vec![1], PhageTerminiConfig::default());
        (0..100)
            .into_par_iter()
            .for_each(|_| assert!(search.lookup(1, "c", &[]).is_empty()));
        assert_eq!(search.unavailable.lock().unwrap().len(), 1);
        assert_eq!(search.cache.lock().unwrap().len(), 1);
        search.check_error().unwrap();
        let hit = parse_hit(
            "1 500 9501 10000 99",
            0,
            10000,
            "existing",
            &PhageTerminiConfig::default(),
        )
        .unwrap()
        .unwrap();
        let first = search.lookup(2, "existing", &[hit]);
        let second = search.lookup(2, "existing", &[]);
        assert!(Arc::ptr_eq(&first, &second));
        assert_eq!(second.len(), 1);
        assert_eq!(search.completed.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn evidence_persistence_keeps_unavailable_distinct_and_is_idempotent() {
        let directory = tempfile::tempdir().unwrap();
        let mut contigs = vec![ContigInfo {
            name: "c".into(),
            aliases: vec![],
            length: 10000,
            sequence: None,
        }];
        let writer = DbWriter::create(
            &directory.path().join("evidence.db"),
            &contigs,
            &[],
            &[],
            false,
            false,
        )
        .unwrap();
        let mut store = SequenceStore::new(directory.path()).unwrap();
        store.push_batch(&mut contigs).unwrap();
        let search = TerminalRepeatSearch::new(store, vec![1], PhageTerminiConfig::default());
        search.lookup(1, "c", &[]);
        search.persist(&writer).unwrap();
        search.persist(&writer).unwrap();
        let conn = writer.lock_conn().unwrap();
        assert_eq!(
            conn.query_row(
                "SELECT Status FROM Terminal_repeat_search WHERE Contig_id = 1",
                [],
                |r| r.get::<_, String>(0)
            )
            .unwrap(),
            "unavailable"
        );
        assert_eq!(
            conn.query_row("SELECT count(*) FROM Terminal_repeat_evidence", [], |r| r
                .get::<_, i64>(
                0
            ))
            .unwrap(),
            0
        );
        assert_eq!(
            conn.query_row("SELECT count(*) FROM Contig_directRepeats", [], |r| r
                .get::<_, i64>(0))
                .unwrap(),
            0
        );
    }

    #[test]
    #[ignore = "requires BLAST+ on PATH"]
    fn adaptive_search_finds_repeat_longer_than_initial_window() {
        // Deterministic non-low-complexity DNA without another dependency.
        let mut state = 1234567u64;
        let mut dna = |len| {
            (0..len)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((state >> 32) & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let repeat = dna(15_000);
        let mut sequence = repeat.clone();
        sequence.extend(dna(8000));
        sequence.extend_from_slice(&repeat);
        let hits = search_ends(&sequence, "long_dtr", &PhageTerminiConfig::default()).unwrap();
        assert!(hits
            .iter()
            .any(|h| h.is_direct && (h.position2 - h.position1).abs() >= 14_999));
    }
}

#[cfg(test)]
mod mode_tests {
    use super::*;
    #[test]
    fn none_and_all_do_not_launch_candidate_searches_and_preserve_assessment_status() {
        let directory = tempfile::tempdir().unwrap();
        for mode in [crate::processing::BlastMode::None,crate::processing::BlastMode::Mags,crate::processing::BlastMode::All] {
            let store = SequenceStore::new(directory.path()).unwrap();
            let search = TerminalRepeatSearch::new(store,vec![1],PhageTerminiConfig::default());
            *search.mode.lock().unwrap() = mode;
            assert!(search.lookup_for_sample(1,"a",&[],"sample-a").is_empty());
            assert!(search.unassessed.lock().unwrap().contains(&1));
            assert!(search.error.lock().unwrap().is_none());
            assert_eq!(search.completed.load(Ordering::Relaxed),0);
            let phases = search.finish_sample("sample-a");
            assert_eq!((phases.search_ns,phases.dna_ns,phases.parse_ns),(0,0,0));
        }
        let search = TerminalRepeatSearch::new(SequenceStore::new(directory.path()).unwrap(),vec![1],PhageTerminiConfig::default());
        *search.mode.lock().unwrap() = crate::processing::BlastMode::All;
        search.mark_general_assessed(&[1]);
        assert!(search.lookup_for_sample(1,"a",&[],"sample-b").is_empty());
        assert!(!search.unassessed.lock().unwrap().contains(&1));
    }
}
