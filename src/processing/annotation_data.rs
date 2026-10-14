//! Module-dependent annotation access with a bounded cache of active CDS batches.
use crate::features::CdsIndex;
use crate::types::FeatureAnnotation;
use anyhow::Result;
use duckdb::{params, Connection};
use std::collections::{HashMap, VecDeque};
use std::ops::Deref;
#[cfg(test)]
use std::path::Path;
use std::sync::{Arc, Mutex};

/// DNA is needed only when a mismatch can produce a codon annotation.
pub(super) fn needs_codon_annotations(
    arrays: &crate::features::FeatureArrays,
    threshold: f64,
) -> bool {
    arrays
        .mismatch_base_counts
        .iter()
        .zip(&arrays.primary_reads)
        .any(|(counts, &depth)| {
            let maximum = counts.iter().copied().max().unwrap_or(0);
            maximum > 0 && depth > 0 && !((maximum as f64 / depth as f64) < threshold)
        })
}

pub enum CdsHandle<'a> {
    Borrowed(&'a CdsIndex),
    Owned(Arc<CdsIndex>),
}
impl Deref for CdsHandle<'_> {
    type Target = CdsIndex;
    fn deref(&self) -> &CdsIndex {
        match self {
            Self::Borrowed(index) => index,
            Self::Owned(index) => index,
        }
    }
}
/// Existing in-memory callers and the calculation pipeline share the same interface.
pub trait CdsSource: Sync {
    fn lookup(&self, contig_id: i64) -> Option<CdsHandle<'_>>;
    fn check_error(&self) -> Result<()> {
        Ok(())
    }
}
impl CdsSource for HashMap<i64, CdsIndex> {
    fn lookup(&self, id: i64) -> Option<CdsHandle<'_>> {
        self.get(&id).map(CdsHandle::Borrowed)
    }
}
struct Batch {
    indexes: HashMap<i64, Arc<CdsIndex>>,
    bytes: usize,
}
struct State {
    batches: HashMap<i64, Batch>,
    order: VecDeque<i64>,
    bytes: usize,
    error: Option<String>,
    loads:u64, load_ns:u64, oversized_loads:u64,
}
pub(super) struct DatabaseCdsSource {
    state: Option<Mutex<State>>,
    connection: Option<Mutex<Connection>>,
}
impl DatabaseCdsSource {
    #[cfg(test)]
    pub fn new(database: &Path, enabled: bool) -> Result<Self> {
        if !enabled { return Ok(Self { state: None, connection: None }); }
        Self::from_connection(Connection::open(database)?, enabled)
    }

    pub fn from_connection(conn: Connection, enabled: bool) -> Result<Self> {
        if !enabled { return Ok(Self { state: None, connection: None }); }
        conn.prepare("SELECT Contig_id FROM Contig LIMIT 0")?;
        let exists: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name='Annotation_sequence')", [], |r| r.get(0))?;
        if !exists {
            return Ok(Self { state: None, connection:None });
        }
        Ok(Self {
            connection: Some(Mutex::new(conn)),
            state: Some(Mutex::new(State {
                batches: HashMap::new(),
                order: VecDeque::new(),
                bytes: 0,
                error: None, loads:0,load_ns:0,oversized_loads:0,
            })),
        })
    }
    pub(super) fn diagnostics(&self)->(u64,u64,u64,usize) {
        match &self.state {Some(state)=>{let state=state.lock().unwrap();(state.loads,state.load_ns,state.oversized_loads,state.bytes)},None=>(0,0,0,0)}
    }
    fn load(conn: &Connection, block: i64, requested_id:i64) -> Result<Batch> {
        // Contig IDs are assigned in import order. Loading a modest ID range
        // amortizes queries while avoiding a global DNA index.
        let first = block * 256 + 1;
        let mut stmt = conn.prepare(
            "WITH candidates AS MATERIALIZED (
                SELECT Annotation_id, Contig_id, \"Start\", \"End\", Strand
                FROM Contig_annotation_core
                WHERE Contig_id BETWEEN ? AND ? AND \"Type\"='CDS'
             )
             SELECT ca.Contig_id, ca.\"Start\", ca.\"End\", ca.Strand, aseq.Nucleotide_sequence
             FROM candidates ca JOIN (
                SELECT Annotation_id, Nucleotide_sequence FROM Annotation_sequence
                WHERE Annotation_id BETWEEN (SELECT min(Annotation_id) FROM candidates)
                                        AND (SELECT max(Annotation_id) FROM candidates)
             ) aseq USING (Annotation_id)
             WHERE aseq.Nucleotide_sequence IS NOT NULL ORDER BY ca.Annotation_id",
        )?;
        let rows = stmt.query_map(params![first, first + 255], |r| {
            Ok(FeatureAnnotation {
                contig_id: r.get(0)?,
                start: r.get(1)?,
                end: r.get(2)?,
                strand: r.get(3)?,
                feature_type: "CDS".into(),
                nucleotide_sequence: r.get(4)?,
                qualifiers: HashMap::new(),
                protein_sequence: None,
                segments: vec![],
                parent_key: None,
                self_key: None,
            })
        })?;
        let mut indexes:HashMap<i64,CdsIndex>=HashMap::new();
        let mut oversized=false;
        let mut sequence_bytes=0usize;
        for row in rows {
            let annotation=row?;
            if oversized && annotation.contig_id!=requested_id {continue;}
            sequence_bytes+=annotation.nucleotide_sequence.as_ref().map_or(0,|sequence|sequence.capacity())+std::mem::size_of::<crate::features::CdsInterval>();
            indexes.entry(annotation.contig_id).or_insert_with(CdsIndex::empty).push_annotation(annotation);
            if !oversized && sequence_bytes>128*1024*1024 {
                // Stop accumulating unrelated DNA as soon as the batch budget
                // is crossed. Only the requested contig remains live.
                indexes.retain(|id,_|*id==requested_id);
                oversized=true;
            }
        }
        for index in indexes.values_mut() {index.sort_intervals();}
        let indexes=indexes.into_iter().map(|(id,index)|(id,Arc::new(index))).collect::<HashMap<_,_>>();
        let bytes=indexes.values().map(|index|index.estimated_bytes()).sum::<usize>()+indexes.capacity()*64;
        // A partial block must never be cached as if all contigs were assessed.
        let bytes=if oversized {usize::MAX} else {bytes};
        Ok(Batch { indexes, bytes })
    }
}
impl CdsSource for DatabaseCdsSource {
    fn lookup(&self, id: i64) -> Option<CdsHandle<'_>> {
        let state = self.state.as_ref()?;
        let block=(id-1)/256;
        {
            let state=state.lock().unwrap();
            if state.error.is_some() {return None;}
            if let Some(batch)=state.batches.get(&block) {return batch.indexes.get(&id).cloned().map(CdsHandle::Owned);}
        }
        // Cache hits remain available while the independent SQL lock is held.
        let connection=self.connection.as_ref()?.lock().unwrap();
        {
            let state=state.lock().unwrap();
            if state.error.is_some() {return None;}
            if let Some(batch)=state.batches.get(&block) {return batch.indexes.get(&id).cloned().map(CdsHandle::Owned);}
        }
        let load_start=std::time::Instant::now();
        let loaded=Self::load(&connection,block,id);
        drop(connection);
        let mut state=state.lock().unwrap();
        state.loads+=1;
        state.load_ns+=load_start.elapsed().as_nanos() as u64;
        let batch=match loaded {
            Ok(batch)=>batch,
            Err(error)=>{state.error=Some(format!("Loading CDS batch for contig {id}: {error:#}"));return None;}
        };
        if batch.bytes==usize::MAX {state.oversized_loads+=1;}
        let result = batch.indexes.get(&id).cloned().map(CdsHandle::Owned);
        const CACHE_BYTES: usize = 128 * 1024 * 1024;
        const CACHE_BATCHES: usize = 64;
        while !state.order.is_empty()
            && (state.bytes.saturating_add(batch.bytes) > CACHE_BYTES || state.order.len() >= CACHE_BATCHES)
        {
            let old = state.order.pop_front().unwrap();
            if let Some(old) = state.batches.remove(&old) {
                state.bytes -= old.bytes;
            }
        }
        // An oversized batch is released after this lookup; only the current
        // contig's Arc remains owned by its calculation worker.
        if batch.bytes <= CACHE_BYTES {
            state.bytes += batch.bytes;
            state.order.push_back(block);
            state.batches.insert(block, batch);
        }
        result
    }
    fn check_error(&self) -> Result<()> {
        if let Some(state) = &self.state {
            let state = state.lock().map_err(|e| anyhow::anyhow!("{e}"))?;
            if let Some(error) = &state.error {
                anyhow::bail!("{error}");
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn database_batches_match_owned_cds_and_reuse_cached_dna() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("cds.db");
        let conn = Connection::open(&path).unwrap();
        conn.execute_batch("CREATE TABLE Contig (Contig_id BIGINT); CREATE TABLE Contig_annotation_core (Annotation_id BIGINT, Contig_id BIGINT, \"Start\" BIGINT, \"End\" BIGINT, Strand BIGINT, \"Type\" VARCHAR);
            CREATE TABLE Annotation_sequence (Annotation_id BIGINT, Nucleotide_sequence VARCHAR);
            INSERT INTO Contig_annotation_core VALUES (1,1,1,6,-1,'CDS'),(2,1,8,10,1,'gene'),(3,257,1,3,1,'CDS');
            INSERT INTO Annotation_sequence VALUES (1,'ATGAAA'),(2,'AAA'),(3,'ATG')").unwrap();
        let source = DatabaseCdsSource::new(&path, true).unwrap();
        let first = source.lookup(1).unwrap();
        let repeated = source.lookup(1).unwrap();
        let actual = first.find_cds(3).unwrap();
        assert_eq!(
            (
                actual.start,
                actual.end,
                actual.strand,
                actual.nucleotide_sequence.as_str()
            ),
            (1, 6, -1, "ATGAAA")
        );
        assert_eq!(
            actual.nucleotide_sequence.as_ptr(),
            repeated.find_cds(3).unwrap().nucleotide_sequence.as_ptr()
        );
        assert!(first.find_cds(9).is_none());
        assert!(source.lookup(2).is_none());
        assert_eq!(
            source
                .lookup(257)
                .unwrap()
                .find_cds(2)
                .unwrap()
                .nucleotide_sequence,
            "ATG"
        );
        source.check_error().unwrap();
        assert_eq!(
            source.state.as_ref().unwrap().lock().unwrap().batches.len(),
            2
        );
        let disabled = DatabaseCdsSource::new(&path, false).unwrap();
        assert!(disabled.state.is_none() && disabled.lookup(1).is_none());
    }
    #[test]
    fn failed_load_is_reported_instead_of_silently_omitting_codon_annotations() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("broken.db");
        Connection::open(&path)
            .unwrap()
            .execute_batch("CREATE TABLE Contig (Contig_id BIGINT); CREATE TABLE Annotation_sequence (Annotation_id BIGINT)")
            .unwrap();
        let source = DatabaseCdsSource::new(&path, true).unwrap();
        assert!(source.lookup(1).is_none());
        assert!(source.check_error().is_err());
    }
}
