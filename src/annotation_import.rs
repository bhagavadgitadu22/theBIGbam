//! File-at-a-time annotation ingestion and disk-backed preprocessing sequences.
//! Persisted annotations are loaded in bounded batches during BAM processing.
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::Write;
#[cfg(not(unix))]
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use crate::db::DbWriter;
use crate::parser::{compute_annotation_sequences, parse_annotations, parse_fasta};
use crate::processing::{ProcessConfig, ViewMode};
use crate::types::{ContigInfo, FeatureAnnotation};

/// A sequential temporary sequence file with a small in-memory offset index.
/// Positional reads let BLAST workers load sequences concurrently on Unix.
pub(crate) struct SequenceStore {
    file: tempfile::NamedTempFile,
    entries: Vec<Option<(u64, usize)>>,
    end: u64,
    database: Option<std::sync::Mutex<duckdb::Connection>>,
    database_ids: Vec<i64>,
}

impl SequenceStore {
    pub(crate) fn new(directory: &Path) -> Result<Self> {
        Ok(Self {
            file: tempfile::Builder::new()
                .prefix(".thebigbam-sequences-")
                .tempfile_in(directory)?,
            entries: Vec::new(),
            end: 0,
            database: None, database_ids: Vec::new(),
        })
    }

    /// Candidate-only sequence reads; no assembly DNA is retained between searches.
    /// Share the writer's database instance rather than independently reopening its path.
    pub(crate) fn from_connection(directory: &Path, conn: duckdb::Connection, ids: Vec<i64>) -> Result<Self> {
        // Validate eagerly so a catalog error stops before processing any samples.
        conn.prepare("SELECT Contig_length FROM Contig WHERE Contig_id=?")?;
        conn.prepare("SELECT Sequence FROM Contig_sequence WHERE Contig_id=?")?;
        let mut store = Self::new(directory)?;
        store.database = Some(std::sync::Mutex::new(conn));
        store.database_ids = ids;
        Ok(store)
    }

    pub(crate) fn push_batch(&mut self, contigs: &mut [ContigInfo]) -> Result<()> {
        let mut writer = std::io::BufWriter::new(self.file.as_file_mut());
        for contig in contigs {
            let entry = if let Some(sequence) = contig.sequence.take() {
                writer
                    .write_all(&sequence)
                    .context("Writing temporary contig sequences (check free disk space)")?;
                let entry = (self.end, sequence.len());
                self.end += sequence.len() as u64;
                Some(entry)
            } else {
                None
            };
            self.entries.push(entry);
        }
        writer.flush()?;
        Ok(())
    }

    pub(crate) fn has_sequence(&self, index: usize) -> bool {
        self.entries[index].is_some()
    }

    pub(crate) fn load(&self, index: usize, contig: &ContigInfo) -> Result<ContigInfo> {
        let mut loaded = contig.clone();
        loaded.sequence = self.load_sequence(index)?;
        Ok(loaded)
    }

    pub(crate) fn general_repeats(&self,index:usize,name:&str)->Result<Vec<crate::db::RepeatsData>> {
        let Some(database)=&self.database else {return Ok(Vec::new());};
        let id=self.database_ids.get(index).context("Missing reference ID")?;
        let conn=database.lock().map_err(|e|anyhow::anyhow!("{e}"))?;
        let mut repeats=Vec::new();
        for (table,direct) in [("Contig_directRepeats",true),("Contig_invertedRepeats",false)] {
            let mut stmt=conn.prepare(&format!("SELECT Position1,Position2,Position1prime,Position2prime,Pident FROM {table} WHERE Contig_id=?"))?;
            for row in stmt.query_map(duckdb::params![id],|r|Ok(crate::db::RepeatsData{contig_name:name.to_owned(),position1:r.get(0)?,position2:r.get(1)?,position1prime:r.get(2)?,position2prime:r.get(3)?,pident:r.get::<_,i32>(4)? as f64/crate::types::get_column_scale(table,"Pident"),is_direct:direct}))? {repeats.push(row?);}
        }
        Ok(repeats)
    }
    pub(crate) fn reference_length(&self,index:usize)->Result<usize> {
        if let Some(database)=&self.database {
            return Ok(database.lock().map_err(|e|anyhow::anyhow!("{e}"))?.query_row("SELECT Contig_length FROM Contig WHERE Contig_id=?",duckdb::params![self.database_ids[index]],|r|r.get::<_,i64>(0))? as usize);
        }
        Ok(self.entries.get(index).and_then(|entry|*entry).map_or(0,|(_,length)|length))
    }

    pub(crate) fn load_sequence(&self, index: usize) -> Result<Option<Vec<u8>>> {
        if let Some(database) = &self.database {
            let id = self.database_ids.get(index).context("Missing database sequence ID")?;
            let conn = database.lock().map_err(|e| anyhow::anyhow!("{e}"))?;
            let mut stmt = conn.prepare("SELECT Sequence FROM Contig_sequence WHERE Contig_id = ?")?;
            let mut rows = stmt.query(duckdb::params![id])?;
            return match rows.next()? {
                Some(row) => Ok(row.get::<_, Option<String>>(0)?.map(String::into_bytes)),
                None => Ok(None),
            };
        }
        if let Some((offset, len)) = self.entries[index] {
            let mut sequence = vec![0; len];
            #[cfg(unix)]
            {
                use std::os::unix::fs::FileExt;
                self.file.as_file().read_exact_at(&mut sequence, offset)?;
            }
            #[cfg(not(unix))]
            {
                let mut file = self.file.reopen()?;
                file.seek(SeekFrom::Start(offset))?;
                file.read_exact(&mut sequence)?;
            }
            Ok(Some(sequence))
        } else {
            Ok(None)
        }
    }
}

pub(crate) struct ImportedAnnotations {
    pub writer: DbWriter,
    pub contigs: Vec<ContigInfo>,
    pub cds: Vec<FeatureAnnotation>,
    pub sequences: SequenceStore,
    pub mags: Vec<(String, Vec<String>)>,
    pub annotation_count: usize,
}

/// Streaming is available for per-MAG inputs and annotation directories with
/// embedded sequences. The global annotation + separate assembly fallback keeps
/// its existing matching behavior.
pub(crate) fn can_stream(genbank: &Path, assembly: &Path, config: &ProcessConfig) -> bool {
    if config.view_mode == ViewMode::Mag {
        !genbank.is_file()
            && !config.mag_manifest.is_empty()
            && config.mag_manifest.iter().any(|m| m.gb_path.is_some())
    } else {
        genbank.is_dir() && assembly.as_os_str().is_empty()
    }
}

pub(crate) fn import_annotations(
    genbank: &Path,
    output: &Path,
    config: &ProcessConfig,
    has_bam: bool,
    keep_cds: bool,
    mut progress: impl FnMut(usize, usize, usize),
) -> Result<ImportedAnnotations> {
    let directory = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::create_dir_all(directory)?;
    let mut imported = ImportedAnnotations {
        writer: DbWriter::create_streaming(output, has_bam, config.view_mode == ViewMode::Mag)?,
        contigs: Vec::new(),
        cds: Vec::new(),
        sequences: SequenceStore::new(directory)?,
        mags: Vec::new(),
        annotation_count: 0,
    };
    if config.enable_timing { imported.writer.enable_size_report(output)?; }
    let inputs = if config.view_mode == ViewMode::Mag {
        config
            .mag_manifest
            .iter()
            .map(|m| (Some(m.name.clone()), m.gb_path.clone(), m.asm_path.clone()))
            .collect::<Vec<_>>()
    } else {
        let mut paths = std::fs::read_dir(genbank)?
            .map(|e| e.map(|e| e.path()))
            .collect::<std::io::Result<Vec<_>>>()?;
        paths.retain(|p| p.is_file() && crate::parser::is_annotation_file(p));
        paths.sort();
        anyhow::ensure!(
            !paths.is_empty(),
            "No annotation files found in {}",
            genbank.display()
        );
        paths.into_iter().map(|p| (None, Some(p), None)).collect()
    };
    let total = inputs.len();
    for (batch, (mag, annotation_path, assembly_path)) in inputs.into_iter().enumerate() {
        let (mut contigs, mut annotations, mut qualifiers) = match annotation_path {
            Some(path) => {
                parse_annotations(&path).with_context(|| format!("Parsing {}", path.display()))?
            }
            None => (Vec::new(), Vec::new(), Vec::new()),
        };
        if let Some(assembly) = assembly_path {
            let by_name: HashMap<String, usize> = contigs
                .iter()
                .enumerate()
                .map(|(i, c)| (c.name.clone(), i))
                .collect();
            for (name, sequence) in parse_fasta(&assembly)? {
                if let Some(&index) = by_name.get(&name) {
                    if contigs[index].sequence.is_none() {
                        contigs[index].sequence = Some(sequence);
                    }
                } else {
                    contigs.push(ContigInfo {
                        name,
                        aliases: vec![],
                        length: sequence.len(),
                        sequence: Some(sequence),
                    });
                }
            }
        }
        compute_annotation_sequences(&contigs, &mut annotations);
        let offset = imported.contigs.len() as i64;
        for ann in &mut annotations {
            ann.contig_id += offset;
        }
        for (id, _) in &mut qualifiers {
            *id += offset;
        }
        imported.writer.append_annotation_batch(
            &contigs,
            &annotations,
            &qualifiers,
            offset,
            imported.annotation_count as i64,
        )?;
        imported.annotation_count += annotations.len();
        if let Some(name) = mag {
            imported
                .mags
                .push((name, contigs.iter().map(|c| c.name.clone()).collect()));
        }
        // Persisted qualifiers and proteins have no consumers in BAM processing.
        if keep_cds {
            imported
                .cds
                .extend(annotations.into_iter().filter_map(|mut ann| {
                    if ann.feature_type != "CDS" || ann.nucleotide_sequence.is_none() {
                        return None;
                    }
                    ann.qualifiers = HashMap::new();
                    ann.protein_sequence = None;
                    ann.parent_key = None;
                    ann.self_key = None;
                    ann.segments = Vec::new();
                    Some(ann)
                }));
        }
        imported.sequences.push_batch(&mut contigs)?;
        imported.contigs.extend(contigs);
        // Progress at file boundaries; avoid forcing a checkpoint per small MAG.
        if (batch + 1) % 100 == 0 || batch + 1 == total {
            progress(batch + 1, total, imported.annotation_count);
        }
    }
    imported.writer.finish_annotation_import()?;
    eprintln!(
        "Annotation import complete: {} contigs, {} annotations",
        imported.contigs.len(),
        imported.annotation_count
    );
    Ok(imported)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_store_releases_dna_and_supports_out_of_order_reads() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = SequenceStore::new(directory.path()).unwrap();
        let mut contigs = vec![
            ContigInfo {
                name: "a".into(),
                aliases: vec!["alias".into()],
                length: 6,
                sequence: Some(b"ATGCAA".to_vec()),
            },
            ContigInfo {
                name: "b".into(),
                aliases: vec![],
                length: 3,
                sequence: None,
            },
            ContigInfo {
                name: "c".into(),
                aliases: vec![],
                length: 3,
                sequence: Some(b"GCC".to_vec()),
            },
        ];
        store.push_batch(&mut contigs[..2]).unwrap();
        store.push_batch(&mut contigs[2..]).unwrap();
        assert!(contigs.iter().all(|c| c.sequence.is_none()));
        assert!(!store.has_sequence(1));
        assert_eq!(
            store.load(2, &contigs[2]).unwrap().sequence.unwrap(),
            b"GCC"
        );
        let first = store.load(0, &contigs[0]).unwrap();
        assert_eq!(first.sequence.unwrap(), b"ATGCAA");
        assert_eq!(first.aliases, vec!["alias"]);
        use rayon::prelude::*;
        (0..100).into_par_iter().for_each(|i| {
            let index = if i % 2 == 0 { 0 } else { 2 };
            let expected: &[u8] = if index == 0 { b"ATGCAA" } else { b"GCC" };
            assert_eq!(
                store
                    .load(index, &contigs[index])
                    .unwrap()
                    .sequence
                    .unwrap(),
                expected
            );
        });
        drop(store);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    }
    fn config(view_mode: ViewMode) -> ProcessConfig {
        ProcessConfig {
            threads: 2,
            min_aligned_fraction: 0.0,
            min_coverage_depth: 0.0,
            bar_ratio: 0.0,
            sequencing_type: None,
            phagetermini_config: Default::default(),
            gc_params: Default::default(),
            variation_percentage: 0.0,
            min_occurrences: 2,
            enable_timing: true,
            blast: crate::processing::BlastMode::None,
            sample_name: String::new(),
            primary_reads_seen: std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0)),
            terminal_repeats: None,
            view_mode,
            mag_manifest: Vec::new(),
            contig_drops_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
        }
    }

    fn fixture(directory: &Path, name: &str) -> std::path::PathBuf {
        let path = directory.join(format!("{name}.gff3"));
        std::fs::write(&path, format!(
            "##gff-version 3\n##sequence-region {name} 1 12\n{name}\ttest\tCDS\t1\t9\t.\t+\t0\tID={name}_cds;product=example;locus_tag={name}_gene\n##FASTA\n>{name}\nATGAAATAAGCC\n"
        )).unwrap();
        path
    }

    #[test]
    fn file_import_retains_only_required_cds_and_preserves_sequences() {
        let directory = tempfile::tempdir().unwrap();
        let inputs = directory.path().join("inputs");
        std::fs::create_dir(&inputs).unwrap();
        fixture(&inputs, "a");
        fixture(&inputs, "b");
        let config = config(ViewMode::Contig);
        let imported = import_annotations(
            &inputs,
            &directory.path().join("with_cds.db"),
            &config,
            true,
            true,
            |_, _, _| {},
        )
        .unwrap();
        assert_eq!(imported.annotation_count, 2);
        assert_eq!(imported.cds.len(), 2);
        assert_eq!(imported.cds[1].contig_id, 2);
        assert_eq!(
            imported.cds[1].nucleotide_sequence.as_deref(),
            Some("ATGAAATAA")
        );
        assert!(imported
            .cds
            .iter()
            .all(|a| a.qualifiers.is_empty() && a.protein_sequence.is_none()));
        assert!(imported.contigs.iter().all(|c| c.sequence.is_none()));
        assert_eq!(
            imported
                .sequences
                .load(1, &imported.contigs[1])
                .unwrap()
                .sequence
                .unwrap(),
            b"ATGAAATAAGCC"
        );
        let conn = imported.writer.lock_conn().unwrap();
        assert_eq!(
            conn.query_row(
                "SELECT Protein_sequence FROM Annotation_sequence WHERE Annotation_id = 2",
                [],
                |r| r.get::<_, String>(0)
            )
            .unwrap(),
            "MK*"
        );
        let without_cds = import_annotations(
            &inputs,
            &directory.path().join("without_cds.db"),
            &config,
            false,
            false,
            |_, _, _| {},
        )
        .unwrap();
        assert!(without_cds.cds.is_empty());
    }

    #[test]
    fn mag_pipeline_reads_spooled_sequences_for_gc_and_cleans_up() {
        let directory = tempfile::tempdir().unwrap();
        let inputs = directory.path().join("inputs");
        std::fs::create_dir(&inputs).unwrap();
        let mut config = config(ViewMode::Mag);
        for name in ["a", "b"] {
            config.mag_manifest.push(crate::processing::MagInput {
                name: name.into(),
                gb_path: Some(fixture(&inputs, name)),
                asm_path: None,
            });
        }
        let output = directory.path().join("mag.db");
        crate::processing::run_all_samples(
            &inputs,
            Path::new(""),
            &[],
            &output,
            &[],
            &config,
            false,
            Path::new(""),
            "test",
        )
        .unwrap();
        let conn = duckdb::Connection::open(&output).unwrap();
        assert_eq!(
            conn.query_row("SELECT count(*) FROM MAG", [], |r| r.get::<_, i64>(0))
                .unwrap(),
            2
        );
        assert_eq!(
            conn.query_row(
                "SELECT count(*) FROM Contig WHERE GC_mean IS NOT NULL",
                [],
                |r| r.get::<_, i64>(0)
            )
            .unwrap(),
            2
        );
        assert_eq!(
            conn.query_row(
                "SELECT count(*) FROM Annotation_sequence WHERE Nucleotide_sequence = 'ATGAAATAA'",
                [],
                |r| r.get::<_, i64>(0)
            )
            .unwrap(),
            2
        );
        assert!(!std::fs::read_dir(directory.path()).unwrap().any(|e| e
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".thebigbam-sequences-")));
    }
}

#[cfg(test)]
mod database_sequence_tests {
    use super::*;
    #[test]
    fn candidate_dna_is_loaded_by_id_without_retaining_assembly_sequences() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("sequences.db");
        let conn = duckdb::Connection::open(&path).unwrap();
        conn.execute_batch("CREATE TABLE Contig (Contig_id BIGINT, Contig_length BIGINT); INSERT INTO Contig VALUES (7,6),(9,0); CREATE TABLE Contig_sequence (Contig_id BIGINT, Sequence VARCHAR);
            INSERT INTO Contig_sequence VALUES (7,'ATGCAA'),(9,NULL)").unwrap();
        let store = SequenceStore::from_connection(directory.path(),conn.try_clone().unwrap(),vec![7,9,10]).unwrap();
        assert!(store.entries.is_empty());
        assert_eq!(store.reference_length(0).unwrap(), 6);
        assert_eq!(store.load_sequence(0).unwrap(),Some(b"ATGCAA".to_vec()));
        assert_eq!(store.load_sequence(1).unwrap(),None);
        assert_eq!(store.load_sequence(2).unwrap(),None);
    }
    #[test]
    fn shared_reader_sees_live_writer_catalog_without_reopening_a_file() {
        let directory = tempfile::tempdir().unwrap();
        let writer = duckdb::Connection::open_in_memory().unwrap();
        writer.execute_batch("CREATE TABLE Contig (Contig_id BIGINT, Contig_length BIGINT);
            CREATE TABLE Contig_sequence (Contig_id BIGINT, Sequence VARCHAR);
            INSERT INTO Contig VALUES (7,6); INSERT INTO Contig_sequence VALUES (7,'ATGCAA')").unwrap();
        let store = SequenceStore::from_connection(directory.path(), writer.try_clone().unwrap(), vec![7]).unwrap();
        assert_eq!(store.reference_length(0).unwrap(), 6);
        assert_eq!(store.load_sequence(0).unwrap(), Some(b"ATGCAA".to_vec()));
        writer.execute("UPDATE Contig SET Contig_length=9 WHERE Contig_id=7", []).unwrap();
        assert_eq!(store.reference_length(0).unwrap(), 9);
    }

}
