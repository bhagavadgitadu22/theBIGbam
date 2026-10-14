//! Global reference all-versus-all BLAST for reference collections.
//! DNA is exported in bounded database batches; query shards and hits stay on disk.
use super::BLAST_ARGS;
use crate::db::{DbWriter, RepeatsData};
use crate::types::ContigInfo;
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

fn numeric_id(label: &str) -> Result<i64> {
    label
        .strip_prefix("lcl|")
        .unwrap_or(label)
        .strip_prefix('c')
        .context("Invalid reference BLAST ID")?
        .parse()
        .map_err(Into::into)
}

pub(super) fn run_global(writer: &DbWriter, contigs: &[ContigInfo]) -> Result<(Vec<i64>, usize)> {
    let temporary = tempfile::Builder::new()
        .prefix("thebigbam-global-blast-")
        .tempdir()?;
    let subject = temporary.path().join("references.fa");
    let mut database_fasta = BufWriter::new(std::fs::File::create(&subject)?);
    let mut query_files = Vec::new();
    let mut query: Option<BufWriter<std::fs::File>> = None;
    let (mut query_bases, mut query_count) = (0usize, 0usize);
    let mut ids = Vec::new();
    let mut after = 0i64;
    eprintln!("Reference BLAST all: exporting the complete reference collection (including existing references)");
    {
        let conn = writer.lock_conn()?;
        loop {
            let batch = {
                let mut stmt = conn.prepare(
                    "SELECT Contig_id FROM Contig WHERE Contig_id>? ORDER BY Contig_id LIMIT 128",
                )?;
                let rows = stmt.query_map(duckdb::params![after], |r| r.get::<_, i64>(0))?;
                let ids = rows.collect::<duckdb::Result<Vec<_>>>()?;
                ids
            };
            if batch.is_empty() {
                break;
            }
            let first = batch[0];
            after = *batch.last().unwrap();
            let mut stmt=conn.prepare("SELECT Contig_id,Sequence FROM Contig_sequence WHERE Contig_id BETWEEN ? AND ? AND Sequence IS NOT NULL ORDER BY Contig_id")?;
            for row in stmt.query_map(duckdb::params![first, after], |r| {
                Ok((r.get::<_, i64>(0)?, r.get::<_, String>(1)?))
            })? {
                let (id, sequence) = row?;
                if sequence.is_empty() {
                    continue;
                }
                if query.is_none() || query_bases >= 32 * 1024 * 1024 || query_count >= 512 {
                    if let Some(mut previous) = query.take() {
                        previous.flush()?;
                    }
                    let path = temporary
                        .path()
                        .join(format!("queries-{}.fa", query_files.len()));
                    query = Some(BufWriter::new(std::fs::File::create(&path)?));
                    query_files.push(path);
                    query_bases = 0;
                    query_count = 0;
                }
                writeln!(database_fasta, ">c{id}")?;
                writeln!(query.as_mut().unwrap(), ">c{id}")?;
                for chunk in sequence.as_bytes().chunks(80) {
                    database_fasta.write_all(chunk)?;
                    database_fasta.write_all(b"\n")?;
                    query.as_mut().unwrap().write_all(chunk)?;
                    query.as_mut().unwrap().write_all(b"\n")?;
                }
                ids.push(id);
                query_bases += sequence.len();
                query_count += 1;
            }
        }
    }
    database_fasta.flush()?;
    drop(database_fasta);
    if let Some(mut query) = query.take() {
        query.flush()?;
    }
    anyhow::ensure!(!ids.is_empty(), "--blast all requires reference DNA in Contig_sequence; provide annotations with sequences or --assembly");
    let output = std::process::Command::new("makeblastdb")
        .arg("-in")
        .arg(&subject)
        .args(["-dbtype", "nucl", "-parse_seqids"])
        .output()
        .context("Building global reference BLAST database")?;
    anyhow::ensure!(
        output.status.success(),
        "Global makeblastdb failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    eprintln!("Reference BLAST all: {} references, {} bounded query shards; comparisons include different contigs",ids.len(),query_files.len());
    let done = AtomicUsize::new(0);
    let outputs = query_files
        .par_iter()
        .enumerate()
        .map(|(index, query)| -> Result<_> {
            let path = temporary.path().join(format!("hits-{index}.tsv"));
            let output = std::process::Command::new("blastn")
                .arg("-query")
                .arg(query)
                .arg("-db")
                .arg(&subject)
                .args(BLAST_ARGS)
                .args(["-outfmt", "6"])
                .arg("-max_target_seqs")
                .arg(ids.len().to_string())
                .arg("-out")
                .arg(&path)
                .output()
                .context("Running global reference BLAST")?;
            anyhow::ensure!(
                output.status.success(),
                "Global BLAST shard {index} failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let completed = done.fetch_add(1, Ordering::Relaxed) + 1;
            if completed % 10 == 0 || completed == query_files.len() {
                eprintln!(
                    "Reference BLAST all: {completed}/{} query shards finished",
                    query_files.len()
                );
            }
            Ok(path)
        })
        .collect::<Result<Vec<_>>>()?;
    // Do not replace existing tracks until every subprocess has succeeded.
    writer.reset_reference_blast()?;
    let names = contigs
        .iter()
        .map(|contig| {
            writer
                .get_contig_id(&contig.name)
                .map(|id| (id, contig.name.as_str()))
                .context("Missing reference contig ID")
        })
        .collect::<Result<HashMap<_, _>>>()?;
    let mut repeats = Vec::new();
    let mut hits = Vec::new();
    let mut hit_count = 0;
    for path in outputs {
        for line in BufReader::new(std::fs::File::open(path)?).lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            anyhow::ensure!(fields.len() >= 12, "Malformed global BLAST row");
            let q = numeric_id(fields[0])?;
            let s = numeric_id(fields[1])?;
            let identity = fields[2].parse::<f64>()?;
            let p1 = fields[6].parse::<i32>()?;
            let p2 = fields[7].parse::<i32>()?;
            let p1p = fields[8].parse::<i32>()?;
            let p2p = fields[9].parse::<i32>()?;
            let query = names.get(&q).context("Unknown query reference ID")?;
            let subject = names.get(&s).context("Unknown subject reference ID")?;
            if q == s {
                if p1 == p1p && p2 == p2p {
                    continue;
                }
                repeats.push(RepeatsData {
                    contig_name: (*query).to_owned(),
                    position1: p1,
                    position2: p2,
                    position1prime: p1p,
                    position2prime: p2p,
                    pident: identity,
                    is_direct: (p1 < p2) == (p1p < p2p),
                });
            } else {
                hits.push((
                    (*query).to_owned(),
                    (*subject).to_owned(),
                    p1,
                    p2,
                    p1p,
                    p2p,
                    identity,
                ));
                hit_count += 1;
            }
            if repeats.len() + hits.len() >= 10_000 {
                writer.write_mag_blast_batch(&repeats, &hits)?;
                repeats.clear();
                hits.clear();
            }
        }
    }
    writer.write_mag_blast_batch(&repeats, &hits)?;
    writer.update_duplication_from_database()?;
    writer.write_blast_parameters()?;
    writer.lock_conn()?.execute_batch("INSERT OR REPLACE INTO Database_metadata VALUES ('Blast_scope','all_reference_collection'),('Blast_state','complete')")?;
    Ok((ids, hit_count))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numeric_labels_accept_blast_local_ids_and_reject_unknown_names() {
        assert_eq!(numeric_id("c7").unwrap(), 7);
        assert_eq!(numeric_id("lcl|c7").unwrap(), 7);
        assert!(numeric_id("other").is_err());
    }
    #[test]
    #[ignore = "requires BLAST+ on PATH"]
    fn global_search_stores_cross_mag_hits_without_polluting_within_mag_tracks() {
        let directory = tempfile::tempdir().unwrap();
        let dna=b"TGAGGGTGAGCTTTTGTGCACCTCGCCTATAGACTTAAGCTAACGCCTGAAGCTCATTTCACGCTGATTAGGGGATTGTAATCAATGCGAAGCGTAACCAATGGCAGCGCTGGTCACGCGGATGGTGCCAATTAACTGTGTGTGTTATACCAACTGCGATAGCACTGCTGCTACTCCCGCGGCTTCATCATAGCGGTCATGCAGAAACAACTATCGATTCGATCGACGGGTAGCATTCGTACGGTGCTCCCCCGGGTCGAGGGCGCGCAACAGCGGACTTCAAGTTGCTATCCGCCTACGCGG";
        let contigs = ["a", "b"]
            .iter()
            .map(|name| ContigInfo {
                name: name.to_string(),
                aliases: vec![],
                length: dna.len(),
                sequence: Some(dna.to_vec()),
            })
            .collect::<Vec<_>>();
        let writer = DbWriter::create(
            &directory.path().join("global.db"),
            &contigs,
            &[],
            &[],
            true,
            true,
        )
        .unwrap();
        writer
            .write_mags(&[
                ("one".into(), vec!["a".into()]),
                ("two".into(), vec!["b".into()]),
            ])
            .unwrap();
        let (ids, hits) = run_global(&writer, &contigs).unwrap();
        assert_eq!(ids, vec![1, 2]);
        assert!(hits >= 2);
        writer.write_intercontig_hit_features().unwrap();
        let conn = writer.lock_conn().unwrap();
        assert!(
            conn.query_row(
                "SELECT count(*) FROM Contig_blast_hits WHERE Contig_id_1=1 AND Contig_id_2=2",
                [],
                |r| r.get::<_, u64>(0)
            )
            .unwrap()
                > 0
        );
        assert_eq!(
            conn.query_row(
                "SELECT count(*) FROM Contig_blob WHERE Feature_id IN (?,?)",
                duckdb::params![
                    crate::types::feature_name_to_id("hit_count_within_mag").unwrap(),
                    crate::types::feature_name_to_id("hit_identity_within_mag").unwrap()
                ],
                |r| r.get::<_, u64>(0)
            )
            .unwrap(),
            0
        );
        assert_eq!(conn.query_row(
            "SELECT count(*) FROM Contig_blob WHERE Feature_id IN (?,?)",
            duckdb::params![crate::types::feature_name_to_id("hit_count").unwrap(), crate::types::feature_name_to_id("hit_identity").unwrap()],
            |r| r.get::<_, u64>(0)).unwrap(), 4);

    }
}
