//! Filtering and persistence shared by collected sample execution paths.
//! Schedulers retain progress reporting, failure counts, and wall-clock ownership.
use super::{get_rss_mb, mag_coverage_fails, SampleTimings};
use crate::db::{DbWriter, MicrodiversityData, MisassemblyData, SideMisassemblyData, TopologyData};
use crate::types::{PackagingData, PresenceData, SequencingType};
use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::fs;

pub(super) struct SampleData {
    pub presences: Vec<PresenceData>,
    pub packaging: Vec<PackagingData>,
    pub misassembly: Vec<MisassemblyData>,
    pub microdiversity: Vec<MicrodiversityData>,
    pub side_misassembly: Vec<SideMisassemblyData>,
    pub topology: Vec<TopologyData>,
    pub feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    pub mag_results: MagResultVec,
}
impl SampleData {
    pub fn filter(
        self,
        is_mag_mode: bool,
        map: &[(String, Vec<String>)],
        min_af: f64,
        min_cov: f64,
    ) -> (Self, usize) {
        if !is_mag_mode || map.is_empty() {
            return (self, 0);
        }
        let (drop, count) =
            filter_failing_mags(&self.presences, map, &self.mag_results, min_af, min_cov);
        if drop.is_empty() {
            return (self, count);
        }
        let (
            presences,
            packaging,
            misassembly,
            microdiversity,
            side_misassembly,
            topology,
            feature_blobs,
            mag_results,
        ) = apply_mag_filter(
            &drop,
            map,
            self.presences,
            self.packaging,
            self.misassembly,
            self.microdiversity,
            self.side_misassembly,
            self.topology,
            self.feature_blobs,
            self.mag_results,
        );
        (
            Self {
                presences,
                packaging,
                misassembly,
                microdiversity,
                side_misassembly,
                topology,
                feature_blobs,
                mag_results,
            },
            count,
        )
    }
}

/// Important sample notice: absence of retained results does not mean reads never mapped.
pub(super) fn report_empty_sample(name: &str, indexed_reads: u64, mag_mode: bool, min_af: f64, min_cov: f64) {
    let unit = if mag_mode { "MAG" } else { "contig" };
    eprintln!("INFO: Sample '{}': {} mapped alignments in BAM index; no {} coverage results retained (min_aligned_fraction={:.1}%, min_coverage_depth={:.1}). Sample recorded in DB. Check coverage thresholds and reference names/MAG membership.", name, indexed_reads, unit, min_af, min_cov);
}

/// Borrowed sample metadata, separate from potentially large result vectors.
pub(super) struct SampleIdentity<'a> {
    pub name: &'a str,
    pub sequencing_type: SequencingType,
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub bai_mapped_reads: u64,
    pub circular: bool,
}
impl SampleIdentity<'_> {
    pub fn insert(&self, writer: &DbWriter, empty: bool, min_af: f64, min_cov: f64) -> Result<i64> {
        if empty {
            report_empty_sample(self.name, self.bai_mapped_reads, writer.is_mag_mode(), min_af, min_cov);
        }
        let (total,mapped)=(self.total_reads,self.mapped_reads);
        writer.insert_sample(
            self.name,
            self.sequencing_type.as_str(),
            total,
            mapped,
            self.circular,
        )
    }
}

#[derive(Default)]
pub(super) struct SampleWriteDetails {
    pub insert_ns: u64,
    pub contig_ns: u64,
    pub mag_ns: u64,
    pub detail: crate::db::SampleWriteTimings,
}
impl SampleWriteDetails {
    pub fn apply(self, timings: &mut SampleTimings, elapsed: std::time::Duration) {
        timings.db_write = elapsed;
        timings.insert_sample_ns = self.insert_ns;
        timings.write_contig_data_ns = self.contig_ns;
        timings.mag_write_ns = self.mag_ns;
        timings.db_detail = self.detail;
    }
}

pub(super) fn write_sample_payload(
    writer: &DbWriter,
    identity: &SampleIdentity<'_>,
    sample_id: i64,
    data: &SampleData,
    diagnose: bool,
    details: &mut SampleWriteDetails,
) -> Result<()> {
    let start = std::time::Instant::now();
    details.detail = writer.write_sample_data_timed(
        identity.name,
        &data.presences,
        &data.packaging,
        &data.misassembly,
        &data.microdiversity,
        &data.side_misassembly,
        &data.topology,
        &data.feature_blobs,
        identity.circular,
        diagnose,
    )?;
    details.contig_ns = start.elapsed().as_nanos() as u64;
    details.mag_ns = write_mag_results(writer, identity.name, sample_id, &data.mag_results)?;
    Ok(())
}

/// Checkpoint failures propagate so incomplete samples remain retryable.
pub(super) fn checkpoint_sample(
    writer: &DbWriter,
    _name: &str,
    diagnose: bool,
    log: &mut Option<fs::File>,
) -> Result<()> {
    let rss_before = if diagnose { get_rss_mb() } else { 0.0 };
    let start = std::time::Instant::now();
    writer.checkpoint()?;
    if diagnose {
        let rss_after = get_rss_mb();
        let seconds = start.elapsed().as_secs_f64();
        if let Some(f) = log {
            use std::io::Write;
            let _ = writeln!(
                f,
                "  checkpoint           : {:>10.3} s  [RSS: {:.0} → {:.0} MB]",
                seconds, rss_before, rss_after
            );
        }
    }
    Ok(())
}

pub(super) fn filter_failing_mags(
    presences: &[PresenceData],
    mag_contig_map: &[(String, Vec<String>)],
    mag_results: &MagResultVec,
    min_aligned_fraction: f64,
    min_coverage_depth: f64,
) -> (HashSet<String>, usize) {
    let s_af = crate::types::get_column_scale("Coverage", "Aligned_fraction_percentage");
    let s_m = crate::types::get_column_scale("Coverage", "Coverage_mean");

    let stats_by_name: HashMap<&str, &crate::mag_blob::MagCoverageStats> = mag_results
        .iter()
        .map(|(name, _, stats)| (name.as_str(), stats))
        .collect();

    let presence_names: HashSet<&str> = presences.iter().map(|p| p.contig_name.as_str()).collect();

    let mut drop_set: HashSet<String> = HashSet::new();
    let mut n_dropped: usize = 0;

    for (mag_name, members) in mag_contig_map {
        let (af_mag, cov_mag) = if let Some(stats) = stats_by_name.get(mag_name.as_str()) {
            (
                stats.aligned_fraction_pct_x10 as f64 / s_af,
                stats.coverage_mean_x10 as f64 / s_m,
            )
        } else {
            let has_any_presence = members.iter().any(|c| presence_names.contains(c.as_str()));
            if has_any_presence {
                continue;
            } else {
                (0.0, 0.0)
            }
        };

        if mag_coverage_fails(af_mag, cov_mag, min_aligned_fraction, min_coverage_depth) {
            n_dropped += 1;
            for cname in members {
                drop_set.insert(cname.clone());
            }
        }
    }

    (drop_set, n_dropped)
}

pub(super) type MagResultVec = Vec<(String, i64, crate::mag_blob::MagCoverageStats)>;

/// Apply a contig drop set to all per-sample data vectors, removing entries
/// whose contig_name is in `drop`. Also filters `mag_results` by resolving
/// which MAGs had their contigs dropped via `mag_contig_map`.
pub(super) fn apply_mag_filter(
    drop: &HashSet<String>,
    mag_contig_map: &[(String, Vec<String>)],
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    mag_results: MagResultVec,
) -> (
    Vec<PresenceData>,
    Vec<PackagingData>,
    Vec<MisassemblyData>,
    Vec<MicrodiversityData>,
    Vec<SideMisassemblyData>,
    Vec<TopologyData>,
    Vec<(String, i64, crate::blob::EncodedBlob)>,
    MagResultVec,
) {
    let drop_ids: HashSet<i64> = presences
        .iter()
        .filter(|p| drop.contains(&p.contig_name))
        .map(|p| p.contig_id)
        .collect();
    let dropped_mag_names: HashSet<&str> = mag_contig_map
        .iter()
        .filter(|(_, members)| members.iter().any(|c| drop.contains(c)))
        .map(|(name, _)| name.as_str())
        .collect();
    (
        presences
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        packaging
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        misassembly
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        microdiversity
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        side_misassembly
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        topology
            .into_iter()
            .filter(|p| !drop.contains(&p.contig_name))
            .collect(),
        feature_blobs
            .into_iter()
            .filter(|(_, cid, _)| !drop_ids.contains(cid))
            .collect(),
        mag_results
            .into_iter()
            .filter(|(name, _, _)| !dropped_mag_names.contains(name.as_str()))
            .collect(),
    )
}

fn write_mag_results(writer:&DbWriter,_name:&str,id:i64,results:&MagResultVec)->Result<u64> {
    if results.is_empty() {return Ok(0);}
    let start=std::time::Instant::now();
    let conn=writer.lock_conn()?;
    writer.write_mag_coverages(&conn,id,results)?;
    Ok(start.elapsed().as_nanos() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::features::FeatureArrays;
    use crate::types::{get_value_scale, ContigInfo};

    fn data() -> SampleData {
        let mut arrays = FeatureArrays::new(100);
        arrays.primary_reads.fill(1);
        SampleData {
            presences: vec![super::super::contig::build_presence(
                &arrays, 1, "a", 100.0, 2, 1.0,
            )],
            packaging: vec![],
            misassembly: vec![],
            microdiversity: vec![],
            side_misassembly: vec![],
            topology: vec![],
            feature_blobs: vec![(
                "primary_reads".into(),
                1,
                crate::blob::encode_dense_blob(
                    &vec![1; 100],
                    get_value_scale("primary_reads"),
                    100,
                ),
            )],
            mag_results: vec![(
                "mag".into(),
                1,
                crate::mag_blob::compute_mag_coverage_stats_from_raw(&vec![1; 100], 100, &[100], 2),
            )],
        }
    }

    #[test]
    fn shared_filter_retains_passes_and_removes_failing_mag_payload() {
        let map = vec![("mag".to_owned(), vec!["a".to_owned()])];
        let (pass, dropped) = data().filter(true, &map, 100.0, 1.0);
        assert_eq!(dropped, 0);
        assert_eq!(
            (
                pass.presences.len(),
                pass.feature_blobs.len(),
                pass.mag_results.len()
            ),
            (1, 1, 1)
        );
        let (fail, dropped) = data().filter(true, &map, 100.0, 1.1);
        assert_eq!(dropped, 1);
        assert!(
            fail.presences.is_empty()
                && fail.feature_blobs.is_empty()
                && fail.mag_results.is_empty()
        );
        let (contig, dropped) = data().filter(false, &map, 100.0, 1.1);
        assert_eq!((dropped, contig.presences.len()), (0, 1));
    }

    #[test]
    fn shared_writer_preserves_empty_counts_payload_and_diagnostics() {
        let directory = tempfile::tempdir().unwrap();
        let writer = DbWriter::create(
            &directory.path().join("writes.db"),
            &[ContigInfo {
                name: "a".into(),
                aliases: vec![],
                length: 100,
                sequence: None,
            }],
            &[],
            &[],
            true,
            false,
        )
        .unwrap();
        let normal = SampleIdentity {
            name: "normal",
            sequencing_type: SequencingType::ShortSingle,
            total_reads: 12,
            mapped_reads: 2,
            bai_mapped_reads: 20,
            circular: false,
        };
        let id = normal.insert(&writer, false, 50.0, 0.0).unwrap();
        let mut payload = data();
        payload.mag_results.clear();
        let mut details = SampleWriteDetails {
            insert_ns: 7,
            ..Default::default()
        };
        write_sample_payload(&writer, &normal, id, &payload, true, &mut details).unwrap();
        let mut timing = SampleTimings::default();
        details.apply(&mut timing, std::time::Duration::from_secs(1));
        assert_eq!(timing.insert_sample_ns, 7);
        assert_eq!(timing.db_write, std::time::Duration::from_secs(1));
        assert_eq!(
            (timing.db_detail.zoom_rows, timing.db_detail.chunk_rows),
            (1, 1)
        );
        let empty = SampleIdentity {
            name: "empty",
            ..normal
        };
        empty.insert(&writer, true, 50.0, 0.0).unwrap();
        let conn = writer.lock_conn().unwrap();
        let mut stmt = conn
            .prepare("SELECT Sample_name, Number_of_reads, Number_of_mapped_reads FROM Sample ORDER BY Sample_id")
            .unwrap();
        let rows = stmt
            .query_map([], |r| {
                Ok((
                    r.get::<_, String>(0)?,
                    r.get::<_, u64>(1)?,
                    r.get::<_, u64>(2)?,
                ))
            })
            .unwrap()
            .collect::<duckdb::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            rows,
            vec![("normal".into(), 12, 2), ("empty".into(), 12, 2)]
        );
        assert_eq!(
            conn.query_row("SELECT count(*) FROM Feature_blob", [], |r| r
                .get::<_, u64>(0))
                .unwrap(),
            1
        );
    }
}
