//! Parallel BAM file processing and feature calculation.
//!
//! Orchestrates parallel processing of BAM files:
//! 1. Each BAM file processed by separate thread (rayon)
//! 2. Results collected in memory
//! 3. Written to DuckDB sequentially (DuckDB is single-writer)
//!
//! Per-contig feature work is parallelized. Large MAG calculations stream
//! bounded batches to a dedicated writer while computation continues.

use anyhow::{Context, Result};
use atty::Stream;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::htslib;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::sync::mpsc::sync_channel;
use std::thread;

use crate::annotation_import::SequenceStore;
use crate::bam_reader::{detect_sequencing_type, process_contig_streaming_counted, BamPhaseTimings};
#[cfg(test)]
use crate::bam_reader::process_contig_streaming;
use crate::compress::{
    Run,
    add_compressed_feature_with_stats,
};
use crate::db::{DbWriter, FinalizeTimings, MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData, GCContentData, RepeatsData};
use crate::gc_content::{compute_gc_content_and_skew, GCParams};
use crate::features::{CdsIndex, FeatureArrays, ModuleFlags, compute_codon_changes_from_summaries};
use crate::parser::{parse_annotations, compute_annotation_sequences};
use crate::types::{
    ContigInfo, FeatureAnnotation, PackagingData, PresenceData, SequencingType
};

use crate::processing_phage_packaging::{
    classify_packaging_areas, filter_and_merge_to_areas_with_diagnostics,
    is_valid_terminal_repeat, DtrRegion, PeakArea, PhageTerminiConfig,
};
use crate::processing_completeness::compute_all_metrics;

pub mod annotation_data;
use annotation_data::{CdsSource, DatabaseCdsSource, needs_codon_annotations};
mod reference_blast;
mod contig;
mod writing;
use writing::{SampleData, SampleIdentity, SampleWriteDetails, write_sample_payload, checkpoint_sample, MagResultVec};
#[cfg(test)]
use writing::{filter_failing_mags, apply_mag_filter};
use contig::build_presence;

/// Aggregation level for the processing pipeline.
///
/// - `Contig` (default): each contig is processed and filtered independently.
/// - `Mag`: contigs are grouped by MAG (one MAG per input file); filters and
///   some metrics are evaluated at MAG level. See `MagInput` + `mag_manifest`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ViewMode {
    Contig,
    Mag,
}

impl ViewMode {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "mag" => ViewMode::Mag,
            _ => ViewMode::Contig,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            ViewMode::Contig => "contig",
            ViewMode::Mag => "mag",
        }
    }
}

/// One entry in the MAG manifest: a MAG name plus the per-MAG files that define
/// its contigs. Empty `gb_path` / `asm_path` means that side is absent for this
/// MAG (e.g. a directory was supplied on only one of `-g` / `-a`).
#[derive(Clone, Debug)]
pub struct MagInput {
    pub name: String,
    pub gb_path: Option<PathBuf>,
    pub asm_path: Option<PathBuf>,
}

/// A group of contigs belonging to one MAG, ready for batch processing.
struct MagJobGroup {
    mag_name: String,
    mag_id: i64,
    mag_length: u32,
    members: Vec<MagContigJob>,
}

/// One contig within a MAG job group.
struct MagContigJob {
    tid: u32,
    contig_id: i64,
    ref_name: String,
    ref_length: usize,
    offset_in_mag: u32,
}

/// Result of processing one MAG group (all member contigs + MAG-level stats).
struct MagGroupResult {
    mag_name: String,
    mag_id: i64,
    contig_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    mag_cov_stats: crate::mag_blob::MagCoverageStats,
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    mapped_reads: u64,
}

/// Public callers may supply positional metadata; pipeline maps carry actual DB IDs.
pub trait ContigLookup: Sync { fn find(&self, name: &str) -> Option<(usize,i64)>; }
impl ContigLookup for HashMap<&str,usize> {
    fn find(&self,name:&str) -> Option<(usize,i64)> { self.get(name).map(|&i|(i,i as i64+1)) }
}
impl ContigLookup for HashMap<&str,(usize,i64)> {
    fn find(&self,name:&str) -> Option<(usize,i64)> { self.get(name).copied() }
}
fn database_contig_lookup<'a>(contigs:&'a [ContigInfo], writer:&DbWriter) -> Result<HashMap<&'a str,(usize,i64)>> {
    let mut map=HashMap::new();
    for (index,contig) in contigs.iter().enumerate() {
        let id=writer.get_contig_id(&contig.name).with_context(||format!("Missing DB contig ID for '{}'",contig.name))?;
        map.insert(contig.name.as_str(),(index,id));
        for alias in &contig.aliases { map.entry(alias.as_str()).or_insert((index,id)); }
    }
    Ok(map)
}

/// Contig→MAG membership mapping, built once per run from mag_contig_map + DB IDs.
/// Maps contig_name → (mag_name, mag_id, offset_in_mag, mag_length).
pub type MagMemberLookup = HashMap<String, (String, i64, u32, u32)>;

/// Shared nucleotide search parameters for reference and candidate scopes.
pub(crate) const BLAST_ARGS: &[&str] = &["-task", "megablast", "-word_size", "28", "-dust", "20 64 1", "-strand", "both", "-evalue", "1e-10", "-num_threads", "1"];
pub(crate) const BLAST_PARAMETERS: &str = "v2;task=megablast;word_size=28;dust=20 64 1;strand=both;evalue=1e-10;threads=1";
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlastMode { None, Termini, Mags, All }
impl BlastMode {
    pub fn parse(value: &str) -> Result<Self> { match value { "none" => Ok(Self::None), "termini" => Ok(Self::Termini), "mags" => Ok(Self::Mags), "all" => Ok(Self::All), _ => anyhow::bail!("blast must be none, termini, mags, or all") } }
    pub fn reference_search(self) -> bool { matches!(self,Self::All|Self::Mags) }
    pub fn as_str(self) -> &'static str { match self { Self::None => "none", Self::Termini => "termini", Self::Mags => "mags", Self::All => "all" } }
}

/// Configuration for processing.
#[derive(Clone)]
pub struct ProcessConfig {
    pub threads: usize,
    pub min_aligned_fraction: f64,
    pub min_coverage_depth: f64,
    pub bar_ratio: f64,
    /// Sequencing type: if None, auto-detect per sample; if Some, use for all samples
    pub sequencing_type: Option<SequencingType>,
    /// Phage termini detection configuration
    pub phagetermini_config: PhageTerminiConfig,
    /// GC content and GC skew parameters
    pub gc_params: GCParams,
    /// Adaptive smoothing percentage for dense features (0.0 = disabled).
    /// Consecutive positions within this % of the run value are collapsed.
    pub variation_percentage: f64,
    /// Minimum absolute event count for a sparse feature position to be kept.
    /// Applied in addition to bar_ratio filtering: position kept only if
    /// `value >= coverage × bar_ratio AND value >= min_occurrences`.
    /// Default: 2.
    pub min_occurrences: u32,
    /// When true, record per-sample phase timings and write a timing log
    /// alongside the output database. Hidden CLI flag (--diagnose).
    pub enable_timing: bool,
    /// When true, run BLAST-based repeat detection (autoblast in contig mode,
    /// per-MAG blast in MAG mode). Off by default.
    pub blast: BlastMode,
    pub sample_name: String,
    pub primary_reads_seen: Arc<AtomicU64>,
    /// Shared on-demand terminal evidence, initialized by run_all_samples.
    pub terminal_repeats: Option<Arc<crate::terminal_repeats::TerminalRepeatSearch>>,
    /// Aggregation level. `Contig` is the historical default.
    pub view_mode: ViewMode,
    /// Per-MAG input manifest. Empty in contig mode; one entry per MAG in MAG
    /// mode. The Rust side treats this as the authoritative source of MAG
    /// membership; `genbank_path` / `assembly_path` remain the original user
    /// inputs (which may be directory strings in MAG mode).
    pub mag_manifest: Vec<MagInput>,
    /// Counter for (Contig, Sample) pairs dropped by the per-contig AF or coverage
    /// thresholds in contig mode. Incremented inside `process_sample` workers;
    /// a summary is printed at the end of `run_all_samples`.
    #[doc(hidden)]
    pub contig_drops_counter: Arc<AtomicUsize>,
}

impl ProcessConfig {
    /// Parse sequencing type string into SequencingType enum.
    /// If empty or invalid, falls back to auto-detection from BAM file.
    pub fn parse_sequencing_type(seq_type_str: &str) -> Option<SequencingType> {
        match seq_type_str.to_lowercase().as_str() {
            "long" => Some(SequencingType::Long),
            "paired-short" => Some(SequencingType::ShortPaired),
            "single-short" => Some(SequencingType::ShortSingle),
            _ => None,
        }
    }
}

/// Result of processing all samples.
pub struct ProcessResult {
    pub samples_processed: usize,
    pub samples_failed: usize,
    /// Samples written to DB with zero contig coverage (all contigs failed filters).
    pub samples_no_coverage: usize,
    pub total_time_secs: f64,
    pub processing_time_secs: f64,
    pub writing_time_secs: f64,
}

/// Per-sample timing breakdown populated when ProcessConfig::enable_timing is set.
#[derive(Default, Clone)]
pub struct SampleTimings {
    pub sample_name: String,
    pub bam_path: String,
    pub n_contigs_in_bam: usize,
    pub n_contigs_processed: usize,
    pub total_reads: u64,
    // Single-threaded phases inside process_sample (wall-clock).
    pub header_read: std::time::Duration,
    pub total_read_count: std::time::Duration,
    pub filter_tids: std::time::Duration,
    pub par_iter_wall: std::time::Duration,
    pub merge_results: std::time::Duration,
    // Per-contig sub-phase sums inside par_iter (nanoseconds, summed across threads).
    // bam_open_ns now counts map_init calls (~once per rayon worker chunk),
    // not once per contig.
    pub bam_reader_workers: usize,
    pub bam_reader_estimated_mib: u64,
    pub bam_reader_start_rss_mib: u64,
    pub rss_before_reader_open_mib: f64,
    pub bam_open_ns: u64,
    pub streaming_ns: u64,
    pub feature_calc_ns: u64,
    pub dense_encode_ns: u64,
    pub sparse_encode_ns: u64,
    pub db_detail: crate::db::SampleWriteTimings,
    pub db_overlapped: bool,
    pub mag_queue_wait_ns: u64,
    pub terminal: crate::terminal_repeats::TerminalTimings,
    // Breakdown of streaming_ns (see BamPhaseTimings in bam_reader.rs).
    pub stream_fetch_seek_ns: u64,
    pub stream_records_ns: u64,
    pub stream_process_read_ns: u64,
    pub stream_finalize_ns: u64,
    // Filled by the sample-level caller, not process_sample itself.
    pub db_write: std::time::Duration,
    pub insert_sample_ns: u64,
    pub write_contig_data_ns: u64,
    pub mag_write_ns: u64,
    pub total: std::time::Duration,
}

/// Thread-safe accumulators for per-contig sub-phase timings inside the par_iter.
struct TimingAccum {
    bam_open_ns: AtomicU64,
    streaming_ns: AtomicU64,
    feature_calc_ns: AtomicU64,
    mag_queue_wait_ns: AtomicU64,
    sparse_encode_ns: AtomicU64,
    dense_encode_ns: AtomicU64,
    stream_fetch_seek_ns: AtomicU64,
    stream_records_ns: AtomicU64,
    stream_process_read_ns: AtomicU64,
    stream_finalize_ns: AtomicU64,
}

impl TimingAccum {
    fn new() -> Self {
        Self {
            bam_open_ns: AtomicU64::new(0),
            streaming_ns: AtomicU64::new(0),
            feature_calc_ns: AtomicU64::new(0),
            mag_queue_wait_ns: AtomicU64::new(0),
            sparse_encode_ns: AtomicU64::new(0),
            dense_encode_ns: AtomicU64::new(0),
            stream_fetch_seek_ns: AtomicU64::new(0),
            stream_records_ns: AtomicU64::new(0),
            stream_process_read_ns: AtomicU64::new(0),
            stream_finalize_ns: AtomicU64::new(0),
        }
    }
}

/// Pre-sample phase timings bundled for passing through to finish_timing_log.
#[derive(Default, Clone)]
struct PreSampleTimings {
    parse_secs: f64,
    parse_rss_mb: f64,
    db_create_secs: f64,
    db_create_rss_mb: f64,
    circ_secs: f64,
    circ_rss_mb: f64,
    autoblast_secs: f64,
    autoblast_rss_mb: f64,
    gc_secs: [f64; 4],
    gc_rss_mb: f64,
    repeat_blob_secs: f64,
    repeat_rss_mb: f64,
    write_mags_secs: f64,
    write_mags_rss_mb: f64,
    cds_build_secs: f64,
    cds_build_rss_mb: f64,
}

/// Check if a BAM file is missing MD tags by sampling first few reads.
/// Returns true if MD tags are missing, false otherwise.
fn check_missing_md_tags(bam: &mut IndexedReader) -> bool {
    let header = bam.header().clone();
    let target_names = header.target_names();
    let first_tid = target_names.first();

    if first_tid.is_none() {
        return false;
    }

    let contig_name = match std::str::from_utf8(first_tid.unwrap()) {
        Ok(s) => s,
        Err(_) => return false,
    };

    if bam.fetch(contig_name).is_err() {
        return false;
    }

    let mut checked = 0;
    let mut missing_md = 0;

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        checked += 1;
        if record.aux(b"MD").is_err() {
            missing_md += 1;
        }

        if checked >= 10 {
            break;
        }
    }

    checked > 0 && missing_md == checked
}

/// Detect whether a sample BAM was mapped circularly (SAM-spec circular BAM).
///
/// Detection strategy:
/// 1. Check BAM `@CO` header for `theBIGbam:circular=true/false` (written by mapping-per-sample)
/// 2. If no tag found, assume linear
fn detect_sample_circularity(
    bam: &IndexedReader,
    _sample_name: &str,
) -> Result<bool> {
    let header = bam.header();

    // 1. Check @CO header for theBIGbam:circular=true/false
    let header_text = String::from_utf8_lossy(header.as_bytes());
    for line in header_text.lines() {
        if let Some(rest) = line.strip_prefix("@CO\t") {
            if let Some(val) = rest.strip_prefix("theBIGbam:circular=") {
                match val.trim() {
                    "true" => {
                        return Ok(true);
                    }
                    "false" => {
                        return Ok(false);
                    }
                    _ => {} // unrecognized value, fall through to length check
                }
            }
        }
    }

    // 2. No header tag found — assume linear
    Ok(false)
}

/// Strict upfront validation of all inputs.
/// Collects all errors across all samples, prints them, and returns Err if any found.
/// This runs as the very first step before any parsing or database creation.
fn validate_inputs(
    bam_files: &[PathBuf],
    modules: &[String],
) -> Result<()> {
    let flags = ModuleFlags::from_modules(modules);
    let needs_md = flags.needs_md();
    let mut errors: Vec<String> = Vec::new();

    // 1. BAM index check — for each BAM file, check .bam.bai or .bai exists on disk
    for bam_path in bam_files {
        let sample_name = bam_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .replace("_with_MD", "");

        let bai_path1 = bam_path.with_extension("bam.bai");
        let bai_path2 = bam_path.with_extension("bai");
        if !bai_path1.exists() && !bai_path2.exists() {
            errors.push(format!(
                "ERROR: Sample '{}' is missing a BAM index (.bai).\n\
                 If the file is coordinate-sorted, create the index with:\n\
                 \x20 samtools index {}\n\
                 Otherwise, sort and index with:\n\
                 \x20 samtools sort {} -o {}_sorted.bam && samtools index {}_sorted.bam",
                sample_name,
                bam_path.display(),
                bam_path.display(),
                sample_name,
                sample_name,
            ));
        }
    }

    // 2. MD tag check (only when modules include Misalignment or Phage termini)
    if needs_md {
        let mut samples_missing_md: Vec<String> = Vec::new();
        for bam_path in bam_files {
            let sample_name = bam_path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .replace("_with_MD", "");

            // Only check BAMs that have an index (otherwise we already reported the error)
            let bai_path1 = bam_path.with_extension("bam.bai");
            let bai_path2 = bam_path.with_extension("bai");
            if !bai_path1.exists() && !bai_path2.exists() {
                continue;
            }

            let mut check_bam = match IndexedReader::from_path(bam_path) {
                Ok(b) => b,
                Err(_) => continue,
            };
            if check_missing_md_tags(&mut check_bam) {
                samples_missing_md.push(sample_name);
            }
        }
        if !samples_missing_md.is_empty() {
            errors.push(format!(
                "ERROR: Sample(s) missing MD tags (required by Misalignment, Phage termini modules): {}\n\
                 MD tags encode reference base information needed for mismatch and soft-clipping analysis.\n\
                 Add MD tags with:\n\
                 \x20 samtools calmd -b <input>.bam reference.fa > <output>_md.bam && samtools index <output>_md.bam\n\
                 Or remove these modules from -m/--modules.",
                samples_missing_md.join(", "),
            ));
        }
    }

    if !errors.is_empty() {
        // Include full error details in the anyhow error so they propagate
        // through PyO3 into the Python exception message (visible without RUST_BACKTRACE)
        let combined = errors.join("\n\n");
        return Err(anyhow::anyhow!(
            "Validation failed with {} error(s):\n\n{}",
            errors.len(),
            combined,
        ));
    }

    Ok(())
}

/// Detect circularity for all samples upfront. Returns a map of BAM path → is_circular.
fn detect_all_sample_circularities(
    bam_files: &[PathBuf],
) -> Result<HashMap<PathBuf, bool>> {
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Auto-detecting circularity per sample...");
    eprintln!();
    let mut result = HashMap::new();
    for bam_path in bam_files {
        let sample_name = bam_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .replace("_with_MD", "");

        let bam = IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM for circularity detection: {}", bam_path.display()))?;

        let is_circular = detect_sample_circularity(&bam, &sample_name)?;
        result.insert(bam_path.clone(), is_circular);
    }
    Ok(result)
}

/// Merge sequences from a FASTA file into existing contig info.
/// For each FASTA record, finds matching contig by name and sets contig.sequence.
/// Contigs present in the FASTA but not in the annotation are added automatically
/// so they still get coverage tracking (they just won't have annotations).
fn merge_sequences_from_fasta(contigs: &mut Vec<ContigInfo>, fasta_path: &Path) -> Result<()> {
    let records = crate::parser::parse_fasta(fasta_path)?;
    let name_to_idx: HashMap<String, usize> = contigs.iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), i))
        .collect();
    let mut matched = 0usize;
    let mut added = 0usize;
    for (name, seq) in &records {
        if let Some(&idx) = name_to_idx.get(name) {
            contigs[idx].sequence = Some(seq.clone());
            matched += 1;
        } else {
            contigs.push(ContigInfo {
                name: name.clone(),
                aliases: vec![],
                length: seq.len(),
                sequence: Some(seq.clone()),
            });
            added += 1;
        }
    }
    if added > 0 {
        contigs.sort_by(|a, b| a.name.cmp(&b.name));
        eprintln!("Added {} contigs from FASTA that had no annotations", added);
    }
    eprintln!("Merged sequences for {}/{} contigs from FASTA", matched + added, contigs.len());
    Ok(())
}

/// Write contigs with sequences to a temp FASTA and return the temp file handle.
fn write_contigs_fasta(contigs: &[&ContigInfo]) -> Result<tempfile::NamedTempFile> {
    let mut temp_file = tempfile::NamedTempFile::new()
        .context("Failed to create temp FASTA")?;
    use std::io::Write;
    for contig in contigs {
        if let Some(seq) = contig.sequence.as_ref() {
            writeln!(temp_file, ">{}", contig.name)?;
            for chunk in seq.chunks(80) {
                temp_file.write_all(chunk)?;
                writeln!(temp_file)?;
            }
        }
    }
    temp_file.flush()?;
    Ok(temp_file)
}

/// Per-contig BLAST self-alignment for intra-contig repeat detection.
/// Each contig is BLASTed against itself independently, parallelized with rayon.
fn run_autoblast(contigs: &[ContigInfo], sequences: Option<&SequenceStore>) -> Result<Vec<RepeatsData>> {
    let contigs_with_seq: Vec<usize> = contigs.iter().enumerate()
        .filter(|(i, c)| sequences.map_or(c.sequence.is_some(), |s| s.has_sequence(*i)))
        .map(|(i, _)| i).collect();
    if contigs_with_seq.is_empty() {
        return Ok(Vec::new());
    }

    match std::process::Command::new("blastn").arg("-version").output() {
        Ok(o) if o.status.success() => {}
        _ => {
            anyhow::bail!("blastn not found on PATH. Install BLAST+ (e.g. `conda install -c bioconda blast`).");
        }
    }

    let total_contigs = contigs_with_seq.len();
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Reference BLAST repeat detection on {} contigs...", total_contigs);

    let done = std::sync::atomic::AtomicUsize::new(0);

    let all_repeats: Vec<Vec<RepeatsData>> = contigs_with_seq
        .par_iter()
        .map(|&i| -> Result<Vec<RepeatsData>> {
            let loaded;
            let contig = if let Some(store) = sequences {
                loaded = store.load(i, &contigs[i])?; &loaded
            } else { &contigs[i] };
            let temp_file = write_contigs_fasta(&[contig])?;
            let temp_path = temp_file.path().to_path_buf();
            let output = std::process::Command::new("blastn")
                .arg("-query").arg(&temp_path).arg("-subject").arg(&temp_path)
                .args(BLAST_ARGS).arg("-outfmt").arg("6")
                .output().with_context(|| format!("Starting reference BLAST for '{}'",contig.name))?;
            anyhow::ensure!(output.status.success(), "Reference BLAST failed for '{}': {}",contig.name,String::from_utf8_lossy(&output.stderr));

            let stdout = String::from_utf8_lossy(&output.stdout);
            let mut repeats = Vec::new();
            for line in stdout.lines() {
                let line = line.trim();
                if line.is_empty() { continue; }
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 10 { continue; }
                let pident: f64 = f[2].parse().unwrap_or(0.0);
                let qstart: i32 = f[6].parse().unwrap_or(0);
                let qend: i32 = f[7].parse().unwrap_or(0);
                let sstart: i32 = f[8].parse().unwrap_or(0);
                let send: i32 = f[9].parse().unwrap_or(0);
                if qstart == sstart && qend == send { continue; }
                let is_direct = (qstart < qend && sstart < send) || (qstart > qend && sstart > send);
                repeats.push(RepeatsData {
                    contig_name: contig.name.clone(),
                    position1: qstart, position2: qend,
                    position1prime: sstart, position2prime: send,
                    pident, is_direct,
                });
            }

            let n = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if n % 500 == 0 || n == total_contigs {
                eprintln!("  Reference BLAST: {}/{} contigs", n, total_contigs);
            }
            Ok(repeats)
        })
        .collect::<Result<_>>()?;

    let total: Vec<RepeatsData> = all_repeats.into_iter().flatten().collect();
    eprintln!("Found {} repeat hits across {} contigs", total.len(), contigs_with_seq.len());
    Ok(total)
}

/// Per-MAG inter-contig BLAST. Writes one FASTA per MAG containing all member contigs,
/// runs blastn query=subject, and returns (repeats, inter_contig_hits).
///
/// - Intra-contig hits (qid == sid, not exact self-coordinates) → RepeatsData
/// - Inter-contig hits (qid != sid) → written to Contig_blast_hits
///
/// In MAG mode this REPLACES the per-contig self-BLAST: intra-contig hits emerge naturally
/// from each MAG's all-vs-all blast output.
fn run_mag_interblast(
    contigs: &[ContigInfo],
    mag_manifest: &[(String, Vec<String>)],
    db_writer: &DbWriter,
    sequences: Option<&SequenceStore>,
) -> Result<(Vec<RepeatsData>, usize)> {
    use std::io::{BufRead, BufReader, BufWriter, Write};
    use std::time::{Duration, Instant};

    for tool in &["blastn", "makeblastdb"] {
        match std::process::Command::new(tool).arg("-version").output() {
            Ok(output) if output.status.success() => {}
            _ => {
                eprintln!("Warning: {} not found on PATH. Skipping MAG inter-contig BLAST.", tool);
                return Ok((Vec::new(), 0));
            }
        }
    }
    let name_to_contig: HashMap<&str, usize> =
        contigs.iter().enumerate().map(|(i, c)| (c.name.as_str(), i)).collect();
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Reference BLAST within MAGs on {} MAGs...", mag_manifest.len());
    eprintln!("  Hit writes: batches of at most 10,000; timing totals are summed worker seconds");

    let started = Instant::now();
    let total_inter_count = AtomicUsize::new(0);
    let mags_done = AtomicUsize::new(0);
    // Sequence loading, FASTA writing, makeblastdb, blastn, parsing, DB wait, DB writing.
    let timing_ns: [AtomicU64; 7] = std::array::from_fn(|_| AtomicU64::new(0));
    let all_repeats: std::sync::Mutex<Vec<RepeatsData>> = std::sync::Mutex::new(Vec::new());

    mag_manifest.par_iter().try_for_each(|(mag_name, member_names)| -> Result<()> {
        let mag_started = Instant::now();
        let t = Instant::now();
        let members: Vec<std::borrow::Cow<'_, ContigInfo>> = member_names.iter()
            .filter_map(|n| name_to_contig.get(n.as_str()).copied())
            .filter(|&i| sequences.map_or(contigs[i].sequence.is_some(), |s| s.has_sequence(i)))
            .map(|i| match sequences {
                Some(store) => store.load(i, &contigs[i]).map(std::borrow::Cow::Owned),
                None => Ok(std::borrow::Cow::Borrowed(&contigs[i])),
            }).collect::<Result<_>>()?;
        let load_time = t.elapsed();
        let mut phases = [Duration::ZERO; 7];
        phases[0] = load_time;
        if !members.is_empty() {
            // One directory owns all BLAST sidecars and output, including on failure.
            let temporary = tempfile::Builder::new().prefix("thebigbam-mag-blast-").tempdir()?;
            let fasta_path = temporary.path().join("members.fasta");
            let hits_path = temporary.path().join("hits.tsv");
            let t = Instant::now();
            {
                let mut fasta = BufWriter::with_capacity(1024 * 1024, fs::File::create(&fasta_path)?);
                for contig in &members {
                    writeln!(fasta, ">{}", contig.name)?;
                    if let Some(sequence) = &contig.sequence {
                        for chunk in sequence.chunks(80) {
                            fasta.write_all(chunk)?;
                            fasta.write_all(b"\n")?;
                        }
                    }
                }
                fasta.flush()?;
            }
            phases[1] = t.elapsed();
            drop(members);

            let t = Instant::now();
            let output = std::process::Command::new("makeblastdb")
                .arg("-in").arg(&fasta_path).arg("-dbtype").arg("nucl")
                .output().with_context(|| format!("Starting makeblastdb for MAG '{mag_name}'"))?;
            anyhow::ensure!(output.status.success(), "makeblastdb failed for MAG '{}': {}",
                mag_name, String::from_utf8_lossy(&output.stderr));
            phases[2] = t.elapsed();

            let t = Instant::now();
            let output = std::process::Command::new("blastn")
                .arg("-query").arg(&fasta_path).arg("-db").arg(&fasta_path)
                .args(BLAST_ARGS).arg("-outfmt").arg("6")
                .arg("-out").arg(&hits_path)
                .output().with_context(|| format!("Starting blastn for MAG '{mag_name}'"))?;
            anyhow::ensure!(output.status.success(), "blastn failed for MAG '{}': {}",
                mag_name, String::from_utf8_lossy(&output.stderr));
            phases[3] = t.elapsed();

            let t = Instant::now();
            let mut repeats = Vec::new();
            let mut inter = Vec::new();
            let mut db_wait = Duration::ZERO;
            let mut db_write = Duration::ZERO;
            let mut flush = |repeats: &mut Vec<RepeatsData>, inter: &mut Vec<(String, String, i32, i32, i32, i32, f64)>| -> Result<()> {
                let (waited, written) = db_writer.write_mag_blast_batch(repeats, inter)
                    .with_context(|| format!("Writing BLAST hits for MAG '{mag_name}'"))?;
                db_wait += waited;
                db_write += written;
                total_inter_count.fetch_add(inter.len(), Ordering::Relaxed);
                inter.clear();
                all_repeats.lock().map_err(|e| anyhow::anyhow!("Repeat collection lock: {e}"))?.extend(repeats.drain(..));
                Ok(())
            };
            for line in BufReader::new(fs::File::open(&hits_path)?).lines() {
                let line = line?;
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 10 { continue; }
                let pident: f64 = f[2].parse().unwrap_or(0.0);
                let qstart: i32 = f[6].parse().unwrap_or(0);
                let qend: i32 = f[7].parse().unwrap_or(0);
                let sstart: i32 = f[8].parse().unwrap_or(0);
                let send: i32 = f[9].parse().unwrap_or(0);
                if f[0] == f[1] {
                    if qstart == sstart && qend == send { continue; }
                    repeats.push(RepeatsData {
                        contig_name: f[0].to_owned(), position1: qstart, position2: qend,
                        position1prime: sstart, position2prime: send, pident,
                        is_direct: (qstart < qend && sstart < send) || (qstart > qend && sstart > send),
                    });
                } else {
                    inter.push((f[0].to_owned(), f[1].to_owned(), qstart, qend, sstart, send, pident));
                }
                if repeats.len() + inter.len() >= 10_000 {
                    flush(&mut repeats, &mut inter)?;
                }
            }
            flush(&mut repeats, &mut inter)?;
            phases[4] = t.elapsed().saturating_sub(db_wait + db_write);
            phases[5] = db_wait;
            phases[6] = db_write;
        }
        for (total, elapsed) in timing_ns.iter().zip(phases) {
            total.fetch_add(elapsed.as_nanos() as u64, Ordering::Relaxed);
        }
        if mag_started.elapsed().as_secs() >= 300 {
            eprintln!("  Slow MAG '{}': load {:.1}s, FASTA {:.1}s, makeblastdb {:.1}s, blastn {:.1}s, parse {:.1}s, DB wait {:.1}s, DB write {:.1}s",
                mag_name, phases[0].as_secs_f64(), phases[1].as_secs_f64(), phases[2].as_secs_f64(),
                phases[3].as_secs_f64(), phases[4].as_secs_f64(), phases[5].as_secs_f64(), phases[6].as_secs_f64());
        }
        let done = mags_done.fetch_add(1, Ordering::Relaxed) + 1;
        if done <= 5 || done % 100 == 0 || done == mag_manifest.len() {
            let seconds = timing_ns.each_ref().map(|v| v.load(Ordering::Relaxed) as f64 / 1e9);
            eprintln!("  Reference MAG BLAST: {}/{} MAGs, {:.1} min elapsed; worker totals: load {:.0}s, FASTA {:.0}s, makeblastdb {:.0}s, blastn {:.0}s, parse {:.0}s, DB wait {:.0}s, DB write {:.0}s",
                done, mag_manifest.len(), started.elapsed().as_secs_f64() / 60.0,
                seconds[0], seconds[1], seconds[2], seconds[3], seconds[4], seconds[5], seconds[6]);
        }
        Ok(())
    })?;

    let reps = all_repeats.into_inner().map_err(|e| anyhow::anyhow!("Repeat collection lock: {e}"))?;
    let inter_count = total_inter_count.load(Ordering::Relaxed);
    eprintln!("MAG BLAST: {} intra-contig repeat hits, {} inter-contig hits", reps.len(), inter_count);
    Ok((reps, inter_count))
}

/// Compute median of event lengths (0 for exact match, clip/insertion length for near-match)
/// across all positions in an area. The event_lengths array already contains the correct
/// values collected during BAM processing.
fn compute_area_median_clippings(
    area: &PeakArea,
    event_lengths: &[Vec<u32>],
) -> f64 {
    let mut lengths: Vec<u32> = Vec::new();
    for pos in area.start_pos..=area.end_pos {
        let idx = (pos - 1) as usize;
        if idx < event_lengths.len() {
            lengths.extend_from_slice(&event_lengths[idx]);
        }
    }
    if lengths.is_empty() {
        return 0.0;
    }
    lengths.sort_unstable();
    let mid = lengths.len() / 2;
    if lengths.len() % 2 == 0 {
        (lengths[mid - 1] as f64 + lengths[mid] as f64) / 2.0
    } else {
        lengths[mid] as f64
    }
}


/// Encode features from FeatureArrays into BLOB format and compute metrics.
/// Returns a tuple of:
/// - Optional packaging data for phagetermini module (includes all termini with filtering metadata)
/// - Optional metric structs (misassembly, microdiversity, side_misassembly, topology)
fn add_features_from_arrays(
    arrays: &mut FeatureArrays,
    contig_name: &str,
    contig_id: i64,
    config: &ProcessConfig,
    is_circular: bool,
    seq_type: SequencingType,
    flags: ModuleFlags,
    primary_count: u64,
    repeats: &[RepeatsData],
    cds_index: Option<&CdsIndex>,
    blob_output: &mut Vec<(String, i64, crate::blob::EncodedBlob)>,
    timing: Option<&TimingAccum>,
) -> (Option<PackagingData>, Option<(MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData)>) {
    use crate::blob::{encode_dense_blob, smooth_dense_values, encode_sparse_blob, EventMeta, MetadataFlags,
                       codon_category_to_id, codon_to_id, aa_to_id};
    use crate::types::{get_encoding, get_value_scale, Encoding, METADATA_STATS_SCALE};
    let encode_dense_blob = |values: &[i32], scale, len| {
        let start = timing.map(|_| std::time::Instant::now());
        let encoded = encode_dense_blob(values, scale, len);
        if let (Some(a), Some(t)) = (timing, start) {
            a.dense_encode_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        encoded
    };
    let encode_sparse_blob = |positions: &[u32], values: &[i32], metadata: Option<&[EventMeta]>, flags, scale, len| {
        let start = timing.map(|_| std::time::Instant::now());
        let encoded = encode_sparse_blob(positions, values, metadata, flags, scale, len);
        if let (Some(a), Some(t)) = (timing, start) {
            a.sparse_encode_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        encoded
    };
    let pt_config = config.phagetermini_config;
    let contig_length = arrays.ref_length();

    // Select candidates from unmodified read evidence before DTR-dependent filters.
    let terminal_evidence;
    let repeats = if flags.phagetermini {
        if let Some(search) = &config.terminal_repeats {
            if crate::terminal_repeats::has_candidate(&arrays.reads_starts, &arrays.reads_ends,
                &arrays.coverage_reduced, &pt_config) {
                let existing: Vec<_> = repeats.iter().filter(|r| r.contig_name == contig_name
                    && r.pident >= pt_config.min_identity_dtr as f64
                    && is_valid_terminal_repeat(r, contig_length, pt_config.max_distance_duplication))
                    .cloned().collect();
                terminal_evidence = search.lookup_for_sample(contig_id, contig_name, &existing, &config.sample_name);
                terminal_evidence.as_slice()
            } else { &[] }
        } else { repeats }
    } else { repeats };

    // Build DTR regions from repeats for this contig (used for both merging and classification)
    // Filter criteria:
    // - Same contig
    // - ≥90% identity
    // - Valid terminal repeat (one region at start, other at end of contig)
    let dtr_regions: Vec<DtrRegion> = if flags.phagetermini && !repeats.is_empty() {
        repeats
            .iter()
            .filter(|d| {
                d.contig_name == contig_name
                    && d.pident >= pt_config.min_identity_dtr as f64
                    && is_valid_terminal_repeat(d, contig_length, pt_config.max_distance_duplication)
            })
            .map(|d| {
                // Determine if direct (DTR) or inverted (ITR)
                let is_direct = (d.position1 < d.position2 && d.position1prime < d.position2prime)
                    || (d.position1 > d.position2 && d.position1prime > d.position2prime);

                // Determine first and second regions (first = lower start position)
                if d.position1 < d.position1prime {
                    DtrRegion {
                        first_start: d.position1.min(d.position2),
                        first_end: d.position1.max(d.position2),
                        second_start: d.position1prime.min(d.position2prime),
                        second_end: d.position1prime.max(d.position2prime),
                        is_direct,
                    }
                } else {
                    DtrRegion {
                        first_start: d.position1prime.min(d.position2prime),
                        first_end: d.position1prime.max(d.position2prime),
                        second_start: d.position1.min(d.position2),
                        second_end: d.position1.max(d.position2),
                        is_direct,
                    }
                }
            })
            .collect()
    } else {
        Vec::new()
    };

    // Coverage (always compress self-referentially)
    let primary_reads_f64: Vec<f64> = arrays.primary_reads.iter().map(|&x| x as f64).collect();
    let cid = contig_id;
    let clen = contig_length as u32;
    let bar_threshold = config.bar_ratio * 0.01;
    let min_occ = config.min_occurrences;
    if flags.coverage {
        // MAPQ - average mapping quality per position (needed for blob encoding below)
        let mapq_f64: Vec<f64> = arrays.sum_mapq.iter()
            .zip(&arrays.primary_reads)
            .map(|(&sum, &count)| if count > 0 { sum as f64 / count as f64 } else { 0.0 })
            .collect();

        // === BLOB encoding for dense coverage features ===
        // primary_reads: raw i32
        let mut pr_i32: Vec<i32> = arrays.primary_reads.iter().map(|&x| x as i32).collect();
        smooth_dense_values(&mut pr_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("primary_reads"), Encoding::Dense);
        blob_output.push(("primary_reads".into(), cid, encode_dense_blob(&pr_i32, get_value_scale("primary_reads"), clen)));

        // plus/minus strand
        let mut pp_i32: Vec<i32> = arrays.primary_reads_plus_only.iter().map(|&x| x as i32).collect();
        smooth_dense_values(&mut pp_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("primary_reads_plus_only"), Encoding::Dense);
        blob_output.push(("primary_reads_plus_only".into(), cid, encode_dense_blob(&pp_i32, get_value_scale("primary_reads_plus_only"), clen)));
        let mut pm_i32: Vec<i32> = arrays.primary_reads_minus_only.iter().map(|&x| x as i32).collect();
        smooth_dense_values(&mut pm_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("primary_reads_minus_only"), Encoding::Dense);
        blob_output.push(("primary_reads_minus_only".into(), cid, encode_dense_blob(&pm_i32, get_value_scale("primary_reads_minus_only"), clen)));

        // secondary, supplementary
        let mut sec_i32: Vec<i32> = arrays.secondary_reads.iter().map(|&x| x as i32).collect();
        smooth_dense_values(&mut sec_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("secondary_reads"), Encoding::Dense);
        blob_output.push(("secondary_reads".into(), cid, encode_dense_blob(&sec_i32, get_value_scale("secondary_reads"), clen)));
        let mut sup_i32: Vec<i32> = arrays.supplementary_reads.iter().map(|&x| x as i32).collect();
        smooth_dense_values(&mut sup_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("supplementary_reads"), Encoding::Dense);
        blob_output.push(("supplementary_reads".into(), cid, encode_dense_blob(&sup_i32, get_value_scale("supplementary_reads"), clen)));

        let mapq_scale = get_value_scale("mapq").multiplier();
        let mut mapq_i32: Vec<i32> = mapq_f64.iter().map(|&x| (x * mapq_scale).round() as i32).collect();
        smooth_dense_values(&mut mapq_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("mapq"), Encoding::Dense);
        blob_output.push(("mapq".into(), cid, encode_dense_blob(&mapq_i32, get_value_scale("mapq"), clen)));
    }

    // Clipping statistics (shared between assemblycheck run-length encoding and sparse blob encoding)
    let left_clip_counts: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| v.len() as f64).collect();
    let left_clip_means: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
        if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
    }).collect();
    let left_clip_medians: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
        if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
    }).collect();
    let left_clip_stds: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
        if v.len() <= 1 { 0.0 } else {
            let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
            let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
            var.sqrt()
        }
    }).collect();
    let right_clip_counts: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| v.len() as f64).collect();
    let right_clip_means: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
        if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
    }).collect();
    let right_clip_medians: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
        if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
    }).collect();
    let right_clip_stds: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
        if v.len() <= 1 { 0.0 } else {
            let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
            let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
            var.sqrt()
        }
    }).collect();

    let mut left_clip_runs: Vec<Run> = Vec::new();
    let mut right_clip_runs: Vec<Run> = Vec::new();
    if flags.mapping_metrics || flags.phagetermini {
        left_clip_runs = add_compressed_feature_with_stats(&left_clip_counts, &left_clip_means, &left_clip_medians, &left_clip_stds,
            Some(&primary_reads_f64), "left_clippings", config);
        right_clip_runs = add_compressed_feature_with_stats(&right_clip_counts, &right_clip_means, &right_clip_medians, &right_clip_stds,
            Some(&primary_reads_f64), "right_clippings", config);
    }
            
    let sparse_scale = get_value_scale("mismatches").multiplier();

    let filter_sparse = |counts: &[f64], coverage: &[f64], scale: f64| -> (Vec<u32>, Vec<i32>) {
        let mut positions = Vec::new();
        let mut values = Vec::new();
        let n = counts.len().min(coverage.len());
        for i in 0..n {
            let val = counts[i];
            let cov = coverage[i];
            if val >= cov * bar_threshold && val >= min_occ as f64 {
                positions.push(i as u32);
                values.push((val / cov * scale).round() as i32);
            }
        }
        (positions, values)
    };

    if flags.mapping_metrics {
        let mismatches_f64: Vec<f64> = arrays.mismatches.iter().map(|&x| x as f64).collect();
        let deletions_f64: Vec<f64> = arrays.deletions.iter().map(|&x| x as f64).collect();
        let insertion_counts: Vec<f64> = arrays.insertion_lengths.iter().map(|v| v.len() as f64).collect();
        let insertion_means: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
        }).collect();
        let insertion_medians: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();
        let insertion_stds: Vec<f64> = arrays.insertion_lengths.iter().map(|v| {
            if v.len() <= 1 { 0.0 } else {
                let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                var.sqrt()
            }
        }).collect();

        // --- Compute dominant bases/sequences (used by blob encoding below) ---
        let dominant_mismatches = arrays.compute_dominant_mismatch_bases(&arrays.primary_reads, bar_threshold);
        let dominant_insertions = FeatureArrays::compute_dominant_sequences(&arrays.insertion_sequences, &arrays.primary_reads, bar_threshold);
        let dominant_left_clips = FeatureArrays::compute_dominant_sequences(&arrays.left_clip_sequences, &arrays.primary_reads, bar_threshold);
        let dominant_right_clips = FeatureArrays::compute_dominant_sequences(&arrays.right_clip_sequences, &arrays.primary_reads, bar_threshold);

        // Compute codon changes from position-level mismatch summaries
        let dominant_codons = compute_codon_changes_from_summaries(
            &arrays.mismatch_base_counts, &arrays.primary_reads,
            cds_index, contig_length, bar_threshold,
        );

        // mismatches (with sequence + codons)
        {
            let (pos, vals) = filter_sparse(&mismatches_f64, &primary_reads_f64, sparse_scale);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta::default();
                if idx < dominant_mismatches.len() {
                    let (base, pct) = dominant_mismatches[idx];
                    if base != 0 {
                        em.sequence = Some(vec![base]);
                        em.prevalence = Some(pct as i16);
                    }
                }
                if let Some((cat, codon, aa)) = dominant_codons.get(&idx) {
                    em.codon_category = Some(codon_category_to_id(cat));
                    em.codon_id = Some(codon_to_id(codon));
                    em.aa_id = Some(aa_to_id(aa));
                }
                em
            }).collect();
            let flags = MetadataFlags { sparse: true, has_stats: false, has_sequence: true, has_codons: true, has_partner: false };
            debug_assert_eq!(get_encoding("mismatches"), Encoding::Sparse);
            blob_output.push(("mismatches".into(), cid,
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, get_value_scale("mismatches"), clen)));

        }

        // deletions (value only)
        {
            let (pos, vals) = filter_sparse(&deletions_f64, &primary_reads_f64, sparse_scale);
            let flags = MetadataFlags { sparse: true, ..Default::default() };
            debug_assert_eq!(get_encoding("deletions"), Encoding::Sparse);
            blob_output.push(("deletions".into(), cid,
                encode_sparse_blob(&pos, &vals, None, flags, get_value_scale("deletions"), clen)));

        }

        // insertions (with stats + sequence)
        {
            let (pos, vals) = filter_sparse(&insertion_counts, &primary_reads_f64, sparse_scale);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((insertion_means.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    median: Some((insertion_medians.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    std: Some((insertion_stds.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    ..Default::default()
                };
                if let Some((seq_str, pct)) = dominant_insertions.get(&idx) {
                    em.sequence = Some(seq_str.as_bytes().to_vec());
                    em.prevalence = Some(*pct as i16);
                }
                em
            }).collect();
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false, has_partner: false };
            debug_assert_eq!(get_encoding("insertions"), Encoding::Sparse);
            blob_output.push(("insertions".into(), cid,
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, get_value_scale("insertions"), clen)));

        }

        // left_clippings (with stats + sequence) — reuses clip stats computed above
        {
            let (pos, vals) = filter_sparse(&left_clip_counts, &primary_reads_f64, sparse_scale);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((left_clip_means.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    median: Some((left_clip_medians.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    std: Some((left_clip_stds.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    ..Default::default()
                };
                if let Some((seq_str, pct)) = dominant_left_clips.get(&idx) {
                    em.sequence = Some(seq_str.as_bytes().to_vec());
                    em.prevalence = Some(*pct as i16);
                }
                em
            }).collect();
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false, has_partner: false };
            debug_assert_eq!(get_encoding("left_clippings"), Encoding::Sparse);
            blob_output.push(("left_clippings".into(), cid,
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, get_value_scale("left_clippings"), clen)));
        }

        // right_clippings (with stats + sequence) — reuses clip stats computed above
        {
            let (pos, vals) = filter_sparse(&right_clip_counts, &primary_reads_f64, sparse_scale);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((right_clip_means.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    median: Some((right_clip_medians.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    std: Some((right_clip_stds.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    ..Default::default()
                };
                if let Some((seq_str, pct)) = dominant_right_clips.get(&idx) {
                    em.sequence = Some(seq_str.as_bytes().to_vec());
                    em.prevalence = Some(*pct as i16);
                }
                em
            }).collect();
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false, has_partner: false };
            debug_assert_eq!(get_encoding("right_clippings"), Encoding::Sparse);
            blob_output.push(("right_clippings".into(), cid,
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, get_value_scale("right_clippings"), clen)));
        }
    }

    // RNA module
    if flags.rna {
        // splicings — absolute per-base count of CIGAR 'N' spans
        // Stored as raw counts (not coverage-relative) because N-correction
        // decrements coverage at intronic positions, making the ratio meaningless.
        let mut spl_pos = Vec::new();
        let mut spl_vals = Vec::new();
        for (i, &count) in arrays.splices.iter().enumerate() {
            if count > 0 {
                spl_pos.push(i as u32);
                spl_vals.push(count as i32);
            }
        }
        let flags = MetadataFlags { sparse: true, ..Default::default() };
        debug_assert_eq!(get_encoding("splicings"), Encoding::Sparse);
        blob_output.push(("splicings".into(), cid,
            encode_sparse_blob(&spl_pos, &spl_vals, None, flags, get_value_scale("splicings"), clen)));

    }

    // Paired-reads module
    if flags.paired_read_metrics && seq_type.is_short_paired() {
        let non_inward_f64: Vec<f64> = arrays.non_inward_pairs.iter().map(|&x| x as f64).collect();
        let mate_unmapped_f64: Vec<f64> = arrays.mate_not_mapped.iter().map(|&x| x as f64).collect();
        let mate_other_contig_f64: Vec<f64> = arrays.mate_on_another_contig.iter().map(|&x| x as f64).collect();

        // Insert sizes (curve for paired reads)
        let values: Vec<f64> = arrays
            .sum_insert_sizes
            .iter()
            .zip(&arrays.count_insert_sizes)
            .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
            .collect();

        // === BLOB encoding for paired-reads features ===
        let sparse_scale_pr = get_value_scale("non_inward_pairs").multiplier();

        let (pos, vals) = filter_sparse(&non_inward_f64, &primary_reads_f64, sparse_scale_pr);
        let sp_flags = MetadataFlags { sparse: true, ..Default::default() };
        debug_assert_eq!(get_encoding("non_inward_pairs"), Encoding::Sparse);
        blob_output.push(("non_inward_pairs".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("non_inward_pairs"), clen)));

        let (pos, vals) = filter_sparse(&mate_unmapped_f64, &primary_reads_f64, sparse_scale_pr);
        debug_assert_eq!(get_encoding("mate_not_mapped"), Encoding::Sparse);
        blob_output.push(("mate_not_mapped".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("mate_not_mapped"), clen)));

        let (pos, vals) = filter_sparse(&mate_other_contig_f64, &primary_reads_f64, sparse_scale_pr);
        debug_assert_eq!(get_encoding("mate_on_another_contig"), Encoding::Sparse);
        blob_output.push(("mate_on_another_contig".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("mate_on_another_contig"), clen)));


        // insert_sizes: dense curve
        let is_scale = get_value_scale("insert_sizes").multiplier();
        let mut is_i32: Vec<i32> = values.iter().map(|&x| (x * is_scale).round() as i32).collect();
        smooth_dense_values(&mut is_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("insert_sizes"), Encoding::Dense);
        blob_output.push(("insert_sizes".into(), cid, encode_dense_blob(&is_i32, get_value_scale("insert_sizes"), clen)));
    }

    // Long-reads module
    if flags.long_read_metrics && seq_type.is_long() {
        let values: Vec<f64> = arrays
            .sum_read_lengths
            .iter()
            .zip(&arrays.count_read_lengths)
            .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
            .collect();

        // === BLOB encoding for long-reads features ===
        let rl_scale = get_value_scale("read_lengths").multiplier();
        let mut rl_i32: Vec<i32> = values.iter().map(|&x| (x * rl_scale).round() as i32).collect();
        smooth_dense_values(&mut rl_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("read_lengths"), Encoding::Dense);
        blob_output.push(("read_lengths".into(), cid, encode_dense_blob(&rl_i32, get_value_scale("read_lengths"), clen)));
    }

    // Phagetermini features
    let packaging_result = if flags.phagetermini {
        // === STEP 1: Save ORIGINAL data to database (before any DTR merging) ===
        let reads_starts_original: Vec<f64> = arrays.reads_starts.iter().map(|&x| x as f64).collect();
        let reads_ends_original: Vec<f64> = arrays.reads_ends.iter().map(|&x| x as f64).collect();
        let start_evt_medians: Vec<f64> = arrays.start_event_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();
        let end_evt_medians: Vec<f64> = arrays.end_event_lengths.iter().map(|v| {
            if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
        }).collect();

        // Compute dominant clip sequences (used by blob encoding below)
        let dominant_start_clips = FeatureArrays::compute_dominant_sequences(&arrays.start_clip_sequences, &arrays.coverage_reduced, bar_threshold);
        let dominant_end_clips = FeatureArrays::compute_dominant_sequences(&arrays.end_clip_sequences, &arrays.coverage_reduced, bar_threshold);

        // === BLOB encoding for phagetermini features ===
        {
            // coverage_reduced: dense curve
            let mut cr_i32: Vec<i32> = arrays.coverage_reduced.iter().map(|&x| x as i32).collect();
            smooth_dense_values(&mut cr_i32, config.variation_percentage);
            debug_assert_eq!(get_encoding("coverage_reduced"), Encoding::Dense);
            blob_output.push(("coverage_reduced".into(), cid, encode_dense_blob(&cr_i32, get_value_scale("coverage_reduced"), clen)));

            // reads_starts: sparse with median + sequence
            let cov_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
            let rs_scale = get_value_scale("reads_starts").multiplier();
            let mut rs_pos = Vec::new();
            let mut rs_vals = Vec::new();
            for i in 0..reads_starts_original.len().min(cov_reduced_f64.len()) {
                let val = reads_starts_original[i];
                let cov = cov_reduced_f64[i];
                if val >= cov * bar_threshold && val >= min_occ as f64 {
                    rs_pos.push(i as u32);
                    rs_vals.push((val / cov * rs_scale).round() as i32);
                }
            }
            let rs_meta: Vec<EventMeta> = rs_pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    median: Some((start_evt_medians.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    ..Default::default()
                };
                if let Some((seq_str, pct)) = dominant_start_clips.get(&idx) {
                    em.sequence = Some(seq_str.as_bytes().to_vec());
                    em.prevalence = Some(*pct as i16);
                }
                em
            }).collect();
            let pt_flags = MetadataFlags { sparse: true, has_stats: false, has_sequence: true, has_codons: false, has_partner: false };
            debug_assert_eq!(get_encoding("reads_starts"), Encoding::Sparse);
            blob_output.push(("reads_starts".into(), cid,
                encode_sparse_blob(&rs_pos, &rs_vals, Some(&rs_meta), pt_flags, get_value_scale("reads_starts"), clen)));


            // reads_ends: sparse with median + sequence
            let re_scale = get_value_scale("reads_ends").multiplier();
            let mut re_pos = Vec::new();
            let mut re_vals = Vec::new();
            for i in 0..reads_ends_original.len().min(cov_reduced_f64.len()) {
                let val = reads_ends_original[i];
                let cov = cov_reduced_f64[i];
                if val >= cov * bar_threshold && val >= min_occ as f64 {
                    re_pos.push(i as u32);
                    re_vals.push((val / cov * re_scale).round() as i32);
                }
            }
            let re_meta: Vec<EventMeta> = re_pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    median: Some((end_evt_medians.get(idx).copied().unwrap_or(0.0) * METADATA_STATS_SCALE).round() as i32),
                    ..Default::default()
                };
                if let Some((seq_str, pct)) = dominant_end_clips.get(&idx) {
                    em.sequence = Some(seq_str.as_bytes().to_vec());
                    em.prevalence = Some(*pct as i16);
                }
                em
            }).collect();
            debug_assert_eq!(get_encoding("reads_ends"), Encoding::Sparse);
            blob_output.push(("reads_ends".into(), cid,
                encode_sparse_blob(&re_pos, &re_vals, Some(&re_meta), pt_flags, get_value_scale("reads_ends"), clen)));

        }

        // === STEP 2: Calculate aligned fraction and global metrics ===
        let aligned_count = arrays.coverage_reduced.iter().filter(|&&x| x > 0).count();
        let aligned_fraction = if arrays.coverage_reduced.is_empty() {
            0.0
        } else {
            (aligned_count as f64 / arrays.coverage_reduced.len() as f64) * 100.0
        };

        if aligned_fraction >= pt_config.min_aligned_fraction as f64 {
            // Compute clipped_ratio = (number of primary reads - reads with clean termini) / reads with clean termini
            // This represents the fraction of reads that are clipped relative to clean reads.
            // For long reads, each read is split in two (start half and end half), so
            // clean_reads_count already counts each clean terminus separately.
            // We double primary_count to match (each read = 2 potential termini).
            let clipped_ratio = if arrays.clean_reads_count > 0 {
                let effective_primary = if seq_type.is_long() {
                    primary_count as f64 * 2.0
                } else {
                    primary_count as f64
                };
                ((effective_primary - arrays.clean_reads_count as f64) / arrays.clean_reads_count as f64).max(0.0)
            } else {
                1.0
            };

            // Get clipping counts as u64 arrays, filtering by min_clipping_length
            let min_clip_len = pt_config.min_clipping_length;
            let left_clip_counts: Vec<u64> = arrays.left_clipping_lengths.iter()
                .map(|v| v.iter().filter(|&&len| len >= min_clip_len).count() as u64)
                .collect();
            let right_clip_counts: Vec<u64> = arrays.right_clipping_lengths.iter()
                .map(|v| v.iter().filter(|&&len| len >= min_clip_len).count() as u64)
                .collect();

            // === STEP 3-5: Filter, statistical test, merge into areas, global filter ===
            // Returns ALL areas (both kept and discarded) with filtering metadata
            let (mut start_areas, mut end_areas) = filter_and_merge_to_areas_with_diagnostics(
                &arrays.reads_starts,
                &arrays.reads_ends,
                &left_clip_counts,
                &right_clip_counts,
                &arrays.coverage_reduced,
                &arrays.primary_reads,
                clipped_ratio,
                &pt_config,
                contig_length,
                is_circular,
                &dtr_regions,
            );

            // Compute median clipping lengths per area using event lengths
            for area in &mut start_areas {
                area.median_clippings = compute_area_median_clippings(area, &arrays.start_event_lengths);
            }
            for area in &mut end_areas {
                area.median_clippings = compute_area_median_clippings(area, &arrays.end_event_lengths);
            }

            // === STEP 6: Classify packaging based on peak-areas ===
            let (mechanism, left_termini, right_termini, duplication, repeat_length, median_left_clippings, median_right_clippings) = classify_packaging_areas(
                &start_areas,
                &end_areas,
                contig_length,
                is_circular,
                &dtr_regions,
                pt_config.max_distance_peaks,
            );

            // Only return PackagingData if there's a detected mechanism (not "No_packaging")
            // All termini (both kept and discarded) are included with filtering metadata
            if mechanism != "No_packaging" {
                Some(PackagingData {
                    contig_id,
                    contig_name: contig_name.to_string(),
                    mechanism,
                    left_termini,
                    right_termini,
                    duplication,
                    repeat_length,
                    median_left_termini_clippings: median_left_clippings,
                    median_right_termini_clippings: median_right_clippings,
                })
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    // Compute all metrics when mapping_metrics is enabled
    let metrics_result = if flags.mapping_metrics {
        let (misassembly, microdiversity, side_misassembly, topology) = compute_all_metrics(
            &arrays.left_clipping_lengths,
            &arrays.right_clipping_lengths,
            &arrays.insertion_lengths,
            &arrays.deletion_lengths,
            &arrays.primary_reads,
            &arrays.mismatches,
            &arrays.deletions,
            contig_name,
            contig_id,
            contig_length,
            arrays.circularising_reads_count,
            arrays.circularising_inserts_count,
            &arrays.circularising_insert_sizes,
            &arrays.all_proper_insert_sizes,
            arrays.contig_end_mates_mapped_on_another_contig,
            arrays.circularising_confirmed,
            &arrays.circularising_min_overlaps,
            is_circular,
            &left_clip_runs,
            &right_clip_runs,
        );
        Some((misassembly, microdiversity, side_misassembly, topology))
    } else {
        None
    };

    (packaging_result, metrics_result)
}

/// Keep MAG coordinates independent of BAM access order.
fn mag_members_in_bam_order(group: &MagJobGroup) -> Vec<&MagContigJob> {
    let mut members: Vec<_> = group.members.iter().collect();
    members.sort_unstable_by_key(|member| member.tid);
    members
}

/// Process one MAG group: all member contigs in sequence with one BAM reader.
/// Keeps primary coverage in MAG coordinates for summary statistics. Other
/// feature arrays are encoded per contig, without duplicate MAG accumulators.
fn process_mag_group(
    bam: &mut IndexedReader,
    group: &MagJobGroup,
    seq_type: SequencingType,
    flags: ModuleFlags,
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    cds_indexes: &dyn CdsSource,
    accum_ref: Option<&TimingAccum>,
) -> Result<Option<MagGroupResult>> {
    use crate::blob::EncodedBlob;

    let mag_len = group.mag_length as usize;
    if mag_len == 0 || group.members.is_empty() { return Ok(None); }

    // Read each member once; MAG thresholds are applied to the completed results.
    // Only primary coverage is needed for MAG stats. All other per-position
    // features are encoded from contig arrays, so avoid duplicate MAG arrays.
    let mut mag_primary_reads = vec![0i32; mag_len];
    let mut contig_blobs: Vec<(String, i64, EncodedBlob)> = Vec::new();
    let mut presences = Vec::new();
    let mut all_packaging = Vec::new();
    let mut all_misassembly = Vec::new();
    let mut all_microdiversity = Vec::new();
    let mut all_side_misassembly = Vec::new();
    let mut all_topology = Vec::new();
    let mut mapped_reads: u64 = 0;
    // Preserve offset-order boundaries for the existing MAG roughness statistic.
    let contig_lengths_ordered: Vec<u32> = group.members.iter().map(|m| m.ref_length as u32).collect();

    for member in mag_members_in_bam_order(group) {
        let off = member.offset_in_mag as usize;
        let clen = member.ref_length;


        // Process contig streaming (no per-contig AF filtering in MAG mode)
        let ts = accum_ref.map(|_| std::time::Instant::now());
        let mut phase_t: Option<BamPhaseTimings> = accum_ref.map(|_| BamPhaseTimings::default());
        let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming_counted(
            bam, member.tid, &member.ref_name, member.ref_length, seq_type, flags,
            is_circular, 0.0, config.phagetermini_config.min_clipping_length, &mut phase_t,
        ).map(|(result,count)| { config.primary_reads_seen.fetch_add(count,Ordering::Relaxed); result }) {
            Ok(Some(result)) => result,
            Ok(None) => {
                if let (Some(a), Some(pt)) = (accum_ref, phase_t.as_ref()) {
                    a.stream_fetch_seek_ns.fetch_add(pt.fetch_seek_ns, Ordering::Relaxed);
                    a.stream_records_ns.fetch_add(pt.records_ns, Ordering::Relaxed);
                    a.stream_process_read_ns.fetch_add(pt.process_read_ns, Ordering::Relaxed);
                    a.stream_finalize_ns.fetch_add(pt.finalize_ns, Ordering::Relaxed);
                }
                if let (Some(a), Some(t)) = (accum_ref, ts) {
                    a.streaming_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                continue;
            }
            Err(e) => {
                return Err(e).with_context(||format!("Processing contig {} in MAG {}",member.ref_name,group.mag_name));
            }
        };
        if let (Some(a), Some(pt)) = (accum_ref, phase_t.as_ref()) {
            a.stream_fetch_seek_ns.fetch_add(pt.fetch_seek_ns, Ordering::Relaxed);
            a.stream_records_ns.fetch_add(pt.records_ns, Ordering::Relaxed);
            a.stream_process_read_ns.fetch_add(pt.process_read_ns, Ordering::Relaxed);
            a.stream_finalize_ns.fetch_add(pt.finalize_ns, Ordering::Relaxed);
        }
        if let (Some(a), Some(t)) = (accum_ref, ts) {
            a.streaming_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }

        let end = (off + clen).min(mag_len);
        let n = end - off;
        for i in 0..n.min(arrays.primary_reads.len()) {
            mag_primary_reads[off + i] = arrays.primary_reads[i] as i32;
        }

        let coverage_mean = arrays.coverage_mean() as f32;
        presences.push(build_presence(&arrays, member.contig_id, &member.ref_name,
            coverage_pct, primary_count, coverage_mean));
        mapped_reads += primary_count;

        let ts = accum_ref.map(|_| std::time::Instant::now());
        let mut blob_output: Vec<(String, i64, EncodedBlob)> = Vec::new();
        let cds_index = if flags.mapping_metrics && needs_codon_annotations(&arrays, config.bar_ratio * 0.01) {
            cds_indexes.lookup(member.contig_id)
        } else { None };
        let (pkg, metrics) = add_features_from_arrays(
            &mut arrays, &member.ref_name, member.contig_id, config, is_circular, seq_type, flags,
            primary_count, repeats, cds_index.as_deref(), &mut blob_output, accum_ref,
        );
        crate::features::recycle_feature_arrays(arrays);
        contig_blobs.extend(blob_output);
        if let Some(p) = pkg { all_packaging.push(p); }
        if let Some((mis, micro, side, topo)) = metrics {
            all_misassembly.push(mis);
            all_microdiversity.push(micro);
            all_side_misassembly.push(side);
            all_topology.push(topo);
        }
        if let (Some(a), Some(t)) = (accum_ref, ts) {
            a.feature_calc_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
    }

    if presences.is_empty() {
        return Ok(None);
    }

    // --- MAG coverage stats from raw primary_reads accumulator ---
    let total_reads_in_mag: u64 = presences.iter().map(|p| p.read_count).sum();
    let mag_cov_stats = crate::mag_blob::compute_mag_coverage_stats_from_raw(
        &mag_primary_reads, group.mag_length, &contig_lengths_ordered, total_reads_in_mag,
    );

    Ok(Some(MagGroupResult {
        mag_name: group.mag_name.clone(),
        mag_id: group.mag_id,
        contig_blobs,
        mag_cov_stats,
        presences,
        packaging: all_packaging,
        misassembly: all_misassembly,
        microdiversity: all_microdiversity,
        side_misassembly: all_side_misassembly,
        topology: all_topology,
        mapped_reads,
    }))
}

/// Process one BAM file using streaming (single-pass, optimized).
/// Parallelizes at the contig level for samples with many contigs.
/// Returns (feature_blobs, presences, packaging, misassembly, microdiversity, side_misassembly, topology, sample_name, seq_type, total_reads, mapped_reads)
/// Per-sample result including optional MAG-level blobs.
pub struct SampleProcessResult {
    pub feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    pub presences: Vec<PresenceData>,
    pub packaging: Vec<PackagingData>,
    pub misassembly: Vec<MisassemblyData>,
    pub microdiversity: Vec<MicrodiversityData>,
    pub side_misassembly: Vec<SideMisassemblyData>,
    pub topology: Vec<TopologyData>,
    pub sample_name: String,
    pub seq_type: SequencingType,
    pub total_reads: u64,
    /// Primary reads from processed contigs (MAG filtering occurs separately).
    pub mapped_reads: u64,
    /// Total mapped reads from BAI index — non-zero even when all contigs are filtered out.
    pub bai_mapped_reads: u64,
    pub is_circular: bool,
    pub timings: Option<SampleTimings>,
    /// MAG-level coverage stats, ready to write. Empty in contig mode.
    pub mag_results: Vec<(String, i64, crate::mag_blob::MagCoverageStats)>,
}

pub fn process_sample(
    bam_path: &Path,
    flags: ModuleFlags,
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    contig_by_name: &dyn ContigLookup,
    cds_indexes: &dyn CdsSource,
    mag_member_lookup: Option<&MagMemberLookup>,
) -> Result<SampleProcessResult> {
    process_sample_with_mag_sink(bam_path, flags, config, is_circular, repeats,
        contig_by_name, cds_indexes, mag_member_lookup, None)
}

fn process_sample_with_mag_sink(
    bam_path: &Path,
    flags: ModuleFlags,
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    contig_by_name: &dyn ContigLookup,
    cds_indexes: &dyn CdsSource,
    mag_member_lookup: Option<&MagMemberLookup>,
    mag_sink: Option<&(dyn Fn(MagStreamMessage) -> Result<()> + Sync)>,
) -> Result<SampleProcessResult> {
    if let Some(search) = &config.terminal_repeats { search.check_error()?; }
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "");

    let mut sample_config = config.clone();
    sample_config.sample_name = sample_name.clone();
    sample_config.primary_reads_seen = Arc::new(AtomicU64::new(0));
    let config = &sample_config;
    let time_on = config.enable_timing;
    let mut st = if time_on { Some(SampleTimings::default()) } else { None };
    if let Some(s) = st.as_mut() {
        s.sample_name = sample_name.clone();
        s.bam_path = bam_path.display().to_string();
    }

    // --- Phase: header_read ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Open the BAM once in the single-threaded prelude. We keep it alive through
    // the pre-filter phase below so the par_iter never needs to re-read the
    // header — a reopen costs ~50-200 ms for BAMs with many references.
    let reader_rss_before = get_rss_mb();
    if let Some(timing) = st.as_mut() { timing.rss_before_reader_open_mib = reader_rss_before; }
    let mut temp_bam = IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;
    let n_refs = temp_bam.header().target_count();

    // Determine sequencing type: use provided value or auto-detect from this BAM file
    let seq_type = match config.sequencing_type {
        Some(stype) => stype,
        None => detect_sequencing_type(bam_path)
            .with_context(|| format!("Failed to detect sequencing type from {}", bam_path.display()))?,
    };

    if let (Some(s), Some(t)) = (st.as_mut(), t0) {
        s.header_read = t.elapsed();
        s.n_contigs_in_bam = n_refs as usize;
    }

    if let Some(sink) = mag_sink {
        sink(MagStreamMessage::Start { sample_name: sample_name.clone(), seq_type, is_circular })?;
    }

    // --- Phase: unmapped_read_count + per-tid mapped counts ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Get per-tid stats from the BAM index in one call (avoids a second BAM open).
    // Each entry: (tid, length, mapped_records, unmapped_records).
    let idx_stats = temp_bam.index_stats()
        .with_context(|| format!("Failed to get index stats from: {}", bam_path.display()))?;
    let unmapped_reads: u64 = idx_stats.iter().map(|(_, _, _, unmapped)| *unmapped).sum();
    let mut mapped_per_tid = vec![0u64;n_refs as usize];
    for &(tid,_,mapped,_) in &idx_stats { if tid>=0 { if let Some(slot)=mapped_per_tid.get_mut(tid as usize) { *slot=mapped; } } }
    let bai_mapped_reads: u64 = mapped_per_tid.iter().sum();
    let measured_reader_mib = (get_rss_mb() - reader_rss_before).max(0.0) as u64;

    if let (Some(s), Some(t)) = (st.as_mut(), t0) {
        s.total_read_count = t.elapsed();
    }

    // --- Phase: filter_tids (single-threaded) ---
    // Walk the BAM header once and materialise only the tids whose reference
    // name is also a known contig AND that have at least one mapped read.
    // Skipping empty contigs avoids millions of wasted rayon jobs for large
    // metagenomic datasets where most contigs have 0 reads in a given sample.
    let t_filter = if time_on { Some(std::time::Instant::now()) } else { None };
    let jobs: Vec<(u32, i64, String, usize)> = (0..n_refs)
        .filter_map(|tid| {
            if mapped_per_tid.get(tid as usize).copied().unwrap_or(0) == 0 {
                return None;
            }
            let ref_name = std::str::from_utf8(temp_bam.header().tid2name(tid)).ok()?;
            let (_, contig_id) = contig_by_name.find(ref_name)?;
            let ref_length = temp_bam.header().target_len(tid).unwrap_or(0) as usize;
            Some((tid, contig_id, ref_name.to_string(), ref_length))
        })
        .collect();
    if jobs.is_empty() && n_refs > 0 {
        eprintln!("WARNING: Sample '{}': 0/{} BAM references matched known contigs. \
                   Check that contig names in your annotation/assembly files match the BAM reference names.",
                  sample_name, n_refs);
    }
    for tid in 0..n_refs {
        if mapped_per_tid[tid as usize] == 0 { continue; }
        let name = std::str::from_utf8(temp_bam.header().tid2name(tid))?;
        if contig_by_name.find(name).is_none() || mag_member_lookup.map_or(false,|lookup|!lookup.contains_key(name)) {
            let length=temp_bam.header().target_len(tid).unwrap_or(0);
            let count=crate::bam_reader::count_reference_primary_reads(&mut temp_bam,tid,length)?;
            config.primary_reads_seen.fetch_add(count,Ordering::Relaxed);
        }
    }
    drop(temp_bam);
    drop(idx_stats);
    drop(mapped_per_tid);
    if let (Some(s), Some(t)) = (st.as_mut(), t_filter) {
        s.filter_tids = t.elapsed();
    }

    // --- Phase: par_iter ---
    let accum = if time_on { Some(TimingAccum::new()) } else { None };
    let accum_ref = accum.as_ref();
    let streamed_contigs = AtomicUsize::new(0);
    let stream_error = std::sync::Mutex::new(None);
    let streamed_contigs_ref = &streamed_contigs;
    let stream_error_ref = &stream_error;
    let t_par = if time_on { Some(std::time::Instant::now()) } else { None };

    // Headers/indexes are replicated per reader. Reserve half the remaining
    // headroom for feature arrays, result batches, and writer allocations.
    let (total_gb, _) = crate::db::detect_available_memory_gb();
    let reader_workers = bam_reader_worker_limit(config.threads, n_refs as u64,
        measured_reader_mib, total_gb * 1024, get_rss_mb() as u64);
    if let Some(timing) = st.as_mut() {
        timing.bam_reader_workers = reader_workers;
        timing.bam_reader_estimated_mib = measured_reader_mib.max((n_refs as u64 * 512).div_ceil(1024 * 1024));
        timing.bam_reader_start_rss_mib = get_rss_mb() as u64;
    }
    {
        eprintln!("Sample '{}': {} BAM readers for {} calculation threads (estimated header/index {} MiB per reader)",
            sample_name, reader_workers, config.threads, measured_reader_mib.max((n_refs as u64 * 512).div_ceil(1024 * 1024)));
    }

    // MAG mode: group jobs by MAG, process MAG-by-MAG (one BAM open per chunk).
    // Contig mode: par_chunks over individual jobs (unchanged).
    let (all_feature_blobs, all_presences, all_packaging, all_misassembly, all_microdiversity, all_side_misassembly, all_topology, _retained_mapped_reads, mag_results) = if let Some(lookup) = mag_member_lookup {
        // Group jobs into MagJobGroups
        let mut mag_groups_map: HashMap<String, MagJobGroup> = HashMap::new();
        for (tid, contig_id, ref_name, ref_length) in &jobs {
            if let Some((mag_name, mag_id, offset, mag_length)) = lookup.get(ref_name.as_str()) {
                let group = mag_groups_map.entry(mag_name.clone()).or_insert_with(|| MagJobGroup {
                    mag_name: mag_name.clone(),
                    mag_id: *mag_id,
                    mag_length: *mag_length,
                    members: Vec::new(),
                });
                group.members.push(MagContigJob {
                    tid: *tid,
                    contig_id: *contig_id,
                    ref_name: ref_name.clone(),
                    ref_length: *ref_length,
                    offset_in_mag: *offset,
                });
            }
        }
        drop(jobs); // Group members now own the selected job names.
        // Sort members within each group by offset so smooth_per_contig works correctly
        let mut mag_groups: Vec<MagJobGroup> = mag_groups_map.into_values().collect();
        for g in &mut mag_groups {
            g.members.sort_by_key(|m| m.offset_in_mag);
        }
        // Nearby reference IDs stay on the same worker instead of hash-map order.
        mag_groups.sort_unstable_by_key(|g| g.members.iter().map(|m| m.tid).min().unwrap_or(u32::MAX));

        let chunk_size = mag_groups.len().div_ceil(reader_workers).max(1);
        let mag_results_vec: Vec<MagGroupResult> = mag_groups
            .par_chunks(chunk_size.max(1))
            .flat_map_iter(|chunk| {
                let ts_open = accum_ref.map(|_| std::time::Instant::now());
                let mut bam = match IndexedReader::from_path(bam_path) {
                    Ok(b) => Some(b),
                    Err(e) => {
                        *stream_error_ref.lock().unwrap()=Some(anyhow::anyhow!("Opening BAM worker for {}: {}",bam_path.display(),e));
                        None
                    }
                };
                if let (Some(a), Some(t)) = (accum_ref, ts_open) {
                    a.bam_open_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }

                chunk.iter().filter_map(move |group| {
                    let bam = bam.as_mut()?;
                    if stream_error_ref.lock().unwrap().is_some() { return None; }
                    let r = match process_mag_group(bam, group, seq_type, flags, config, is_circular, repeats, cds_indexes, accum_ref) {
                        Ok(Some(result))=>result, Ok(None)=>return None,
                        Err(error)=>{ *stream_error_ref.lock().unwrap()=Some(error); return None; }
                    };
                    if let Some(sink) = mag_sink {
                        streamed_contigs_ref.fetch_add(r.presences.len(), Ordering::Relaxed);
                        // Keep only cheap summaries in the sample result. The bounded
                        // channel owns and releases the large per-contig payloads.
                        let summary = MagGroupResult {
                            mag_name: r.mag_name.clone(), mag_id: r.mag_id,
                            mag_cov_stats: r.mag_cov_stats.clone(), mapped_reads: r.mapped_reads,
                            contig_blobs: vec![], presences: vec![], packaging: vec![],
                            misassembly: vec![], microdiversity: vec![], side_misassembly: vec![], topology: vec![],
                        };
                        let queue_start = accum_ref.map(|_| std::time::Instant::now());
                        let sent = sink(MagStreamMessage::Group(r));
                        if let (Some(a), Some(start)) = (accum_ref, queue_start) {
                            a.mag_queue_wait_ns.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                        }
                        if let Err(error) = sent {
                            *stream_error_ref.lock().unwrap() = Some(error);
                            return None;
                        }
                        Some(summary)
                    } else { Some(r) }
                })
            })
            .collect();

        if let (Some(s), Some(t)) = (st.as_mut(), t_par) {
            s.par_iter_wall = t.elapsed();
            s.n_contigs_processed = if mag_sink.is_some() { streamed_contigs.load(Ordering::Relaxed) }
                else { mag_results_vec.iter().map(|r| r.presences.len()).sum() };
        }

        let t0 = if time_on { Some(std::time::Instant::now()) } else { None };
        let mut all_blobs = Vec::new();
        let mut all_pres = Vec::new();
        let mut all_pkg = Vec::new();
        let mut all_mis = Vec::new();
        let mut all_micro = Vec::new();
        let mut all_side = Vec::new();
        let mut all_topo = Vec::new();
        let mut mapped: u64 = 0;
        let mut mag_res = Vec::new();

        for r in mag_results_vec {
            all_blobs.extend(r.contig_blobs);
            all_pres.extend(r.presences);
            all_pkg.extend(r.packaging);
            all_mis.extend(r.misassembly);
            all_micro.extend(r.microdiversity);
            all_side.extend(r.side_misassembly);
            all_topo.extend(r.topology);
            mapped += r.mapped_reads;
            mag_res.push((r.mag_name, r.mag_id, r.mag_cov_stats));
        }

        if let (Some(s), Some(t)) = (st.as_mut(), t0) {
            s.merge_results = t.elapsed();
        }

        (all_blobs, all_pres, all_pkg, all_mis, all_micro, all_side, all_topo, mapped, mag_res)
    } else {
        // Contig mode: par_chunks over individual contigs
        let chunk_size = jobs.len().div_ceil(reader_workers).max(1);
        let results: Vec<_> = jobs
            .par_chunks(chunk_size.max(1))
            .flat_map_iter(|chunk| {
                let ts_open = accum_ref.map(|_| std::time::Instant::now());
                let mut bam = match IndexedReader::from_path(bam_path) {
                    Ok(b) => Some(b),
                    Err(e) => {
                        *stream_error_ref.lock().unwrap()=Some(anyhow::anyhow!("Opening BAM worker for {}: {}",bam_path.display(),e));
                        None
                    }
                };
                if let (Some(a), Some(t)) = (accum_ref, ts_open) {
                    a.bam_open_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }

                chunk.iter().filter_map(move |(tid, contig_id, ref_name, ref_length)| {
                    if stream_error_ref.lock().unwrap().is_some() {return None;}
                    let bam = bam.as_mut()?;
                    let tid = *tid;
                    let contig_id = *contig_id;
                    let ref_length = *ref_length;


                    let ts = accum_ref.map(|_| std::time::Instant::now());
                    let mut phase_t: Option<BamPhaseTimings> = accum_ref.map(|_| BamPhaseTimings::default());
                    let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming_counted(bam, tid, ref_name, ref_length, seq_type, flags, is_circular, config.min_aligned_fraction, config.phagetermini_config.min_clipping_length, &mut phase_t).map(|(result,count)| { config.primary_reads_seen.fetch_add(count,Ordering::Relaxed); result }) {
                        Ok(Some(result)) => result,
                        Ok(None) => {
                            if let (Some(a), Some(pt)) = (accum_ref, phase_t.as_ref()) {
                                a.stream_fetch_seek_ns.fetch_add(pt.fetch_seek_ns, Ordering::Relaxed);
                                a.stream_records_ns.fetch_add(pt.records_ns, Ordering::Relaxed);
                                a.stream_process_read_ns.fetch_add(pt.process_read_ns, Ordering::Relaxed);
                                a.stream_finalize_ns.fetch_add(pt.finalize_ns, Ordering::Relaxed);
                            }
                            if let (Some(a), Some(t)) = (accum_ref, ts) {
                                a.streaming_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                            }
                            if config.min_aligned_fraction > 0.0 {
                                config.contig_drops_counter.fetch_add(1, Ordering::Relaxed);
                            }
                            return None;
                        }
                        Err(e) => {
                            *stream_error_ref.lock().unwrap()=Some(anyhow::anyhow!("Processing contig {} in {}: {}",ref_name,bam_path.display(),e));
                            return None;
                        }
                    };
                    if let (Some(a), Some(pt)) = (accum_ref, phase_t.as_ref()) {
                        a.stream_fetch_seek_ns.fetch_add(pt.fetch_seek_ns, Ordering::Relaxed);
                        a.stream_records_ns.fetch_add(pt.records_ns, Ordering::Relaxed);
                        a.stream_process_read_ns.fetch_add(pt.process_read_ns, Ordering::Relaxed);
                        a.stream_finalize_ns.fetch_add(pt.finalize_ns, Ordering::Relaxed);
                    }
                    if let (Some(a), Some(t)) = (accum_ref, ts) {
                        a.streaming_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    }

                    let coverage_mean = arrays.coverage_mean() as f32;
                    if config.min_coverage_depth > 0.0
                        && (coverage_mean as f64) < config.min_coverage_depth
                    {
                        config.contig_drops_counter.fetch_add(1, Ordering::Relaxed);
                        crate::features::recycle_feature_arrays(arrays);
                        return None;
                    }

                    let ts = accum_ref.map(|_| std::time::Instant::now());
                    let mut feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)> = Vec::new();
                    let cds_index = if flags.mapping_metrics && needs_codon_annotations(&arrays, config.bar_ratio * 0.01) {
                        cds_indexes.lookup(contig_id)
                    } else { None };
                    let (packaging_info, metrics_info) = add_features_from_arrays(&mut arrays, ref_name, contig_id, config, is_circular, seq_type, flags, primary_count, repeats, cds_index.as_deref(), &mut feature_blobs, accum_ref);
                    let presence = build_presence(&arrays, contig_id, ref_name,
                        coverage_pct, primary_count, coverage_mean);

                    if let (Some(a), Some(t)) = (accum_ref, ts) {
                        a.feature_calc_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    }

                    crate::features::recycle_feature_arrays(arrays);
                    Some((feature_blobs, presence, packaging_info, metrics_info, primary_count))
                })
            })
            .collect();

        if let (Some(s), Some(t)) = (st.as_mut(), t_par) {
            s.par_iter_wall = t.elapsed();
            s.n_contigs_processed = results.len();
        }

        let t0 = if time_on { Some(std::time::Instant::now()) } else { None };
        let mut all_blobs: Vec<(String, i64, crate::blob::EncodedBlob)> = Vec::new();
        let mut all_pres = Vec::new();
        let mut all_pkg = Vec::new();
        let mut all_mis = Vec::new();
        let mut all_micro = Vec::new();
        let mut all_side = Vec::new();
        let mut all_topo = Vec::new();
        let mut mapped: u64 = 0;

        for (blobs, presence, packaging, metrics, primary_count) in results {
            all_blobs.extend(blobs);
            all_pres.push(presence);
            mapped += primary_count;
            if let Some(pkg) = packaging {
                all_pkg.push(pkg);
            }
            if let Some((mis, micro, side, topo)) = metrics {
                all_mis.push(mis);
                all_micro.push(micro);
                all_side.push(side);
                all_topo.push(topo);
            }
        }

        if let (Some(s), Some(t)) = (st.as_mut(), t0) {
            s.merge_results = t.elapsed();
        }

        (all_blobs, all_pres, all_pkg, all_mis, all_micro, all_side, all_topo, mapped, vec![])
    };

    for (_,_,blob) in &all_feature_blobs { blob.ensure_valid()?; }
    if let Some(error) = stream_error.into_inner().unwrap() { return Err(error); }
    let mapped_reads = config.primary_reads_seen.load(Ordering::Relaxed);
    let total_reads = mapped_reads + unmapped_reads;

    if let Some(s) = st.as_mut() {
        s.total_reads = total_reads;
        if let Some(a) = accum.as_ref() {
            s.bam_open_ns = a.bam_open_ns.load(Ordering::Relaxed);
            s.streaming_ns = a.streaming_ns.load(Ordering::Relaxed);
            s.feature_calc_ns = a.feature_calc_ns.load(Ordering::Relaxed);
            s.mag_queue_wait_ns = a.mag_queue_wait_ns.load(Ordering::Relaxed);
            s.sparse_encode_ns = a.sparse_encode_ns.load(Ordering::Relaxed);
            s.dense_encode_ns = a.dense_encode_ns.load(Ordering::Relaxed);
            s.stream_fetch_seek_ns = a.stream_fetch_seek_ns.load(Ordering::Relaxed);
            s.stream_records_ns = a.stream_records_ns.load(Ordering::Relaxed);
            s.stream_process_read_ns = a.stream_process_read_ns.load(Ordering::Relaxed);
            s.stream_finalize_ns = a.stream_finalize_ns.load(Ordering::Relaxed);
        }
    }

    if let Some(search) = &config.terminal_repeats {
        search.check_error()?;
        let phases = search.finish_sample(&sample_name);
        if let Some(timing) = st.as_mut() { timing.terminal = phases; }
    }
    cds_indexes.check_error()?;
    Ok(SampleProcessResult {
        feature_blobs: all_feature_blobs,
        presences: all_presences,
        packaging: all_packaging,
        misassembly: all_misassembly,
        microdiversity: all_microdiversity,
        side_misassembly: all_side_misassembly,
        topology: all_topology,
        sample_name,
        seq_type,
        total_reads,
        mapped_reads,
        bai_mapped_reads,
        is_circular,
        timings: st,
        mag_results,
    })
}

/// Extract contig information from BAM file headers.
/// Deduplicates contigs across all BAM files.
fn extract_contigs_from_bams(bam_files: &[PathBuf]) -> Result<Vec<ContigInfo>> {
    let mut contig_map: HashMap<String, usize> = HashMap::new();

    // Scan all BAM files to collect unique contigs
    for bam_path in bam_files {
        let bam = IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path.display()))?;
        let header = bam.header();

        for tid in 0..header.target_count() {
            let ref_name = std::str::from_utf8(header.tid2name(tid))
                .context("Invalid UTF-8 in reference name")?
                .to_string();
            let actual_length = header.target_len(tid).unwrap_or(0) as usize;

            // Use the contig if not seen, or verify length matches
            contig_map.entry(ref_name.clone())
                .and_modify(|existing_len| {
                    if *existing_len != actual_length {
                        eprintln!("Warning: Contig '{}' has different lengths across BAM files ({} vs {})",
                                  ref_name, *existing_len, actual_length);
                    }
                })
                .or_insert(actual_length);
        }
    }
    
    // Convert to sorted vector of ContigInfo
    let mut contigs: Vec<ContigInfo> = contig_map
        .into_iter()
        .map(|(name, length)| ContigInfo {
            name,
            aliases: vec![],
            length,
            sequence: None, // BAM headers don't contain sequence data
        })
        .collect();
    
    // Sort by name for consistent ordering
    contigs.sort_by(|a, b| a.name.cmp(&b.name));
    
    Ok(contigs)
}

/// Run processing on all BAM files.
pub fn run_all_samples(
    genbank_path: &Path,
    assembly_path: &Path,
    bam_files: &[PathBuf],
    output_db: &Path,
    modules: &[String],
    config: &ProcessConfig,
    create_indexes: bool,
    extend_db: &Path,
    command_line: &str,
) -> Result<ProcessResult> {
    anyhow::ensure!(config.blast != BlastMode::Mags || config.view_mode == ViewMode::Mag,
        "--blast mags requires --view mag");
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let t_run_start = std::time::Instant::now();

    // 1. Strict upfront validation — fail early with actionable errors
    if !bam_files.is_empty() {
        validate_inputs(bam_files, modules)?;
    }

    // Open timing log as early as possible so all phases are captured
    let mut timing_file: Option<fs::File> = if config.enable_timing {
        match open_timing_log(output_db, config.threads, bam_files.len(), command_line) {
            Ok(f) => Some(f),
            Err(e) => { eprintln!("Warning: could not open timing log: {}", e); None }
        }
    } else {
        None
    };

    // 2. Parse annotations (GenBank/GFF3) or extract contigs from BAMs
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Parsing input files...");
    eprintln!();
    let t_parse = if config.enable_timing { Some(std::time::Instant::now()) } else { None };

    let has_genbank = !genbank_path.as_os_str().is_empty();
    let is_extending = !extend_db.as_os_str().is_empty();

    // In MAG mode, per-MAG membership (mag_name → member contig names) is collected during parsing.
    let mut mag_contig_map: Vec<(String, Vec<String>)> = Vec::new();

    let mut streamed_writer = None;
    let mut sequence_store = None;
    let mut streamed_annotation_count = None;
    let streaming = !is_extending && crate::annotation_import::can_stream(genbank_path, assembly_path, config);
    let (mut contigs, mut annotations, mut contig_qualifiers) = if streaming {
        eprintln!("Streaming annotations one file at a time; temporary sequences stored beside the output database");
        let imported = crate::annotation_import::import_annotations(
            genbank_path, output_db, config, !bam_files.is_empty(),
            false,
            |done, total, count| {
                let rss = get_rss_mb();
                eprintln!("Annotation import: {done}/{total} files, {count} annotations [RSS: {rss:.0} MiB]");
                if let Some(ref mut f) = timing_file {
                    use std::io::Write;
                    let _ = writeln!(f, "  Annotation import: {done}/{total} files, {count} annotations [RSS: {rss:.0} MiB]");
                    let _ = f.flush();
                }
            },
        )?;
        streamed_writer = Some(imported.writer);
        sequence_store = Some(imported.sequences);
        mag_contig_map = imported.mags;
        streamed_annotation_count = Some(imported.annotation_count);
        (imported.contigs, imported.cds, Vec::new())
    } else if config.view_mode == ViewMode::Mag {
        // Parse each MAG's annotation and/or FASTA file individually. Contig names have
        // already been validated as globally unique by the Python side.
        let mut all_contigs: Vec<ContigInfo> = Vec::new();
        let mut all_annotations: Vec<FeatureAnnotation> = Vec::new();
        let mut all_quals: Vec<(i64, HashMap<String, String>)> = Vec::new();

        for mag in &config.mag_manifest {
            let mut member_names: Vec<String> = Vec::new();

            if let Some(gb) = &mag.gb_path {
                let (mag_contigs, mag_anns, mag_quals) = parse_annotations(gb)?;
                // Remap parser-local contig_ids to their global position in all_contigs.
                let base_global = all_contigs.len() as i64;
                for c in &mag_contigs {
                    member_names.push(c.name.clone());
                }
                for mut ann in mag_anns {
                    ann.contig_id += base_global;
                    all_annotations.push(ann);
                }
                for (qid, qmap) in mag_quals {
                    all_quals.push((qid + base_global, qmap));
                }
                all_contigs.extend(mag_contigs);
            }

            if let Some(asm) = &mag.asm_path {
                // Merge FASTA sequences into the contigs that came from gb, or create new
                // contigs if there was no gb for this MAG.
                if mag.gb_path.is_some() {
                    // member_names already populated from annotation; merge sequences.
                    let fasta_records = crate::parser::parse_fasta(asm)?;
                    let name_to_idx: HashMap<String, usize> = all_contigs.iter()
                        .enumerate()
                        .map(|(i, c)| (c.name.clone(), i))
                        .collect();
                    let mut unmatched_count = 0usize;
                    for (name, seq) in fasta_records {
                        if let Some(&idx) = name_to_idx.get(&name) {
                            if all_contigs[idx].sequence.is_none() {
                                all_contigs[idx].sequence = Some(seq);
                            }
                        } else {
                            // FASTA contig name doesn't match any annotation contig (e.g.
                            // annotation used the MAG bin name while the assembly uses
                            // assembler node names). Add as a real contig and include in
                            // this MAG's member list so BAM coverage is attributed to the MAG.
                            member_names.push(name.clone());
                            all_contigs.push(ContigInfo {
                                name,
                                aliases: vec![],
                                length: seq.len(),
                                sequence: Some(seq),
                            });
                            unmatched_count += 1;
                        }
                    }
                    if unmatched_count > 0 {
                        eprintln!(
                            "WARNING: MAG '{}': {} assembly FASTA contigs not found in annotation \
                             file — added as contigs for BAM matching. Annotation seqids and \
                             assembly contig names differ (check your annotation pipeline).",
                            mag.name, unmatched_count
                        );
                    }
                } else {
                    let fasta_records = crate::parser::parse_fasta(asm)?;
                    for (name, seq) in fasta_records {
                        member_names.push(name.clone());
                        all_contigs.push(ContigInfo {
                            name,
                            aliases: vec![],
                            length: seq.len(),
                            sequence: Some(seq),
                        });
                    }
                }
            }

            mag_contig_map.push((mag.name.clone(), member_names));
        }

        // Fallback: single annotation file covers all MAGs (e.g. -g one_big.gff -a mags_dir/).
        if all_annotations.is_empty() && has_genbank && genbank_path.is_file() {
            eprintln!("No per-MAG annotation files; parsing global annotation file for all MAGs...");
            let (parser_contigs, global_anns, global_quals) = parse_annotations(genbank_path)?;

            let local_to_name: HashMap<i64, &str> = parser_contigs
                .iter()
                .enumerate()
                .map(|(i, c)| ((i + 1) as i64, c.name.as_str()))
                .collect();

            let name_to_global: HashMap<&str, i64> = all_contigs
                .iter()
                .enumerate()
                .map(|(i, c)| (c.name.as_str(), (i + 1) as i64))
                .collect();

            for mut ann in global_anns {
                if let Some(name) = local_to_name.get(&ann.contig_id) {
                    if let Some(&global_id) = name_to_global.get(name) {
                        ann.contig_id = global_id;
                        all_annotations.push(ann);
                    }
                }
            }
            for (qid, qmap) in global_quals {
                if let Some(name) = local_to_name.get(&qid) {
                    if let Some(&global_id) = name_to_global.get(name) {
                        all_quals.push((global_id, qmap));
                    }
                }
            }
            eprintln!("Distributed {} annotations across {} contigs from global annotation file",
                all_annotations.len(), all_contigs.len());
        }

        (all_contigs, all_annotations, all_quals)
    } else if !has_genbank {
        eprintln!("No GenBank file provided - extracting contigs from BAM headers");
        let contigs = extract_contigs_from_bams(bam_files)?;
        eprintln!("Found {} contigs from BAM files", contigs.len());
        (contigs, Vec::new(), Vec::new())
    } else if has_genbank && genbank_path.is_dir() {
        // Directory of annotation files in contig mode — parse each file and merge.
        let mut all_contigs: Vec<ContigInfo> = Vec::new();
        let mut all_annotations: Vec<FeatureAnnotation> = Vec::new();
        let mut all_quals: Vec<(i64, HashMap<String, String>)> = Vec::new();

        let mut entries: Vec<_> = std::fs::read_dir(genbank_path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.is_file() && crate::parser::is_annotation_file(p))
            .collect();
        entries.sort();

        if entries.is_empty() {
            return Err(anyhow::anyhow!(
                "No annotation files found in directory: {}",
                genbank_path.display()
            ));
        }

        for path in &entries {
            let base_global = all_contigs.len() as i64;
            let (contigs, anns, quals) = parse_annotations(path)?;
            for mut ann in anns {
                ann.contig_id += base_global;
                all_annotations.push(ann);
            }
            for (qid, qmap) in quals {
                all_quals.push((qid + base_global, qmap));
            }
            all_contigs.extend(contigs);
        }
        (all_contigs, all_annotations, all_quals)
    } else {
        parse_annotations(genbank_path)?
    };
    eprintln!(
        "Found {} contigs with {} annotations",
        contigs.len(),
        streamed_annotation_count.unwrap_or(annotations.len())
    );

    eprintln!("Found {} BAM files", bam_files.len());

    // 3. If assembly_path provided, parse FASTA and merge sequences into contigs.
    //    Skipped in MAG mode — sequences are already merged per-MAG above.
    if config.view_mode != ViewMode::Mag && !assembly_path.as_os_str().is_empty() {
        if assembly_path.is_file() {
            eprintln!();
            eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            eprintln!("### Merging sequences from assembly FASTA...");
            merge_sequences_from_fasta(&mut contigs, assembly_path)?;
        } else if assembly_path.is_dir() {
            eprintln!();
            eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            eprintln!("### Merging sequences from assembly FASTA directory...");
            let mut fasta_entries: Vec<_> = std::fs::read_dir(assembly_path)?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| p.is_file() && crate::parser::is_fasta_file(p))
                .collect();
            fasta_entries.sort();
            for path in &fasta_entries {
                merge_sequences_from_fasta(&mut contigs, path)?;
            }
        }
    }

    // 3b. Compute nucleotide/protein sequences for CDS annotations
    if !streaming && !annotations.is_empty() {
        compute_annotation_sequences(&contigs, &mut annotations);
    }
    let parse_secs = t_parse.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let parse_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f);
        let _ = writeln!(f, "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        let _ = writeln!(f, "### Pre-sample phases");
        let _ = writeln!(f);
        let _ = writeln!(f, "  Annotation+FASTA parse: {:>9.3} s  [RSS: {:.0} MB]", parse_secs, parse_rss_mb);
        let _ = f.flush();
    }

    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent).context("Failed to create output directory")?;
    }

    std::thread::sleep(std::time::Duration::from_millis(50));

    // 4. Create or open database
    let t_db_create = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    // For extend mode: read existing contigs, determine new contigs, open existing DB
    let (mut db_writer, mut new_contigs_for_blast, repeats_from_db) = if is_extending {
        // Read existing contig names from the database
        let existing_conn = duckdb::Connection::open(extend_db)
            .with_context(|| format!("Failed to open existing database: {}", extend_db.display()))?;
        crate::db::configure_for_bulk_writes(&existing_conn)?;

        // View_mode compatibility check: reject cross-mode extends early.
        let existing_view_mode: String = existing_conn
            .query_row(
                "SELECT Value FROM Database_metadata WHERE Key = 'View_mode'",
                [],
                |row| row.get::<_, String>(0),
            )
            .unwrap_or_else(|_| "contig".to_string());
        let existing_is_mag = existing_view_mode == "mag";
        let new_is_mag = config.view_mode == ViewMode::Mag;
        if existing_is_mag && !new_is_mag {
            return Err(anyhow::anyhow!(
                "existing database was built in MAG mode; --extend requires --view mag"
            ));
        }
        if !existing_is_mag && new_is_mag {
            return Err(anyhow::anyhow!(
                "cannot switch a contig-mode database to MAG mode; rebuild from scratch"
            ));
        }

        let mut existing_contig_names: HashSet<String> = HashSet::new();
        let mut existing_contigs: Vec<ContigInfo> = Vec::new();
        {
            let mut stmt = existing_conn.prepare(
                "SELECT Contig_name, Contig_length FROM Contig"
            )?;
            let rows = stmt.query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?))
            })?;
            for row in rows {
                let (name, length) = row?;
                existing_contig_names.insert(name.clone());
                existing_contigs.push(ContigInfo {
                    name,
                    aliases: vec![],
                    length: length as usize,
                    sequence: None,
                });
            }
        }

        // Load existing MAG names (for MAG mode extend validation)
        let existing_mag_names: HashSet<String> = {
            let has_mag_table: bool = existing_conn.query_row(
                "SELECT COUNT(*) > 0 FROM information_schema.tables WHERE table_name = 'MAG'",
                [], |row| row.get(0),
            ).unwrap_or(false);
            if has_mag_table {
                let mut names = HashSet::new();
                let mut stmt = existing_conn.prepare("SELECT MAG_name FROM MAG")?;
                let rows = stmt.query_map([], |row| row.get::<_, String>(0))?;
                for row in rows {
                    names.insert(row?);
                }
                names
            } else {
                HashSet::new()
            }
        };

        // Detect if DB was created with annotations (genbank mode)
        let db_has_annotations: bool = existing_conn
            .query_row("SELECT COUNT(*) FROM Contig_annotation", [], |row| row.get::<_, i64>(0))
            .unwrap_or(0) > 0;

        // Global search recomputes reference evidence without retaining old repeat DNA metadata.
        let mut repeats_from_db: Vec<crate::db::RepeatsData> = Vec::new();
        if !(config.blast==BlastMode::All) {
            // Build Contig_id → name mapping for repeat loading
            let mut id_to_name: HashMap<i64, String> = HashMap::new();
            {
                let mut stmt = existing_conn.prepare("SELECT Contig_id, Contig_name FROM Contig")?;
                let rows = stmt.query_map([], |row| {
                    Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
                })?;
                for row in rows {
                    let (id, name) = row?;
                    id_to_name.insert(id, name);
                }
            }
            for (table, is_direct) in &[("Contig_directRepeats", true), ("Contig_invertedRepeats", false)] {
                let has_table: bool = existing_conn.query_row(
                    "SELECT COUNT(*) > 0 FROM information_schema.tables WHERE table_name = ?",
                    duckdb::params![table], |row| row.get(0),
                ).unwrap_or(false);
                if !has_table { continue; }
                let mut stmt = existing_conn.prepare(&format!(
                    "SELECT Contig_id, Position1, Position2, Position1prime, Position2prime, Pident FROM {}", table
                ))?;
                let rows = stmt.query_map([], |row| {
                    Ok((
                        row.get::<_, i64>(0)?,
                        row.get::<_, i32>(1)?,
                        row.get::<_, i32>(2)?,
                        row.get::<_, i32>(3)?,
                        row.get::<_, i32>(4)?,
                        row.get::<_, i32>(5)?,
                    ))
                })?;
                for row in rows {
                    let (cid, p1, p2, p1p, p2p, pident_int) = row?;
                    if let Some(name) = id_to_name.get(&cid) {
                        repeats_from_db.push(crate::db::RepeatsData {
                            contig_name: name.clone(),
                            position1: p1,
                            position2: p2,
                            position1prime: p1p,
                            position2prime: p2p,
                            pident: pident_int as f64 / crate::types::get_column_scale("Contig_directRepeats", "Pident"),
                            is_direct: *is_direct,
                        });
                    }
                }
            }
        }

        drop(existing_conn);

        // Determine new contigs based on DB origin mode
        let new_contigs: Vec<ContigInfo>;
        let mut new_annotations: Vec<FeatureAnnotation>;
        let new_contig_qualifiers: Vec<(i64, HashMap<String, String>)>;

        if has_genbank {
            // GenBank provided: filter to only NEW contig names
            // Build parser contig_id → contig_name mapping (parser uses 1-based IDs)
            let parser_id_to_name: HashMap<i64, String> = contigs.iter()
                .enumerate()
                .map(|(i, c)| ((i + 1) as i64, c.name.clone()))
                .collect();
            let new_contig_names: HashSet<String> = contigs.iter()
                .filter(|c| !existing_contig_names.contains(&c.name))
                .map(|c| c.name.clone())
                .collect();
            new_contigs = contigs.into_iter()
                .filter(|c| new_contig_names.contains(&c.name))
                .collect();
            // Filter annotations to only new contigs, and remap contig_ids
            // DbWriter::open will assign new IDs starting from max_contig_id + 1
            // Build a mapping from new contig name → position in new_contigs vec
            let new_contig_name_to_pos: HashMap<String, usize> = new_contigs.iter()
                .enumerate()
                .map(|(i, c)| (c.name.clone(), i))
                .collect();
            new_annotations = annotations.into_iter()
                .filter(|a| {
                    parser_id_to_name.get(&a.contig_id)
                        .map(|name| new_contig_names.contains(name))
                        .unwrap_or(false)
                })
                .collect();
            // Remap contig_ids: use position+1 in new_contigs vec as temp IDs
            // DbWriter::open will add max_contig_id offset
            for ann in &mut new_annotations {
                if let Some(name) = parser_id_to_name.get(&ann.contig_id) {
                    if let Some(&pos) = new_contig_name_to_pos.get(name) {
                        ann.contig_id = (pos + 1) as i64;
                    }
                }
            }
            // Filter and remap contig qualifiers to new contigs only
            new_contig_qualifiers = std::mem::take(&mut contig_qualifiers).into_iter()
                .filter_map(|(cid, quals)| {
                    parser_id_to_name.get(&cid)
                        .and_then(|name| new_contig_name_to_pos.get(name))
                        .map(|&pos| ((pos + 1) as i64, quals))
                })
                .collect();
            if !new_contigs.is_empty() {
                eprintln!("Found {} new contigs to add (filtered from {} total)", new_contigs.len(), new_contigs.len() + existing_contig_names.len());
            }
        } else if !db_has_annotations {
            // BAM-only mode, DB was also BAM-only: auto-add new contigs from BAM headers
            new_contigs = contigs.into_iter()
                .filter(|c| !existing_contig_names.contains(&c.name))
                .collect();
            new_annotations = Vec::new();
            new_contig_qualifiers = Vec::new();
            if !new_contigs.is_empty() {
                eprintln!("Found {} new contigs from BAM headers", new_contigs.len());
            }
        } else {
            // BAM-only extend on genbank DB: no new contigs, process against existing only
            new_contigs = Vec::new();
            new_annotations = Vec::new();
            new_contig_qualifiers = Vec::new();
            eprintln!("No new contigs (genbank-mode database, no -g/-a provided)");
        }

        // Build all_contigs for BAM processing: existing + new
        // Replace contigs vec with the combined list
        contigs = existing_contigs;
        for nc in &new_contigs {
            contigs.push(nc.clone());
        }
        annotations = Vec::new();

        // In MAG mode, validate that new contigs only belong to new MAGs.
        // Adding contigs to an existing MAG would invalidate its aggregates.
        if config.view_mode == ViewMode::Mag && !mag_contig_map.is_empty() {
            let new_contig_set: HashSet<&str> = new_contigs.iter().map(|c| c.name.as_str()).collect();
            for (mag_name, members) in &mag_contig_map {
                if existing_mag_names.contains(mag_name) {
                    let new_in_existing: Vec<&str> = members.iter()
                        .filter(|m| new_contig_set.contains(m.as_str()))
                        .map(|m| m.as_str())
                        .collect();
                    if !new_in_existing.is_empty() {
                        return Err(anyhow::anyhow!(
                            "Cannot add contigs to existing MAG '{}': {:?}. \
                             Adding contigs to an existing MAG would invalidate its aggregates. \
                             Rebuild from scratch instead.",
                            mag_name, new_in_existing
                        ));
                    }
                }
            }
            // Filter mag_contig_map to only new MAGs
            mag_contig_map.retain(|(name, _)| !existing_mag_names.contains(name));
        }

        let writer = DbWriter::open(extend_db, &new_contigs, &new_annotations, &new_contig_qualifiers, !bam_files.is_empty())?;
        if config.enable_timing { writer.enable_size_report(extend_db)?; writer.record_import_sizes(&new_contigs, &new_annotations, &new_contig_qualifiers)?; }
        writer.update_metadata_modification()?;
        writer.write_compression_metadata(
            config.min_aligned_fraction,
            config.min_coverage_depth,
            config.bar_ratio,
            config.variation_percentage,
            config.min_occurrences,
        )?;
        (writer, Some(new_contigs), repeats_from_db)
    } else {
        let writer = if let Some(writer) = streamed_writer.take() {
            writer
        } else { DbWriter::create(
            output_db,
            &contigs,
            &annotations,
            &contig_qualifiers,
            !bam_files.is_empty(),
            config.view_mode == ViewMode::Mag,
        )? };
        if config.enable_timing && writer.enable_size_report(output_db)? { writer.record_import_sizes(&contigs, &annotations, &contig_qualifiers)?; }
        writer.write_metadata(
            modules,
            config.min_aligned_fraction,
            config.min_coverage_depth,
            config.bar_ratio,
            config.variation_percentage,
            config.min_occurrences,
            config.view_mode.as_str(),
            config.blast.reference_search(),
        )?;
        (writer, None, Vec::new())
    };

    let persisted_annotation_count = streamed_annotation_count.unwrap_or(annotations.len());
    annotations = Vec::new();
    drop(contig_qualifiers);
    db_writer.prepare_terminase_cache(ModuleFlags::from_modules(modules).phagetermini)?;
    db_writer.write_blast_mode(config.blast.as_str())?;
    let db_create_secs = t_db_create.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let db_create_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  DB create/open       : {:>10.3} s  [RSS: {:.0} MB]", db_create_secs, db_create_rss_mb);
        let _ = f.flush();
    }

    // 5. Auto-detect circularity per sample (now with contigs available for length comparison)
    let t_circ = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let circularity_map = if !bam_files.is_empty() {
        detect_all_sample_circularities(bam_files)?
    } else {
        HashMap::new()
    };
    let circ_secs = t_circ.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let circ_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  Circularity detect   : {:>10.3} s  [RSS: {:.0} MB]", circ_secs, circ_rss_mb);
        for bam_path in bam_files {
            let sample = bam_path.file_stem().unwrap_or_default().to_string_lossy().replace("_with_MD", "");
            let _ = writeln!(f, "  Sample '{}': circular={}", sample, circularity_map.get(bam_path).copied().unwrap_or(false));
        }
        let _ = f.flush();
    }


    // 6. If --blast enabled and sequences available, run repeat detection.
    // Contig mode: per-contig BLAST self-alignment (autoblast).
    // MAG mode: per-MAG BLAST handles both intra-contig repeats and inter-contig hits.
    // Non-streaming imports also release assembly DNA once it is persisted.
    // GC/general BLAST consume disk-backed batches; termini use database candidates.
    if sequence_store.is_none() {
        let directory = output_db.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let mut store = SequenceStore::new(directory)?;
        match new_contigs_for_blast.as_mut() {
            Some(new_contigs) => {
                store.push_batch(new_contigs)?;
                for contig in &mut contigs { contig.sequence = None; }
            }
            None => store.push_batch(&mut contigs)?,
        }
        sequence_store = Some(store);
    }
    let blast_contigs = new_contigs_for_blast.as_deref().unwrap_or(&contigs);
    let has_sequence = |i: usize, c: &ContigInfo| sequence_store.as_ref()
        .map_or(c.sequence.is_some(), |store| store.has_sequence(i));
    let n_contigs_with_seq = blast_contigs.iter().enumerate().filter(|(i, c)| has_sequence(*i, c)).count();
    let has_sequences = n_contigs_with_seq > 0;
    let mut general_assessed_ids = if config.blast.reference_search() {
        blast_contigs.iter().enumerate().filter(|(i,c)| has_sequence(*i,c))
            .filter_map(|(_,c)| db_writer.get_contig_id(&c.name)).collect::<Vec<_>>()
    } else { Vec::new() };

    let t_autoblast = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let global_collection=config.blast==BlastMode::All;
    let repeats = if global_collection {
        let (assessed,hit_count)=reference_blast::run_global(&db_writer,&contigs)?;
        general_assessed_ids=assessed;
        eprintln!("Reference BLAST all: {hit_count} inter-contig hits read (including matches between contigs)");
        Vec::new()
    } else if config.blast.reference_search() && has_sequences {
        if config.view_mode == ViewMode::Mag {
            let blast_names: HashSet<&str> = blast_contigs.iter().map(|c| c.name.as_str()).collect();
            let filtered_manifest: Vec<(String, Vec<String>)> = mag_contig_map
                .iter()
                .map(|(mn, members)| {
                    let kept: Vec<String> = members.iter()
                        .filter(|n| blast_names.contains(n.as_str()))
                        .cloned()
                        .collect();
                    (mn.clone(), kept)
                })
                .filter(|(_, m)| !m.is_empty())
                .collect();
            let (reps, inter_count) = run_mag_interblast(blast_contigs, &filtered_manifest, &db_writer, sequence_store.as_ref())?;
            db_writer.compute_duplication_percentages(&reps)?;
            eprintln!("Reference MAG BLAST: {inter_count} inter-contig hits read");
            reps
        } else {
            let reps = run_autoblast(blast_contigs, sequence_store.as_ref())?;
            db_writer.write_repeats(&reps)?;
            db_writer.compute_duplication_percentages(&reps)?;
            reps
        }
    } else {
        if !config.blast.reference_search() && has_sequences {
            eprintln!();
            eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            eprintln!("### BLAST repeat detection: skipped (use --blast all to enable)");
        }
        Vec::new()
    };
    let autoblast_secs = t_autoblast.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let autoblast_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  Contigs with sequence: {}", n_contigs_with_seq);
        let _ = writeln!(f, "  Autoblast            : {:>10.3} s  [RSS: {:.0} MB]", autoblast_secs, autoblast_rss_mb);
        let _ = f.flush();
    }

    let repeats = if repeats.is_empty() && !repeats_from_db.is_empty() {
        repeats_from_db
    } else {
        repeats
    };

    if config.blast.reference_search() && has_sequences {
        let compatible: bool = db_writer.lock_conn()?.query_row(
            "SELECT Value=? FROM Database_metadata WHERE Key='Blast_parameters'",
            duckdb::params![BLAST_PARAMETERS], |r| r.get(0)).unwrap_or(false);
        if !is_extending || compatible { db_writer.write_blast_parameters()?; }
    }

    // 7. Compute and write GC content and GC skew from sequence data (only for new contigs)
    // Process in batches of 10k to bound memory and replace millions of individual
    // UPDATE statements with one bulk UPDATE FROM per batch.
    let t_gc = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let contigs_with_seq: Vec<usize> = blast_contigs.iter().enumerate()
        .filter(|(i, c)| has_sequence(*i, c)).map(|(i, _)| i).collect();
    let n_gc_contigs = contigs_with_seq.len();
    if n_gc_contigs > 0 {
        eprintln!();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("### Computing GC content for {} contigs...", n_gc_contigs);
    }
    let mut gc_compute_ns: u64 = 0;
    let mut gc_write_blob_ns: u64 = 0;
    let mut gc_update_stats_ns: u64 = 0;
    const GC_BATCH: usize = 10_000;
    let mut batch_start = 0;
    while batch_start < n_gc_contigs {
        // Bound both record count and bases; one exceptionally large contig is
        // processed alone. GC output arrays are substantially larger than DNA.
        let mut batch_end = batch_start;
        let mut bases = 0usize;
        while batch_end < n_gc_contigs && batch_end - batch_start < GC_BATCH {
            let len = blast_contigs[contigs_with_seq[batch_end]].length;
            if batch_end > batch_start && bases.saturating_add(len) > 8 * 1024 * 1024 { break; }
            bases = bases.saturating_add(len);
            batch_end += 1;
        }
        let batch = &contigs_with_seq[batch_start..batch_end];
        let loaded: Vec<std::borrow::Cow<'_, ContigInfo>> = batch.iter().map(|&i| {
            match sequence_store.as_ref() {
                Some(store) => store.load(i, &blast_contigs[i]).map(std::borrow::Cow::Owned),
                None => Ok(std::borrow::Cow::Borrowed(&blast_contigs[i])),
            }
        }).collect::<Result<_>>()?;

        let t0 = std::time::Instant::now();
        let gc_data: Vec<GCContentData> = loaded
            .par_iter()
            .filter_map(|contig| {
                contig.sequence.as_ref().map(|seq| {
                    let (gc_values, stats, skew_values, skew_stats) = compute_gc_content_and_skew(
                        seq,
                        config.gc_params.gc_content_window_size,
                        config.gc_params.gc_skew_window_size,
                    );
                    GCContentData {
                        contig_name: contig.name.clone(),
                        gc_values,
                        skew_values,
                        stats,
                        skew_stats,
                    }
                })
            })
            .collect();
        gc_compute_ns += t0.elapsed().as_nanos() as u64;

        if !gc_data.is_empty() {
            let t1 = std::time::Instant::now();
            db_writer.write_gc_content(&gc_data)?;
            db_writer.write_gc_skew(&gc_data)?;
            gc_write_blob_ns += t1.elapsed().as_nanos() as u64;

            let t2 = std::time::Instant::now();
            db_writer.stage_contig_gc_stats(&gc_data)?;
            gc_update_stats_ns += t2.elapsed().as_nanos() as u64;


        }
        batch_start = batch_end;
    }
    if n_gc_contigs > 0 {
        let start = std::time::Instant::now();
        db_writer.flush_contig_gc_stats()?;
        gc_update_stats_ns += start.elapsed().as_nanos() as u64;
    }
    // Termini candidates are discovered later, while reading samples. Keep
    // disk-backed DNA until that stage finishes when the module is enabled.
    if !ModuleFlags::from_modules(modules).phagetermini {
        drop(sequence_store.take());
    }
    let gc_secs: [f64; 4] = [
        t_gc.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0),
        gc_compute_ns as f64 / 1e9,
        gc_write_blob_ns as f64 / 1e9,
        gc_update_stats_ns as f64 / 1e9,
    ];
    let gc_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  GC content + skew    : {:>10.3} s  [RSS: {:.0} MB]", gc_secs[0], gc_rss_mb);
        let _ = writeln!(f, "    compute            : {:>10.3} s", gc_secs[1]);
        let _ = writeln!(f, "    write blobs        : {:>10.3} s", gc_secs[2]);
        let _ = writeln!(f, "    update stats       : {:>10.3} s", gc_secs[3]);
        let _ = f.flush();
    }

    {
        let rss_before = if config.enable_timing { get_rss_mb() } else { 0.0 };
        let t_ckpt = std::time::Instant::now();
        if let Err(e) = db_writer.checkpoint() {
            eprintln!("WARNING: checkpoint after GC failed: {}", e);
        }
        if config.enable_timing {
            let rss_after = get_rss_mb();
            let ckpt_secs = t_ckpt.elapsed().as_secs_f64();
            if let Some(ref mut f) = timing_file {
                use std::io::Write;
                let _ = writeln!(f, "  checkpoint (post-GC) : {:>10.3} s  [RSS: {:.0} → {:.0} MB]",
                    ckpt_secs, rss_before, rss_after);
            }
        }
    }

    // Convert repeat detections into per-contig blobs (only for new contigs).
    let t_repeat_blob = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let new_contig_ids: Vec<i64> = blast_contigs.iter()
        .filter_map(|c| db_writer.get_contig_id(&c.name))
        .collect();
    db_writer.convert_repeat_blobs(if global_collection {&general_assessed_ids} else {&new_contig_ids})?;
    let repeat_blob_secs = t_repeat_blob.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let repeat_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  Repeat blob conversion: {:>9.3} s  [RSS: {:.0} MB]", repeat_blob_secs, repeat_rss_mb);
        let _ = f.flush();
    }

    {
        let rss_before = if config.enable_timing { get_rss_mb() } else { 0.0 };
        let t_ckpt = std::time::Instant::now();
        if let Err(e) = db_writer.checkpoint() {
            eprintln!("WARNING: checkpoint after repeat blobs failed: {}", e);
        }
        if config.enable_timing {
            let rss_after = get_rss_mb();
            let ckpt_secs = t_ckpt.elapsed().as_secs_f64();
            if let Some(ref mut f) = timing_file {
                use std::io::Write;
                let _ = writeln!(f, "  checkpoint (post-rep): {:>10.3} s  [RSS: {:.0} → {:.0} MB]",
                    ckpt_secs, rss_before, rss_after);
            }
        }
    }

    // Write MAG rows after per-contig GC/duplication stats are finalized.
    let t_write_mags = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    if config.view_mode == ViewMode::Mag && !mag_contig_map.is_empty() {
        eprintln!();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("### Writing {} MAGs ({} associations)...",
            mag_contig_map.len(),
            mag_contig_map.iter().map(|(_, v)| v.len()).sum::<usize>(),
        );
        db_writer.write_mags(&mag_contig_map)?;
    }
    // Membership must exist before deriving explicitly within-MAG tracks.
    if config.blast.reference_search()
        && (global_collection || has_sequences) { db_writer.write_intercontig_hit_features()?; }
    let write_mags_secs = t_write_mags.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
    let write_mags_rss_mb = if config.enable_timing { get_rss_mb() } else { 0.0 };
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  Write MAGs           : {:>10.3} s  [RSS: {:.0} MB]", write_mags_secs, write_mags_rss_mb);
        let _ = f.flush();
    }

    // If no BAM files provided, skip sample processing (genbank-only mode)
    if bam_files.is_empty() {
        eprintln!();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("### No BAM files provided - skipping sample processing");
        eprintln!();
        eprintln!("Database populated with {} contigs and {} annotations", contigs.len(), persisted_annotation_count);

        let index_secs = if create_indexes {
            db_writer.create_serve_indexes()?
        } else {
            0.0
        };

        let t_finalize = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
        let ft = db_writer.finalize()?;
        if let Some(t) = t_finalize {
            if let Some(ref mut f) = timing_file {
                use std::io::Write;
                let _ = writeln!(f);
                let _ = writeln!(f, "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
                let _ = writeln!(f, "### Post-sample phases");
                let _ = writeln!(f);
                let _ = writeln!(f, "  Build indexes        : {:>10.3} s", index_secs);
                let _ = writeln!(f, "  DB finalize          : {:>10.3} s", t.elapsed().as_secs_f64());
                let _ = writeln!(f, "    create_views       : {:>10.3} s", ft.create_views_secs);
                let _ = writeln!(f, "    cleanup_vars       : {:>10.3} s", ft.cleanup_vars_secs);
                let _ = writeln!(f, "    drop_empty_tables  : {:>10.3} s", ft.drop_empty_secs);
                let _ = writeln!(f, "    update_counts      : {:>10.3} s", ft.update_counts_secs);
                let _ = writeln!(f, "    CHECKPOINT         : {:>10.3} s", ft.checkpoint_secs);
                let _ = f.flush();
            }
        }

        return Ok(ProcessResult {
            samples_processed: 0,
            samples_failed: 0,
            samples_no_coverage: 0,
            total_time_secs: 0.0,
            processing_time_secs: 0.0,
            writing_time_secs: 0.0,
        });
    }

    let processing_database = if is_extending { extend_db } else { output_db };

    // 8. Process samples
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!();
    eprintln!("Modules: {}", modules.join(", "));

    let flags = ModuleFlags::from_modules(modules);
    let mut runtime_config = config.clone();
    if flags.phagetermini {
        drop(sequence_store.take());
        let ids = contigs.iter().map(|c| db_writer.get_contig_id(&c.name)
            .with_context(|| format!("Missing database ID for '{}'", c.name))).collect::<Result<Vec<_>>>()?;
        let directory = output_db.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let store = SequenceStore::from_connection(directory, db_writer.lock_conn()?.try_clone()?, ids.clone())?;
        runtime_config.terminal_repeats = Some(Arc::new(crate::terminal_repeats::TerminalRepeatSearch::new(
            store, ids, config.phagetermini_config,
        )));
        eprintln!("Phage termini: BLAST mode {}", config.blast.as_str());
    }
    db_writer.write_blast_mode(config.blast.as_str())?;
    if let Some(search) = &runtime_config.terminal_repeats {
        search.configure_mode(config.blast, &db_writer)?;
        if config.blast.reference_search() && !general_assessed_ids.is_empty() {
            search.mark_general_assessed(&general_assessed_ids);
        }
    }
    let pre_timings = PreSampleTimings {
        parse_secs,
        parse_rss_mb,
        db_create_secs,
        db_create_rss_mb,
        circ_secs,
        circ_rss_mb,
        autoblast_secs,
        autoblast_rss_mb,
        gc_secs,
        gc_rss_mb,
        repeat_blob_secs,
        repeat_rss_mb,
        write_mags_secs,
        write_mags_rss_mb,
        cds_build_secs: 0.0, // filled in process_samples_parallel/sequential
        cds_build_rss_mb: 0.0,
    };
    drop(new_contigs_for_blast);
    for c in &mut contigs { c.sequence = None; }
    let result = process_samples_parallel(&bam_files, &contigs, flags, &runtime_config, &circularity_map, db_writer, &repeats, annotations, processing_database, timing_file, pre_timings, &mag_contig_map, t_run_start, create_indexes)?;

    if let Some(search) = &runtime_config.terminal_repeats { search.check_error()?; }
    print_summary(&result, output_db);

    Ok(result)
}

/// Bound BAM readers by requested threads and available memory. The reference
/// estimate protects against RSS deltas undercounting reused header allocations.
fn bam_reader_worker_limit(threads: usize, references: u64, measured_mib: u64, total_mib: u64, rss_mib: u64) -> usize {
    let reader_mib = measured_mib.max((references.saturating_mul(512)).div_ceil(1024 * 1024)).max(64);
    let reserve = (total_mib / 5).max(256);
    let available = total_mib.saturating_sub(rss_mib).saturating_sub(reserve) / 2;
    let affordable = (available / reader_mib.saturating_add(256)).max(1) as usize;
    threads.max(1).min(affordable)
}

fn mag_coverage_fails(af: f64, mean: f64, min_af: f64, min_mean: f64) -> bool {
    (min_af > 0.0 && af < min_af) || (min_mean > 0.0 && mean < min_mean)
}

/// Holds processed sample data ready for database writing.
struct SampleResult {
    sample_name: String,
    sequencing_type: SequencingType,
    total_reads: u64,
    mapped_reads: u64,
    bai_mapped_reads: u64,
    is_circular: bool,
    /// Compressed BLOB data: Vec<(feature_name, contig_name, encoded_blob)>
    feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    timings: Option<SampleTimings>,
    /// MAG-level coverage stats ready to write. Empty in contig mode.
    mag_results: MagResultVec,
}

fn process_samples_parallel(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    flags: ModuleFlags,
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    repeats: &[RepeatsData],
    annotations: Vec<FeatureAnnotation>,
    output_db: &Path,
    mut timing_file: Option<fs::File>,
    mut pre_timings: PreSampleTimings,
    mag_contig_map: &[(String, Vec<String>)],
    run_start: std::time::Instant,
    create_indexes: bool,
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let is_tty = atty::is(Stream::Stderr);
    let start_time = std::time::Instant::now();

    // Sequential sample processing when:
    // - Single thread: no parallelism needed
    // - Many contigs (≥500): inner per-contig parallelism already saturates threads,
    //   outer sample parallelism causes all samples to progress simultaneously
    //   without any completing. Sequential ensures samples finish one-by-one.
    if config.threads == 1 || contigs.len() >= 500 {
        return process_samples_sequential(bam_files, contigs, flags, config, circularity_map, db_writer, is_tty, repeats, annotations, output_db, timing_file, pre_timings, mag_contig_map, run_start, create_indexes);
    }

    // Adaptive channel size based on memory pressure from contig count
    // - Metagenomic (≥500 contigs): channel_size=1 → process one sample at a time (prevents OOM)
    // - Normal (<500 contigs): channel_size=threads → full parallelism
    let channel_size = if contigs.len() >= 500 {
        1  // Large samples → process one at a time to prevent OOM
    } else {
        config.threads  // Normal samples → full parallelism
    };

    // Sample-parallel mode: producer-consumer with bounded channel
    let mp = MultiProgress::new();

    // For TTY: use interactive progress bars
    // For non-TTY (SLURM logs): use hidden bars and print log messages
    let process_pb = if is_tty {
        let pb = mp.add(ProgressBar::new(total as u64));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        mp.add(ProgressBar::hidden())
    };

    let write_pb = if is_tty {
        let pb = mp.add(ProgressBar::new(total as u64));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.yellow} Writing:     [{elapsed_precise}] [{bar:40.yellow/red}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        mp.add(ProgressBar::hidden())
    };

    // Bounded channel limits memory based on dataset characteristics
    let (tx, rx) = sync_channel::<SampleResult>(channel_size);

    let completed_count = AtomicUsize::new(0);
    let failed_count = AtomicUsize::new(0);

    // Clone values needed by writer thread
    let write_pb_clone = write_pb.clone();

    let timing_threads = config.threads;
    let enable_timing_writer = config.enable_timing;
    let terminal_repeats_writer = config.terminal_repeats.clone();

    // MAG threshold filtering data (moved into writer thread)
    let is_mag_mode = config.view_mode == ViewMode::Mag;
    let min_af = config.min_aligned_fraction;
    let min_cov = config.min_coverage_depth;
    let mag_map_owned: Vec<(String, Vec<String>)> = mag_contig_map.to_vec();
    // Build MagMemberLookup before writer thread takes db_writer
    let mag_member_lookup: Option<MagMemberLookup> = if is_mag_mode {
        Some(crate::mag_blob::build_mag_member_lookup(&db_writer)?)
    } else {
        None
    };

    // Build shared lookups once — avoids rebuilding per sample (critical for 2M+ contigs).
    let contig_by_name = database_contig_lookup(contigs, &db_writer)?;
    let t_cds = std::time::Instant::now();
    drop(annotations);
    // Finish fallible annotation-loader setup before spawning the database writer.
    let cds_indexes = DatabaseCdsSource::from_connection(db_writer.lock_conn()?.try_clone()?, flags.mapping_metrics)?;
    pre_timings.cds_build_secs = t_cds.elapsed().as_secs_f64();
    pre_timings.cds_build_rss_mb = get_rss_mb();
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  CDS index build      : {:>10.3} s  [RSS: {:.0} MB]", pre_timings.cds_build_secs, pre_timings.cds_build_rss_mb);
        let total_pre = pre_timings.parse_secs + pre_timings.db_create_secs + pre_timings.circ_secs
            + pre_timings.autoblast_secs + pre_timings.gc_secs[0] + pre_timings.repeat_blob_secs
            + pre_timings.write_mags_secs + pre_timings.cds_build_secs;
        let _ = writeln!(f, "  TOTAL pre-sample     : {:>10.3} s  [RSS: {:.0} MB]", total_pre, get_rss_mb());
        let _ = f.flush();
    }

    // Spawn dedicated writer thread
    let writer_handle = thread::spawn(move || -> Result<(usize, std::time::Duration, Vec<SampleTimings>, usize, usize, FinalizeTimings, f64, f64, f64, f64)> {
        let write_start = std::time::Instant::now();
        let mut written_count = 0usize;
        let mut mag_drops: usize = 0;
        let mut samples_no_coverage: usize = 0;
        let mut all_timings: Vec<SampleTimings> = Vec::new();
        let mut timing_file = timing_file;

        // Receive and write samples as they arrive
        for result in rx {

            let identity = SampleIdentity { name: &result.sample_name, sequencing_type: result.sequencing_type,
                total_reads: result.total_reads, mapped_reads: result.mapped_reads, bai_mapped_reads: result.bai_mapped_reads, circular: result.is_circular };
            let (data, dropped) = SampleData {
                presences: result.presences,
                packaging: result.packaging,
                misassembly: result.misassembly,
                microdiversity: result.microdiversity,
                side_misassembly: result.side_misassembly,
                topology: result.topology,
                feature_blobs: result.feature_blobs,
                mag_results: result.mag_results,
            }.filter(is_mag_mode, &mag_map_owned, min_af, min_cov);
            mag_drops += dropped;

            if data.presences.is_empty() {
                samples_no_coverage += 1;
                let write_start = std::time::Instant::now();
                identity.insert(&db_writer,true,min_af,min_cov)?;
                db_writer.sample_complete(&result.sample_name)?;
                if let Some(mut t) = result.timings {
                    SampleWriteDetails { insert_ns: write_start.elapsed().as_nanos() as u64, ..Default::default() }.apply(&mut t, write_start.elapsed());
                    if let Some(ref mut f) = timing_file {
                        write_sample_timing(f, &t, timing_threads)?;
                    }
                    all_timings.push(t);
                }
                drop(data);
                log_sample_cleanup(&result.sample_name, enable_timing_writer, &mut timing_file);
                write_pb_clone.inc(1);
                continue;
            }

            // Insert sample
            let sample_write_start = std::time::Instant::now();
            let t_insert = std::time::Instant::now();
            let sample_id=identity.insert(&db_writer,false,min_af,min_cov)?;
            let insert_ns = t_insert.elapsed().as_nanos() as u64;

            // Write all data for this sample
            let mut details = SampleWriteDetails { insert_ns, ..Default::default() };
            if let Err(error)=write_sample_payload(&db_writer,&identity,sample_id,&data,enable_timing_writer,&mut details) {
                db_writer.sample_failed(&result.sample_name,&error.to_string())?;
                return Err(error);
            }
            written_count+=1;
            write_pb_clone.set_message(result.sample_name.clone());
            write_pb_clone.inc(1);

            // Write per-sample timings incrementally to the log file
            if let Some(mut t) = result.timings {
                details.apply(&mut t, sample_write_start.elapsed());
                if let Some(ref mut f) = timing_file {
                    let _ = write_sample_timing(f, &t, timing_threads);
                }
                all_timings.push(t);
            }

            drop(data);
            checkpoint_sample(&db_writer, &result.sample_name, enable_timing_writer, &mut timing_file)?;
            db_writer.sample_complete(&result.sample_name)?;
            log_sample_cleanup(&result.sample_name, enable_timing_writer, &mut timing_file);

        }

        // Build query indexes
        let index_secs = if create_indexes {
            db_writer.create_serve_indexes()?
        } else {
            0.0
        };

        // Finalize database
        let rss_before_finalize = get_rss_mb();
        let t_finalize = std::time::Instant::now();
        if let Some(search) = &terminal_repeats_writer { search.persist(&db_writer)?; }
        let ft = db_writer.finalize()?;
        let finalize_secs = t_finalize.elapsed().as_secs_f64();
        let rss_after_finalize = get_rss_mb();
        write_pb_clone.finish();

        Ok((written_count, write_start.elapsed(), all_timings, mag_drops, samples_no_coverage, ft, index_secs, finalize_secs, rss_before_finalize, rss_after_finalize))
    });

    // Process samples in parallel, sending to channel immediately
    // send() blocks if channel is full (backpressure)
    let mag_member_lookup_ref = mag_member_lookup.as_ref();
    bam_files.par_iter().for_each(|bam_path| {
        let sample_start = std::time::Instant::now();

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        match process_sample(bam_path, flags, config, sample_circular, repeats, &contig_by_name, &cds_indexes, mag_member_lookup_ref) {
            Ok(mut r) => {
                let sample_time = sample_start.elapsed().as_secs_f64();
                if let Some(ref mut t) = r.timings {
                    t.total = sample_start.elapsed();
                }
                completed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("{} ({:.2}s)", r.sample_name, sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }

                if tx.send(SampleResult {
                    sample_name: r.sample_name,
                    sequencing_type: r.seq_type,
                    total_reads: r.total_reads,
                    mapped_reads: r.mapped_reads,
                    bai_mapped_reads: r.bai_mapped_reads,
                    is_circular: r.is_circular,
                    feature_blobs: r.feature_blobs,
                    presences: r.presences,
                    packaging: r.packaging,
                    misassembly: r.misassembly,
                    microdiversity: r.microdiversity,
                    side_misassembly: r.side_misassembly,
                    topology: r.topology,
                    timings: r.timings,
                    mag_results: r.mag_results,
                }).is_err() {failed_count.fetch_add(1,Ordering::SeqCst);}
            }
            Err(e) => {
                eprintln!("\nError processing {}: {}", bam_path.display(), e);
                let sample_time = sample_start.elapsed().as_secs_f64();
                completed_count.fetch_add(1, Ordering::SeqCst);
                failed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("{} (failed after {:.2}s)", bam_path.file_stem().unwrap_or_default().to_string_lossy(), sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
        }
    });

    // Drop sender to signal writer thread that processing is complete
    drop(tx);

    process_pb.finish();
    let processing_time = start_time.elapsed();

    // Wait for writer thread to finish
    let (written_count, writing_time, all_timings, mag_drops, samples_no_coverage, ft, index_secs, finalize_secs, rss_before_finalize, rss_after_finalize) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    if samples_no_coverage > 0 {
        eprintln!("{} samples recorded with no coverage results retained after filtering", samples_no_coverage);
    }
    let dropped_contigs = config.contig_drops_counter.load(Ordering::Relaxed);
    if dropped_contigs > 0 {
        eprintln!("{} contig/sample pairs failed coverage thresholds", dropped_contigs);
    }
    if mag_drops > 0 {
        eprintln!("{} MAG/sample pairs failed coverage thresholds", mag_drops);
    }
    if ft.zero_coverage_contigs > 0 {
        eprintln!("INFO: {} contigs had no coverage in any sample (Contig rows retained, Number_of_samples=0)", ft.zero_coverage_contigs);
    }
    if ft.zero_coverage_mags > 0 {
        eprintln!("INFO: {} MAGs had no coverage in any sample (MAG rows retained, Number_of_samples=0)", ft.zero_coverage_mags);
    }

    // Clear MultiProgress to avoid duplicate bar display
    mp.clear().ok();

    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);

    // Append grand totals to the timing log
    if config.enable_timing && !all_timings.is_empty() {
        let log_path = timing_log_path(output_db);
        match fs::OpenOptions::new().append(true).open(&log_path) {
            Ok(mut f) => {
                if let Err(e) = finish_timing_log(&mut f, &all_timings, elapsed, &pre_timings, run_start, &ft, index_secs, finalize_secs, rss_before_finalize, rss_after_finalize) {
                    eprintln!("Warning: failed to write timing totals: {}", e);
                }
            }
            Err(e) => eprintln!("Warning: could not reopen timing log: {}", e),
        }
    }

    Ok(ProcessResult {
        samples_processed: written_count,
        samples_failed: failed,
        samples_no_coverage,
        total_time_secs: elapsed.as_secs_f64(),
        processing_time_secs: processing_time.as_secs_f64(),
        writing_time_secs: writing_time.as_secs_f64(),
    })
}

enum MagStreamMessage {
    Start { sample_name: String, seq_type: SequencingType, is_circular: bool },
    Group(MagGroupResult),
}

#[derive(Default)]
struct MagWriteBatch {
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    mags: MagResultVec,
    payload_bytes: usize,
}

impl MagWriteBatch {
    fn push(&mut self, r: MagGroupResult) {
        self.payload_bytes += r.contig_blobs.iter().map(|(_, _, b)|
            b.zoom.len() + b.chunks.iter().map(Vec::len).sum::<usize>()).sum::<usize>();
        self.presences.extend(r.presences);
        self.packaging.extend(r.packaging);
        self.misassembly.extend(r.misassembly);
        self.microdiversity.extend(r.microdiversity);
        self.side_misassembly.extend(r.side_misassembly);
        self.topology.extend(r.topology);
        self.blobs.extend(r.contig_blobs);
        self.mags.push((r.mag_name, r.mag_id, r.mag_cov_stats));
    }
    fn is_full(&self) -> bool {
        self.mags.len() >= 64 || self.presences.len() >= 32768 || self.payload_bytes >= 128 * 1024 * 1024
    }
}

#[derive(Default)]
struct MagStreamWrites {
    sample_id: Option<i64>,
    insert_ns: u64,
    contig_ns: u64,
    mag_ns: u64,
    detail: crate::db::SampleWriteTimings,
    write_time: std::time::Duration,
    coverage_rows: usize,
    mag_drops: usize,
}

fn checkpoint_byte_limit()->Result<usize> {
    let mib=std::env::var("THEBIGBAM_CHECKPOINT_MIB").ok().map(|s|s.parse::<usize>()).transpose()?.unwrap_or(1024);
    anyhow::ensure!(mib>0,"THEBIGBAM_CHECKPOINT_MIB must be positive");
    mib.checked_mul(1024*1024).context("Checkpoint byte budget overflow")
}

fn flush_mag_batch(writer: &DbWriter, sample_name: &str, sample_id: i64, circular: bool,
    timing: bool, batch: &mut MagWriteBatch, writes: &mut MagStreamWrites) -> Result<()> {
    if batch.mags.is_empty() { return Ok(()); }
    let start = std::time::Instant::now();
    let detail = writer.write_sample_batch_timed(sample_name, &batch.presences, &batch.packaging,
        &batch.misassembly, &batch.microdiversity, &batch.side_misassembly, &batch.topology,
        &batch.blobs, circular, timing, false)?;
    writes.contig_ns += start.elapsed().as_nanos() as u64;
    writes.detail.add(&detail);
    let start = std::time::Instant::now();
    {
        let conn = writer.lock_conn()?;
        writer.write_mag_coverages(&conn, sample_id, &batch.mags)?;
    }
    writes.mag_ns += start.elapsed().as_nanos() as u64;
    writes.coverage_rows += batch.presences.len();
    *batch = MagWriteBatch::default();
    Ok(())
}

/// A bounded queue overlaps per-MAG processing and database writing. Appenders
/// flush batches; checkpoints are periodic rather than once per small batch.
fn process_mag_sample_streamed(
    bam_path: &Path, flags: ModuleFlags, config: &ProcessConfig, circular: bool,
    repeats: &[RepeatsData], contig_by_name: &dyn ContigLookup,
    cds_indexes: &dyn CdsSource, lookup: &MagMemberLookup, writer: &DbWriter,
) -> Result<(SampleProcessResult, MagStreamWrites)> {
    let outcome: Result<(SampleProcessResult, MagStreamWrites)> = thread::scope(|scope| {
        let (tx, rx) = sync_channel::<MagStreamMessage>(2);
        let write_handle = scope.spawn(move || -> Result<MagStreamWrites> {
            let mut writes = MagStreamWrites::default();
            let mut sample_name = String::new();
            let mut sample_circular = circular;
            let mut batch = MagWriteBatch::default();
            let checkpoint_limit=checkpoint_byte_limit()?;
            let mut checkpoint_bytes = 0usize;
            let mut checkpoint_rows = 0usize;
            for message in rx {
                match message {
                    MagStreamMessage::Start { sample_name: name, seq_type, is_circular } => {
                        let start = std::time::Instant::now();
                        let id = writer.insert_sample(&name, seq_type.as_str(), 0, 0, is_circular)?;
                        writes.insert_ns = start.elapsed().as_nanos() as u64;
                        writes.write_time += start.elapsed();
                        writes.sample_id = Some(id);
                        sample_name = name;
                        sample_circular = is_circular;
                    }
                    MagStreamMessage::Group(group) => {
                        let af = group.mag_cov_stats.aligned_fraction_pct_x10 as f64 /
                            crate::types::get_column_scale("Coverage", "Aligned_fraction_percentage");
                        let mean = group.mag_cov_stats.coverage_mean_x10 as f64 /
                            crate::types::get_column_scale("Coverage", "Coverage_mean");
                        if mag_coverage_fails(af, mean, config.min_aligned_fraction, config.min_coverage_depth) {
                            writes.mag_drops += 1;
                            continue;
                        }
                        if group.presences.is_empty() { continue; }
                        batch.push(group);
                        if batch.is_full() {
                            let start = std::time::Instant::now();
                            checkpoint_bytes += batch.payload_bytes;
                            checkpoint_rows += batch.presences.len();
                            flush_mag_batch(writer, &sample_name, writes.sample_id.context("Missing sample start")?,
                                sample_circular, config.enable_timing, &mut batch, &mut writes)?;
                            if checkpoint_bytes >= checkpoint_limit || checkpoint_rows >= 1_000_000 {
                                let checkpoint = std::time::Instant::now();
                                writer.checkpoint()?;
                                if config.enable_timing { writes.detail.checkpoint_ns += checkpoint.elapsed().as_nanos() as u64; }
                                checkpoint_bytes = 0;
                                checkpoint_rows = 0;
                            }
                            writes.write_time += start.elapsed();
                        }
                    }
                }
            }
            if let Some(id) = writes.sample_id {
                let start = std::time::Instant::now();
                flush_mag_batch(writer, &sample_name, id, sample_circular, config.enable_timing, &mut batch, &mut writes)?;
                writes.write_time += start.elapsed();
            }
            Ok(writes)
        });
        let sink = move |message| tx.send(message).map_err(|_| anyhow::anyhow!("MAG database writer stopped"));
        let result = process_sample_with_mag_sink(bam_path, flags, config, circular, repeats,
            contig_by_name, cds_indexes, Some(lookup), Some(&sink));
        drop(sink);
        let mut writes = write_handle.join().map_err(|_| anyhow::anyhow!("MAG database writer panicked"))??;
        let result = result?;
        let id = writes.sample_id.context("MAG stream did not start")?;
        let (total,mapped)=(result.total_reads,result.mapped_reads);
        let start = std::time::Instant::now();
        writer.update_sample_read_counts(id, total, mapped)?;
        writes.write_time += start.elapsed();
        Ok((result, writes))
    });
    if let Err(error)=&outcome {
        let name=bam_path.file_stem().unwrap_or_default().to_string_lossy().replace("_with_MD", "");
        writer.sample_failed(&name,&error.to_string())?;
    }
    outcome
}

/// Process samples one at a time. Contig view writes after processing; MAG
/// view overlaps per-MAG work with a bounded database-writing queue.
fn process_samples_sequential(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    flags: ModuleFlags,
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    is_tty: bool,
    repeats: &[RepeatsData],
    annotations: Vec<FeatureAnnotation>,
    _output_db: &Path,
    timing_file: Option<fs::File>,
    mut pre_timings: PreSampleTimings,
    mag_contig_map: &[(String, Vec<String>)],
    run_start: std::time::Instant,
    create_indexes: bool,
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let start_time = std::time::Instant::now();
    let is_mag_mode = config.view_mode == ViewMode::Mag;
    let mut mag_drops: usize = 0;
    let mut samples_no_coverage: usize = 0;

    let pb = if is_tty {
        let pb = ProgressBar::new(total as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} Processing:  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        pb
    } else {
        ProgressBar::hidden()
    };

    let mut processed = 0usize;
    let mut failed = 0usize;
    let mut processing_time_total = std::time::Duration::ZERO;
    let mut writing_time_total = std::time::Duration::ZERO;
    let mut all_timings: Vec<SampleTimings> = Vec::new();
    let mut timing_file = timing_file;

    let mag_member_lookup: Option<MagMemberLookup> = if is_mag_mode {
        Some(crate::mag_blob::build_mag_member_lookup(&db_writer)?)
    } else {
        None
    };

    // Build shared lookups once — avoids rebuilding per sample (critical for 2M+ contigs).
    let contig_by_name = database_contig_lookup(contigs, &db_writer)?;
    let t_cds = std::time::Instant::now();
    drop(annotations);
    let cds_indexes = DatabaseCdsSource::from_connection(db_writer.lock_conn()?.try_clone()?, flags.mapping_metrics)?;
    pre_timings.cds_build_secs = t_cds.elapsed().as_secs_f64();
    pre_timings.cds_build_rss_mb = get_rss_mb();
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "  CDS index build      : {:>10.3} s  [RSS: {:.0} MB]", pre_timings.cds_build_secs, pre_timings.cds_build_rss_mb);
        let total_pre = pre_timings.parse_secs + pre_timings.db_create_secs + pre_timings.circ_secs
            + pre_timings.autoblast_secs + pre_timings.gc_secs[0] + pre_timings.repeat_blob_secs
            + pre_timings.write_mags_secs + pre_timings.cds_build_secs;
        let _ = writeln!(f, "  TOTAL pre-sample     : {:>10.3} s  [RSS: {:.0} MB]", total_pre, get_rss_mb());
        let _ = f.flush();
    }

    for bam_path in bam_files {
        let sample_start = std::time::Instant::now();
        let sample_memory_name = bam_path.file_stem().unwrap_or_default().to_string_lossy().replace("_with_MD", "");

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        if let Some(lookup) = mag_member_lookup.as_ref() {
            let outcome = process_mag_sample_streamed(bam_path, flags, config, sample_circular,
                repeats, &contig_by_name, &cds_indexes, lookup, &db_writer);
            match outcome {
                Ok((mut result, writes)) => {
                    processing_time_total += sample_start.elapsed();
                    writing_time_total += writes.write_time;
                    // Absent groups have zero coverage; retain existing drop counts.
                    let known: HashSet<&str> = result.mag_results.iter().map(|(name, _, _)| name.as_str()).collect();
                    let missing_drops = if mag_coverage_fails(0.0, 0.0, config.min_aligned_fraction, config.min_coverage_depth) {
                        mag_contig_map.iter().filter(|(name, _)| !known.contains(name.as_str())).count()
                    } else { 0 };
                    mag_drops += writes.mag_drops + missing_drops;
                    if writes.coverage_rows == 0 {
                        samples_no_coverage += 1;
                        writing::report_empty_sample(&result.sample_name, result.bai_mapped_reads, true,
                            config.min_aligned_fraction, config.min_coverage_depth);
                    }
                    else { processed += 1; }
                    if let Some(ref mut t) = result.timings {
                        t.db_overlapped = true;
                        t.db_write = writes.write_time;
                        t.insert_sample_ns = writes.insert_ns;
                        t.write_contig_data_ns = writes.contig_ns;
                        t.mag_write_ns = writes.mag_ns;
                        t.db_detail = writes.detail;
                        t.total = sample_start.elapsed();
                        if let Some(ref mut f) = timing_file { write_sample_timing(f, t, config.threads)?; }
                        all_timings.push(t.clone());
                    }
                    // Bound WAL memory at sample completion as before.
                    let checkpoint_start = std::time::Instant::now();
                    let rss_before = if config.enable_timing { get_rss_mb() } else { 0.0 };
                    db_writer.checkpoint()?;
                    db_writer.sample_complete(&result.sample_name)?;
                    if config.enable_timing {
                        use std::io::Write;
                        if let Some(ref mut f) = timing_file {
                            writeln!(f, "  checkpoint           : {:>10.3} s  [RSS: {:.0} → {:.0} MB]",
                                checkpoint_start.elapsed().as_secs_f64(), rss_before, get_rss_mb())?;
                        }
                    }
                    let msg = format!("{} ({:.2}s)", result.sample_name, sample_start.elapsed().as_secs_f64());
                    pb.set_message(msg.clone());
                    pb.inc(1);
                    if !is_tty { eprintln!("{}", msg); }
                }
                Err(error) => {
                    failed += 1;
                    processing_time_total += sample_start.elapsed();
                    pb.inc(1);
                    let _=db_writer.sample_failed(&sample_memory_name,&error.to_string());
                    eprintln!("Error processing {}: {}", bam_path.display(), error);
                }
            }
            log_cds_cache(&cds_indexes,&sample_memory_name,config.enable_timing,&mut timing_file);
            log_sample_cleanup(&sample_memory_name, config.enable_timing, &mut timing_file);
            continue;
        }
        match process_sample(bam_path, flags, config, sample_circular, repeats, &contig_by_name, &cds_indexes, mag_member_lookup.as_ref()) {
            Ok(mut r) => {
                let process_elapsed = sample_start.elapsed();
                processing_time_total += process_elapsed;

                let identity = SampleIdentity { name: &r.sample_name, sequencing_type: r.seq_type,
                    total_reads: r.total_reads, mapped_reads: r.mapped_reads, bai_mapped_reads: r.bai_mapped_reads, circular: r.is_circular };
                let (data, dropped) = SampleData {
                    presences: r.presences,
                    packaging: r.packaging,
                    misassembly: r.misassembly,
                    microdiversity: r.microdiversity,
                    side_misassembly: r.side_misassembly,
                    topology: r.topology,
                    feature_blobs: r.feature_blobs,
                    mag_results: r.mag_results,
                }.filter(is_mag_mode, mag_contig_map, config.min_aligned_fraction, config.min_coverage_depth);
                mag_drops += dropped;

                if data.presences.is_empty() {
                    samples_no_coverage += 1;
                    let write_start = std::time::Instant::now();
                    identity.insert(&db_writer,true,config.min_aligned_fraction,config.min_coverage_depth)?;
                    db_writer.sample_complete(&r.sample_name)?;
                    writing_time_total += write_start.elapsed();
                    if let Some(ref mut t) = r.timings {
                        SampleWriteDetails { insert_ns: write_start.elapsed().as_nanos() as u64, ..Default::default() }.apply(t, write_start.elapsed());
                        t.total = sample_start.elapsed();
                        if let Some(ref mut f) = timing_file {
                            write_sample_timing(f, t, config.threads)?;
                        }
                        all_timings.push(t.clone());
                    }
                    drop(data);
                    log_sample_cleanup(&sample_memory_name, config.enable_timing, &mut timing_file);
                    let msg = format!("{} ({:.2}s)", r.sample_name, sample_start.elapsed().as_secs_f64());
                    pb.set_message(msg.clone());
                    pb.inc(1);
                    if !is_tty { eprintln!("{}", msg); }
                    continue;
                }

                let write_start = std::time::Instant::now();

                // Insert sample
                let t_insert = std::time::Instant::now();
                let sample_id = match identity.insert(&db_writer, false, config.min_aligned_fraction, config.min_coverage_depth) {
                    Ok(id) => id,
                    Err(e) => {
                        eprintln!("\nError inserting sample {}: {}", r.sample_name, e);
                        failed += 1;
                        drop(data);
                        log_sample_cleanup(&sample_memory_name, config.enable_timing, &mut timing_file);
                        pb.inc(1);
                        continue;
                    }
                };
                let insert_ns = t_insert.elapsed().as_nanos() as u64;

                // Write sample data
                let mut write_succeeded=false;
                let mut details = SampleWriteDetails { insert_ns, ..Default::default() };
                if let Err(e) = write_sample_payload(&db_writer, &identity, sample_id, &data, config.enable_timing, &mut details) {
                    eprintln!("\nError writing data for {}: {}", r.sample_name, e);
                    db_writer.sample_failed(&r.sample_name,&e.to_string())?;
                    failed += 1;
                } else {
                    write_succeeded=true;
                    processed += 1;
                }

                writing_time_total += write_start.elapsed();

                // Write per-sample timings incrementally to the log file
                if let Some(ref mut t) = r.timings {
                    details.apply(t, write_start.elapsed());
                    t.total = sample_start.elapsed();
                    if let Some(ref mut f) = timing_file {
                        let _ = write_sample_timing(f, t, config.threads);
                    }
                    all_timings.push(t.clone());
                }

                drop(data);
                checkpoint_sample(&db_writer, &r.sample_name, config.enable_timing, &mut timing_file)?;
                if write_succeeded {db_writer.sample_complete(&r.sample_name)?;}
                let total_sample_time = sample_start.elapsed().as_secs_f64();
                let msg = format!("{} ({:.2}s)", r.sample_name, total_sample_time);
                pb.set_message(msg.clone());
                pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
            Err(e) => {
                db_writer.sample_failed(&sample_memory_name,&e.to_string())?;
                eprintln!("\nError processing {}: {}", bam_path.display(), e);
                failed += 1;
                processing_time_total += sample_start.elapsed();
                let msg = "ERR".to_string();
                pb.set_message(msg.clone());
                pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
        }
        log_cds_cache(&cds_indexes,&sample_memory_name,config.enable_timing,&mut timing_file);
        log_sample_cleanup(&sample_memory_name, config.enable_timing, &mut timing_file);
    }

    if samples_no_coverage > 0 {
        eprintln!("{} samples recorded with no coverage results retained after filtering", samples_no_coverage);
    }
    let dropped_contigs = config.contig_drops_counter.load(Ordering::Relaxed);
    if dropped_contigs > 0 {
        eprintln!("{} contig/sample pairs failed coverage thresholds", dropped_contigs);
    }
    if mag_drops > 0 {
        eprintln!("{} MAG/sample pairs failed coverage thresholds", mag_drops);
    }

    // Build query indexes
    let index_secs = if create_indexes {
        db_writer.create_serve_indexes()?
    } else {
        0.0
    };

    // Finalize database
    let rss_before_finalize = get_rss_mb();
    let t_finalize = std::time::Instant::now();
    if let Some(search) = &config.terminal_repeats { search.persist(&db_writer)?; }
    let ft = db_writer.finalize()?;
    let finalize_secs = t_finalize.elapsed().as_secs_f64();
    if ft.zero_coverage_contigs > 0 {
        eprintln!("INFO: {} contigs had no coverage in any sample (Contig rows retained, Number_of_samples=0)", ft.zero_coverage_contigs);
    }
    if ft.zero_coverage_mags > 0 {
        eprintln!("INFO: {} MAGs had no coverage in any sample (MAG rows retained, Number_of_samples=0)", ft.zero_coverage_mags);
    }
    let rss_after_finalize = get_rss_mb();
    pb.finish_with_message("Done");
    if !is_tty {
        eprintln!("Done ({:.2}s)", start_time.elapsed().as_secs_f64());
    }

    // Append grand totals to the timing log
    if config.enable_timing && !all_timings.is_empty() {
        if let Some(ref mut f) = timing_file {
            if let Err(e) = finish_timing_log(f, &all_timings, start_time.elapsed(), &pre_timings, run_start, &ft, index_secs, finalize_secs, rss_before_finalize, rss_after_finalize) {
                eprintln!("Warning: failed to write timing totals: {}", e);
            }
        }
    }

    Ok(ProcessResult {
        samples_processed: processed,
        samples_failed: failed,
        samples_no_coverage,
        total_time_secs: start_time.elapsed().as_secs_f64(),
        processing_time_secs: processing_time_total.as_secs_f64(),
        writing_time_secs: writing_time_total.as_secs_f64(),
    })
}

/// Snapshot after sample-local payloads are released. RSS includes reusable
/// buffers, bounded caches, DuckDB residency, and concurrent work in parallel mode.
fn log_cds_cache(source:&DatabaseCdsSource,sample:&str,enabled:bool,log:&mut Option<fs::File>) {
    if !enabled {return;}
    let (loads,ns,oversized,bytes)=source.diagnostics();
    let line=format!("[Sample {sample}] CDS cache cumulative: {loads} SQL loads, {:.3}s loading, {oversized} oversized batches, {:.1} MiB retained",ns as f64/1e9,bytes as f64/(1024.0*1024.0));
    if let Some(file)=log {use std::io::Write;let _=writeln!(file,"{line}");let _=file.flush();}
}

fn log_sample_cleanup(sample: &str, enabled: bool, log: &mut Option<fs::File>) {
    if !enabled { return; }
    let (rss, peak) = get_rss_and_peak_mb();
    if let Some(file) = log {
        use std::io::Write;
        let _ = writeln!(file, "  Memory after sample payload cleanup [{}]: RSS {:.0} MiB; process peak {:.0} MiB", sample, rss, peak);
        let _ = file.flush();
    }
}

fn ns_to_s(ns: u64) -> f64 {
    ns as f64 / 1e9
}

fn get_rss_mb() -> f64 {
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    if let Some(kb_str) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb_str.parse::<u64>() {
                            return kb as f64 / 1024.0;
                        }
                    }
                }
            }
        }
        0.0
    }
    #[cfg(not(target_os = "linux"))]
    { 0.0 }
}

fn get_rss_and_peak_mb() -> (f64, f64) {
    #[cfg(target_os = "linux")]
    {
        let mut rss = 0.0;
        let mut hwm = 0.0;
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    if let Some(kb_str) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb_str.parse::<u64>() {
                            rss = kb as f64 / 1024.0;
                        }
                    }
                } else if line.starts_with("VmHWM:") {
                    if let Some(kb_str) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb_str.parse::<u64>() {
                            hwm = kb as f64 / 1024.0;
                        }
                    }
                }
            }
        }
        (rss, hwm)
    }
    #[cfg(not(target_os = "linux"))]
    { (0.0, 0.0) }
}

fn timing_log_path(output_db: &Path) -> PathBuf {
    let mut p = output_db.as_os_str().to_owned();
    p.push(".timings.log");
    PathBuf::from(p)
}

/// Open the timing log and write the header. Returns the open file handle
/// so callers can append per-sample sections as they finish.
fn open_timing_log(output_db: &Path, threads: usize, n_samples: usize, command_line: &str) -> Result<fs::File> {
    use std::io::Write;
    let log_path = timing_log_path(output_db);
    let mut f = fs::File::create(&log_path)
        .with_context(|| format!("Failed to create timing log: {}", log_path.display()))?;

    writeln!(f, "Timing log for command: {}", command_line)?;
    writeln!(f)?;
    writeln!(f, "  Output DB    : {}", output_db.display())?;
    writeln!(f, "  Threads      : {}", threads)?;
    writeln!(f, "  Samples      : {}", n_samples)?;
    writeln!(f, "  Memory       : RSS tracked via /proc/self/status")?;
    let (duckdb_mib, detected_gb, mem_source) = crate::db::duckdb_memory_budget()?;
    writeln!(f, "  DuckDB budget: {} MiB at startup ({} GiB from {}; recalculated on connection open)", duckdb_mib, detected_gb, mem_source)?;
    let _ = f.flush();

    eprintln!("Timing log: {}", log_path.display());
    Ok(f)
}

/// Append one sample's timing breakdown to an already-open timing log.
fn write_sample_timing(f: &mut fs::File, t: &SampleTimings, threads: usize) -> Result<()> {
    use std::io::Write;
    writeln!(f)?;
    writeln!(f, "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")?;
    writeln!(f, "### Sample: {}  [RSS: {:.0} MB]", t.sample_name, get_rss_mb())?;
    writeln!(f)?;
    writeln!(f, "  RSS before reader open: {:>10.0} MiB", t.rss_before_reader_open_mib)?;
    writeln!(f, "  BAM path             : {}", t.bam_path)?;
    writeln!(f, "  Contigs in BAM       : {}", t.n_contigs_in_bam)?;
    writeln!(f, "  Contigs processed    : {}", t.n_contigs_processed)?;
    writeln!(f, "  Total reads          : {}", t.total_reads)?;
    writeln!(f, "  ---- Wall-clock phases ----")?;
    writeln!(f, "  Header read          : {:>10.3} s", t.header_read.as_secs_f64())?;
    writeln!(f, "  Total read count     : {:>10.3} s", t.total_read_count.as_secs_f64())?;
    writeln!(f, "  Filter tids (pre-par): {:>10.3} s", t.filter_tids.as_secs_f64())?;
    writeln!(f, "  Par-iter wall        : {:>10.3} s", t.par_iter_wall.as_secs_f64())?;
    writeln!(f, "  Merge contig results : {:>10.3} s", t.merge_results.as_secs_f64())?;
    writeln!(f, "  DB write             : {:>10.3} s", t.db_write.as_secs_f64())?;
    writeln!(f, "  TOTAL per sample     : {:>10.3} s", t.total.as_secs_f64())?;
    if t.db_overlapped {
        writeln!(f, "  DB writes overlap processing; DB write is active writer time, not an additional wall phase.")?;
        writeln!(f, "  MAG queue send/wait  : {:>10.3} s (summed worker-seconds, included in par-iter wall)", ns_to_s(t.mag_queue_wait_ns))?;
    }
    writeln!(f, "  ---- Par-iter sub-phase sums (worker-seconds; {} readers, {} configured threads) ----", t.bam_reader_workers, threads)?;
    writeln!(f, "  BAM reader budget    : {} readers; estimated {} MiB each; start RSS {} MiB", t.bam_reader_workers, t.bam_reader_estimated_mib, t.bam_reader_start_rss_mib)?;
    writeln!(f, "  BAM reader open      : {:>10.3} s", ns_to_s(t.bam_open_ns))?;
    writeln!(f, "  process_contig_stream: {:>10.3} s", ns_to_s(t.streaming_ns))?;
    writeln!(f, "    . fetch seek       : {:>10.3} s", ns_to_s(t.stream_fetch_seek_ns))?;
    writeln!(f, "    . record loop (htslib+prep): {:>10.3} s", ns_to_s(t.stream_records_ns))?;
    writeln!(f, "    . process_read     : {:>10.3} s", ns_to_s(t.stream_process_read_ns))?;
    writeln!(f, "    . finalize+cov     : {:>10.3} s", ns_to_s(t.stream_finalize_ns))?;
    writeln!(f, "  Terminal search      : {:.3} s (worker-seconds, included in feature calculation)", ns_to_s(t.terminal.search_ns))?;
    writeln!(f, "    . DNA loading      : {:.3} s", ns_to_s(t.terminal.dna_ns))?;
    writeln!(f, "    . cache waits      : {:.3} s", ns_to_s(t.terminal.wait_ns))?;
    writeln!(f, "    . output parsing   : {:.3} s", ns_to_s(t.terminal.parse_ns))?;
    writeln!(f, "  Feature calculation  : {:>10.3} s", ns_to_s(t.feature_calc_ns))?;
    writeln!(f, "    . dense encoding   : {:>10.3} s", ns_to_s(t.dense_encode_ns))?;
    writeln!(f, "    . sparse encoding  : {:>10.3} s", ns_to_s(t.sparse_encode_ns))?;
    writeln!(f, "    . preparation + metrics: {:>7.3} s", ns_to_s(t.feature_calc_ns.saturating_sub(t.dense_encode_ns + t.sparse_encode_ns)))?;
    let sub_sum_ns = t.bam_open_ns + t.streaming_ns + t.feature_calc_ns;
    writeln!(f, "  (sum of above)       : {:>10.3} worker-seconds (overlapping work; not a wall-time estimate)", ns_to_s(sub_sum_ns))?;
    writeln!(f, "  ---- DB write sub-phases ----")?;
    writeln!(f, "  Insert sample        : {:>10.3} s", ns_to_s(t.insert_sample_ns))?;
    writeln!(f, "  Write contig data    : {:>10.3} s", ns_to_s(t.write_contig_data_ns))?;
    writeln!(f, "  Write MAG data       : {:>10.3} s", ns_to_s(t.mag_write_ns))?;
    write_db_detail_timing(f, &t.db_detail)?;
    let _ = f.flush();
    Ok(())
}

fn write_db_detail_timing(f: &mut fs::File, t: &crate::db::SampleWriteTimings) -> Result<()> {
    use std::io::Write;
    for (name, ns) in [
        ("coverage rows", t.coverage_ns), ("metric rows", t.metrics_ns),
        ("write checkpoints", t.checkpoint_ns), ("blob append", t.blob_append_ns),
        ("zoom flush", t.zoom_flush_ns), ("chunk flush", t.chunk_flush_ns),
    ] {
        writeln!(f, "    . {:<19}: {:>10.3} s", name, ns_to_s(ns))?;
    }
    writeln!(f, "  Blobs written        : {} zoom rows, {} chunk rows; {} zoom bytes, {} chunk bytes (encoded payload)",
        t.zoom_rows, t.chunk_rows, t.zoom_bytes, t.chunk_bytes)?;
    Ok(())
}

/// Append grand totals to the timing log and close it.
fn finish_timing_log(
    f: &mut fs::File,
    timings: &[SampleTimings],
    total_wall: std::time::Duration,
    pt: &PreSampleTimings,
    run_start: std::time::Instant,
    ft: &FinalizeTimings,
    index_secs: f64,
    finalize_secs: f64,
    rss_before_finalize: f64,
    rss_after_finalize: f64,
) -> Result<()> {
    use std::io::Write;
    writeln!(f)?;
    writeln!(f, "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")?;
    writeln!(f, "### Grand totals")?;
    writeln!(f)?;
    writeln!(f, "  ---- Pre-sample phases ----")?;
    writeln!(f, "  Annotation+FASTA parse: {:>9.3} s  [RSS: {:.0} MB]", pt.parse_secs, pt.parse_rss_mb)?;
    writeln!(f, "  DB create/open       : {:>10.3} s  [RSS: {:.0} MB]", pt.db_create_secs, pt.db_create_rss_mb)?;
    writeln!(f, "  Circularity detect   : {:>10.3} s  [RSS: {:.0} MB]", pt.circ_secs, pt.circ_rss_mb)?;
    writeln!(f, "  Autoblast            : {:>10.3} s  [RSS: {:.0} MB]", pt.autoblast_secs, pt.autoblast_rss_mb)?;
    writeln!(f, "  GC content + skew    : {:>10.3} s  [RSS: {:.0} MB]", pt.gc_secs[0], pt.gc_rss_mb)?;
    writeln!(f, "    compute            : {:>10.3} s", pt.gc_secs[1])?;
    writeln!(f, "    write blobs        : {:>10.3} s", pt.gc_secs[2])?;
    writeln!(f, "    update stats       : {:>10.3} s", pt.gc_secs[3])?;
    writeln!(f, "  Repeat blob conversion: {:>9.3} s  [RSS: {:.0} MB]", pt.repeat_blob_secs, pt.repeat_rss_mb)?;
    writeln!(f, "  Write MAGs           : {:>10.3} s  [RSS: {:.0} MB]", pt.write_mags_secs, pt.write_mags_rss_mb)?;
    writeln!(f, "  CDS index build      : {:>10.3} s  [RSS: {:.0} MB]", pt.cds_build_secs, pt.cds_build_rss_mb)?;
    writeln!(f, "  ---- Per-sample phases ----")?;
    let sum_dur = |g: fn(&SampleTimings) -> std::time::Duration| -> f64 {
        timings.iter().map(|t| g(t).as_secs_f64()).sum()
    };
    let sum_ns = |g: fn(&SampleTimings) -> u64| -> f64 {
        timings.iter().map(|t| ns_to_s(g(t))).sum()
    };
    writeln!(f, "  Header read          : {:>10.3} s", sum_dur(|t| t.header_read))?;
    writeln!(f, "  Total read count     : {:>10.3} s", sum_dur(|t| t.total_read_count))?;
    writeln!(f, "  Filter tids (pre-par): {:>10.3} s", sum_dur(|t| t.filter_tids))?;
    writeln!(f, "  Par-iter wall        : {:>10.3} s", sum_dur(|t| t.par_iter_wall))?;
    writeln!(f, "  Merge contig results : {:>10.3} s", sum_dur(|t| t.merge_results))?;
    writeln!(f, "  DB write             : {:>10.3} s", sum_dur(|t| t.db_write))?;
    writeln!(f, "  Total per-sample sum : {:>10.3} s", sum_dur(|t| t.total))?;
    writeln!(f, "  Wall-clock (calc)    : {:>10.3} s", total_wall.as_secs_f64())?;
    if timings.iter().any(|t| t.db_overlapped) {
        writeln!(f, "  Processing and DB write times overlap for streamed MAG samples.")?;
        writeln!(f, "  MAG queue send/wait  : {:>10.3} s (summed worker-seconds)", sum_ns(|t| t.mag_queue_wait_ns))?;
    }
    writeln!(f, "  ---- Par-iter sub-phase sums ----")?;
    writeln!(f, "  BAM reader open      : {:>10.3} s", sum_ns(|t| t.bam_open_ns))?;
    writeln!(f, "  process_contig_stream: {:>10.3} s", sum_ns(|t| t.streaming_ns))?;
    writeln!(f, "    . fetch seek       : {:>10.3} s", sum_ns(|t| t.stream_fetch_seek_ns))?;
    writeln!(f, "    . record loop (htslib+prep): {:>10.3} s", sum_ns(|t| t.stream_records_ns))?;
    writeln!(f, "    . process_read     : {:>10.3} s", sum_ns(|t| t.stream_process_read_ns))?;
    writeln!(f, "    . finalize+cov     : {:>10.3} s", sum_ns(|t| t.stream_finalize_ns))?;
    writeln!(f, "  Terminal BLAST       : {:.3} worker-seconds", ns_to_s(timings.iter().map(|t| t.terminal.search_ns).sum()))?;
    writeln!(f, "  Terminal DNA loading : {:.3} worker-seconds", ns_to_s(timings.iter().map(|t| t.terminal.dna_ns).sum()))?;
    writeln!(f, "  Terminal cache waits : {:.3} worker-seconds", ns_to_s(timings.iter().map(|t| t.terminal.wait_ns).sum()))?;
    writeln!(f, "  Terminal parsing     : {:.3} worker-seconds", ns_to_s(timings.iter().map(|t| t.terminal.parse_ns).sum()))?;
    writeln!(f, "  Feature calculation  : {:>10.3} s", sum_ns(|t| t.feature_calc_ns))?;
    writeln!(f, "    . dense encoding   : {:>10.3} s", sum_ns(|t| t.dense_encode_ns))?;
    writeln!(f, "    . sparse encoding  : {:>10.3} s", sum_ns(|t| t.sparse_encode_ns))?;
    writeln!(f, "  ---- DB write sub-phase sums ----")?;
    writeln!(f, "  Insert sample        : {:>10.3} s", sum_ns(|t| t.insert_sample_ns))?;
    writeln!(f, "  Write contig data    : {:>10.3} s", sum_ns(|t| t.write_contig_data_ns))?;
    writeln!(f, "  Write MAG data       : {:>10.3} s", sum_ns(|t| t.mag_write_ns))?;
    let db_total = crate::db::SampleWriteTimings {
        coverage_ns: timings.iter().map(|t| t.db_detail.coverage_ns).sum(),
        metrics_ns: timings.iter().map(|t| t.db_detail.metrics_ns).sum(),
        checkpoint_ns: timings.iter().map(|t| t.db_detail.checkpoint_ns).sum(),
        blob_append_ns: timings.iter().map(|t| t.db_detail.blob_append_ns).sum(),
        zoom_flush_ns: timings.iter().map(|t| t.db_detail.zoom_flush_ns).sum(),
        chunk_flush_ns: timings.iter().map(|t| t.db_detail.chunk_flush_ns).sum(),
        zoom_rows: timings.iter().map(|t| t.db_detail.zoom_rows).sum(),
        chunk_rows: timings.iter().map(|t| t.db_detail.chunk_rows).sum(),
        zoom_bytes: timings.iter().map(|t| t.db_detail.zoom_bytes).sum(),
        chunk_bytes: timings.iter().map(|t| t.db_detail.chunk_bytes).sum(),
    };
    write_db_detail_timing(f, &db_total)?;
    writeln!(f, "  ---- Post-sample phases ----")?;
    writeln!(f, "  RSS before finalize  : {:>10.0} MB", rss_before_finalize)?;
    writeln!(f, "  Build indexes        : {:>10.3} s", index_secs)?;
    writeln!(f, "  DB finalize          : {:>10.3} s", finalize_secs)?;
    writeln!(f, "    create_views       : {:>10.3} s", ft.create_views_secs)?;
    writeln!(f, "    cleanup_vars       : {:>10.3} s", ft.cleanup_vars_secs)?;
    writeln!(f, "    drop_empty_tables  : {:>10.3} s", ft.drop_empty_secs)?;
    writeln!(f, "    update_counts      : {:>10.3} s", ft.update_counts_secs)?;
    writeln!(f, "    CHECKPOINT         : {:>10.3} s", ft.checkpoint_secs)?;
    writeln!(f, "  RSS after finalize   : {:>10.0} MB", rss_after_finalize)?;
    writeln!(f, "  ---- Memory ----")?;
    let (current_rss, peak_rss) = get_rss_and_peak_mb();
    writeln!(f, "  Peak RSS (VmHWM)     : {:>10.0} MB", peak_rss)?;
    writeln!(f, "  Current RSS          : {:>10.0} MB", current_rss)?;
    writeln!(f)?;
    let pre_sample_timed = pt.parse_secs + pt.db_create_secs + pt.circ_secs
        + pt.autoblast_secs + pt.gc_secs[0] + pt.repeat_blob_secs
        + pt.write_mags_secs;
    let sample_phase_wall = total_wall.as_secs_f64();
    let total_timed = pre_sample_timed + sample_phase_wall;
    let real_wall_secs = run_start.elapsed().as_secs_f64();
    let untimed = real_wall_secs - total_timed;
    writeln!(f, "  ---- Timed vs real ----")?;
    writeln!(f, "  Pre-sample timed     : {:>10.3} s", pre_sample_timed)?;
    writeln!(f, "  Sample phase (wall)  : {:>10.3} s  (CDS build + samples + finalize)", sample_phase_wall)?;
    writeln!(f, "  Total accounted      : {:>10.3} s", total_timed)?;
    writeln!(f, "  Wall-clock real      : {:>10.3} s", real_wall_secs)?;
    writeln!(f, "  Untimed (gap)        : {:>10.3} s  ({:.1}%)", untimed, if real_wall_secs > 0.0 { untimed / real_wall_secs * 100.0 } else { 0.0 })?;
    Ok(())
}

fn print_summary(result: &ProcessResult, output_db: &Path) {
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!("  Samples processed: {}/{}", result.samples_processed, result.samples_processed + result.samples_failed + result.samples_no_coverage);
    if result.samples_no_coverage > 0 {
        eprintln!("  Samples with no retained coverage: {} (recorded in DB after filtering)", result.samples_no_coverage);
    }
    eprintln!("  Total time:        {:.2}s", result.total_time_secs);
    if let Ok(meta) = std::fs::metadata(output_db) {
        eprintln!("  Total DB size:     {:.1} MB", meta.len() as f64 / (1024.0 * 1024.0));
    }
    eprintln!();
    eprintln!("  Output: {:?}", output_db);
}

#[cfg(test)]
mod mag_processing_tests {
    use super::*;
    use rust_htslib::bam::{self, header::HeaderRecord};

    fn config(min_af: f64, min_cov: f64) -> ProcessConfig {
        ProcessConfig {
            threads: 1, min_aligned_fraction: min_af, min_coverage_depth: min_cov,
            bar_ratio: 10.0, sequencing_type: Some(SequencingType::ShortSingle),
            phagetermini_config: PhageTerminiConfig::default(), gc_params: GCParams::default(),
            variation_percentage: 0.0, min_occurrences: 2, enable_timing: true,
            blast: BlastMode::None, sample_name: String::new(), primary_reads_seen: Arc::new(AtomicU64::new(0)), terminal_repeats: None, view_mode: ViewMode::Mag,
            mag_manifest: vec![], contig_drops_counter: Arc::new(AtomicUsize::new(0)),
        }
    }

    fn fixture() -> (tempfile::TempDir, PathBuf, MagJobGroup) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let mut header = bam::Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
        for name in ["a", "b", "absent"] {
            header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", name).push_tag(b"LN", 100));
        }
        let view = bam::HeaderView::from_header(&header);
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        for line in [
            "r1\t0\ta\t1\t60\t10M5N5M\t*\t0\t0\tAAAAAAAAAAAAAAA\t*\tMD:Z:15",
            "sec\t256\ta\t1\t60\t50M\t*\t0\t0\t*\t*\tMD:Z:50",
            "supp\t2048\ta\t1\t60\t50M\t*\t0\t0\t*\t*\tMD:Z:50",
            "r2\t0\tb\t1\t60\t10M2D8M\t*\t0\t0\tAAAAAAAAAAAAAAAAAA\t*\tMD:Z:10^AA8",
        ] {
            writer.write(&bam::Record::from_sam(&view, line.as_bytes()).unwrap()).unwrap();
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        let group = MagJobGroup {
            mag_name: "mag".into(), mag_id: 1, mag_length: 300,
            // Absent member has no job, but contributes to the full MAG denominator.
            members: ["a", "b"].iter().enumerate().map(|(tid, name)| MagContigJob {
                tid: tid as u32, contig_id: tid as i64 + 1, ref_name: name.to_string(),
                ref_length: 100, offset_in_mag: tid as u32 * 100,
            }).collect(),
        };
        (directory, path, group)
    }

    #[test]
    fn whole_sample_primary_counts_survive_coverage_filtering() {
        let (_directory,path,_group)=fixture();
        let mut configuration=config(99.0,0.0);configuration.view_mode=ViewMode::Contig;
        let names=HashMap::from([("a",0usize),("b",1usize),("absent",2usize)]);
        let result=process_sample(&path,ModuleFlags::from_modules(&["Coverage".into()]),&configuration,false,&[],&names,&HashMap::new(),None).unwrap();
        assert!(result.presences.is_empty());
        assert_eq!((result.total_reads,result.mapped_reads),(2,2));
        assert_eq!(result.bai_mapped_reads,4);
    }

    #[test]
    fn database_lookup_resolves_actual_ids_independent_of_vector_order() {
        let directory=tempfile::tempdir().unwrap();let path=directory.path().join("ids.db");
        let contigs=vec![ContigInfo{name:"a".into(),aliases:vec![],length:100,sequence:None},ContigInfo{name:"b".into(),aliases:vec![],length:100,sequence:None}];
        let writer=DbWriter::create(&path,&contigs,&[],&[],true,false).unwrap();
        writer.lock_conn().unwrap().execute("UPDATE Contig SET Contig_id=7 WHERE Contig_name='b'",[]).unwrap();drop(writer);
        let writer=DbWriter::open(&path,&[],&[],&[],true).unwrap();
        let reordered=vec![contigs[1].clone(),contigs[0].clone()];
        let lookup=database_contig_lookup(&reordered,&writer).unwrap();
        assert_eq!(lookup.find("a"),Some((1,1)));assert_eq!(lookup.find("b"),Some((0,7)));
    }

    #[test]
    fn accumulator_matches_indexed_reader_when_records_cross_contig_boundaries() {
        use rust_htslib::bam::Read;
        let (_directory, path, group) = fixture();
        let flags = ModuleFlags::from_modules(&["Coverage".into(), "Misalignment".into(), "RNA".into()]);
        let mut sequential = bam::Reader::from_path(&path).unwrap();
        let records = sequential.records().collect::<std::result::Result<Vec<_>, _>>().unwrap();
        let mut indexed = IndexedReader::from_path(&path).unwrap();
        let (mut cigar, mut sequence) = (Vec::new(), Vec::new());
        for tid in 0..3 {
            let mut accumulator = crate::bam_reader::ContigAccumulator::new(100,
                SequencingType::ShortSingle, flags, false, 5, false);
            for record in records.iter().filter(|r| r.tid() == tid) {
                accumulator.consume_record(record, &mut cigar, &mut sequence);
            }
            let actual = accumulator.finish(0.0, &mut None);
            let expected = process_contig_streaming(&mut indexed, tid as u32, "fixture", 100,
                SequencingType::ShortSingle, flags, false, 0.0, 5, &mut None).unwrap();
            match (actual, expected) {
                (Some((mut actual, actual_af, actual_reads)), Some((mut expected, expected_af, expected_reads))) => {
                    assert_eq!((actual_af, actual_reads), (expected_af, expected_reads));
                    let (mut actual_blobs, mut expected_blobs) = (Vec::new(), Vec::new());
                    for (arrays, output) in [(&mut actual, &mut actual_blobs), (&mut expected, &mut expected_blobs)] {
                        add_features_from_arrays(arrays, &group.members[tid as usize].ref_name,
                            tid as i64 + 1, &config(0.0, 0.0), false, SequencingType::ShortSingle,
                            flags, actual_reads, &[], None, output, None);
                    }
                    assert_eq!(actual_blobs.len(), expected_blobs.len());
                    for (a, e) in actual_blobs.iter().zip(&expected_blobs) {
                        assert_eq!((&a.0, a.1, &a.2.zoom, &a.2.chunks), (&e.0, e.1, &e.2.zoom, &e.2.chunks));
                    }
                    crate::features::recycle_feature_arrays(actual);
                    crate::features::recycle_feature_arrays(expected);
                }
                (None, None) => {}
                _ => panic!("Accumulator and indexed reader disagree for reference {tid}"),
            }
        }
    }

    #[test]
    fn selective_recycled_arrays_preserve_encoded_features() {
        for modules in [
            vec!["Coverage"], vec!["Coverage", "Misalignment"], vec!["Coverage", "RNA"],
            vec!["Coverage", "Phage termini"], vec!["Coverage", "Paired-reads"], vec!["Coverage", "Long-reads"],
        ] {
            let flags = ModuleFlags::from_modules(&modules.iter().map(|s| s.to_string()).collect::<Vec<_>>());
            let seq_type = if flags.long_read_metrics { SequencingType::Long }
                else if flags.paired_read_metrics { SequencingType::ShortPaired }
                else { SequencingType::ShortSingle };
            let mut selected = FeatureArrays::for_modules(100, flags);
            for len in [100, 120, 80] {
                selected.reset(len, flags);
                let mut full = FeatureArrays::new(len);
                let cigar = vec![(b'S' as u32, 5), (b'M' as u32, 10), (b'I' as u32, 2),
                    (b'M' as u32, 3), (b'D' as u32, 2), (b'M' as u32, 5), (b'S' as u32, 4)];
                for arrays in [&mut full, &mut selected] {
                    for _ in 0..2 {
                        crate::features::process_read(arrays, 3, 23, 29, 20, true, true, false,
                            false, false, false, false, &cigar, Some(b"13^AA5"), &vec![b'A'; 29],
                            60, seq_type, flags, false, 5);
                    }
                    if flags.phagetermini { arrays.finalize_strands(seq_type); }
                }
                assert_eq!(full.primary_reads, selected.primary_reads);
                assert_eq!(selected.coverage_median_and_trimmed_mean(0.05),
                    (full.coverage_median(), full.coverage_trimmed_mean(0.05)));
                let mut expected = Vec::new();
                let mut actual = Vec::new();
                add_features_from_arrays(&mut full, "a", 1, &config(0.0, 0.0), false,
                    seq_type, flags, 2, &[], None, &mut expected, None);
                add_features_from_arrays(&mut selected, "a", 1, &config(0.0, 0.0), false,
                    seq_type, flags, 2, &[], None, &mut actual, None);
                assert_eq!(actual.len(), expected.len());
                for (actual, expected) in actual.iter().zip(&expected) {
                    assert_eq!((&actual.0, actual.1, &actual.2.zoom, &actual.2.chunks),
                        (&expected.0, expected.1, &expected.2.zoom, &expected.2.chunks));
                }
                if !flags.phagetermini { assert!(selected.reads_starts.is_empty()); }
                if !flags.paired_read_metrics { assert!(selected.sum_insert_sizes.is_empty()); }
            }
        }
    }

    #[test]
    fn failed_mag_sink_is_reported_to_the_sample_caller() {
        let (_directory, path, _group) = fixture();
        let flags = ModuleFlags::from_modules(&["Coverage".into()]);
        let names = HashMap::from([("a", 0), ("b", 1), ("absent", 2)]);
        let lookup = HashMap::from([
            ("a".into(), ("mag".into(), 1, 0, 300)),
            ("b".into(), ("mag".into(), 1, 100, 300)),
        ]);
        let sink = |message| match message {
            MagStreamMessage::Start { .. } => Ok(()),
            MagStreamMessage::Group(_) => anyhow::bail!("injected writer failure"),
        };
        let result = process_sample_with_mag_sink(&path, flags, &config(0.0, 0.0), false,
            &[], &names, &HashMap::new(), Some(&lookup), Some(&sink));
        assert!(result.err().unwrap().to_string().contains("injected writer failure"));
    }

    #[test]
    fn streamed_mag_database_matches_collected_database() {
        for threshold in [0.0, 11.7, 20.0] {
            let (directory, path, _group) = fixture();
            let contigs: Vec<_> = ["a", "b", "absent"].iter().map(|name| ContigInfo {
                name: name.to_string(), aliases: vec![], length: 100, sequence: None,
            }).collect();
            let baseline = DbWriter::create(&directory.path().join("baseline.db"), &contigs, &[], &[], true, true).unwrap();
            let streamed = DbWriter::create(&directory.path().join("streamed.db"), &contigs, &[], &[], true, true).unwrap();
            let mags = vec![("mag".into(), vec!["a".into(), "b".into(), "absent".into()])];
            baseline.write_mags(&mags).unwrap();
            streamed.write_mags(&mags).unwrap();
            let lookup = crate::mag_blob::build_mag_member_lookup(&baseline).unwrap();
            let names = HashMap::from([("a", 0), ("b", 1), ("absent", 2)]);
            let flags = ModuleFlags::from_modules(&["Coverage".into(), "Misalignment".into()]);
            let cfg = config(threshold, 0.0);
            let expected = process_sample(&path, flags, &cfg, false, &[], &names, &HashMap::new(), Some(&lookup)).unwrap();
            let (drop, _) = filter_failing_mags(&expected.presences, &mags, &expected.mag_results, threshold, 0.0);
            let data = apply_mag_filter(&drop, &mags, expected.presences, expected.packaging,
                expected.misassembly, expected.microdiversity, expected.side_misassembly, expected.topology,
                expected.feature_blobs, expected.mag_results);
            let (total,mapped)=(expected.total_reads,expected.mapped_reads);
            let id = baseline.insert_sample(&expected.sample_name, expected.seq_type.as_str(), total, mapped, false).unwrap();
            if !data.0.is_empty() {
                baseline.write_sample_data(&expected.sample_name, &data.0, &data.1, &data.2, &data.3,
                    &data.4, &data.5, &data.6, false).unwrap();
                let conn = baseline.lock_conn().unwrap();
                baseline.write_mag_coverages(&conn, id, &data.7).unwrap();
            }
            let (actual, writes) = process_mag_sample_streamed(&path, flags, &cfg, false, &[],
                &names, &HashMap::new(), &lookup, &streamed).unwrap();
            assert!(actual.presences.is_empty() && actual.feature_blobs.is_empty());
            assert_eq!(writes.coverage_rows, data.0.len());
            for table in ["Sample", "Coverage", "Misassembly", "Microdiversity", "Side_misassembly",
                "Topology", "Feature_blob", "Feature_blob_chunk", "MAG_coverage"] {
                let rows = |writer: &DbWriter| {
                    let conn = writer.lock_conn().unwrap();
                    let mut stmt = conn.prepare(&format!("SELECT CAST(to_json(t) AS VARCHAR) AS row FROM {table} t ORDER BY row")).unwrap();
                    stmt.query_map([], |row| row.get::<_, String>(0)).unwrap()
                        .collect::<duckdb::Result<Vec<_>>>().unwrap()
                };
                assert_eq!(rows(&baseline), rows(&streamed), "table {table}, threshold {threshold}");
            }
        }
    }

    #[test]
    fn bam_order_does_not_change_mag_coordinate_offsets() {
        let (_directory, path, mut group) = fixture();
        group.members[0].offset_in_mag = 100;
        group.members[1].offset_in_mag = 0;
        group.members.sort_by_key(|m| m.offset_in_mag);
        assert_eq!(group.members[0].tid, 1);
        let order = mag_members_in_bam_order(&group);
        assert_eq!((order[0].tid, order[0].offset_in_mag), (0, 100));
        let mut bam = IndexedReader::from_path(&path).unwrap();
        let flags = ModuleFlags::from_modules(&["Coverage".into()]);
        let result = process_mag_group(&mut bam, &group, SequencingType::ShortSingle, flags,
            &config(0.0, 0.0), false, &[], &HashMap::new(), None).unwrap().unwrap();
        assert_eq!(result.mag_cov_stats.aligned_fraction_pct_x10, 117);
        assert_eq!(result.mag_cov_stats.coverage_mean_x10, 1);
    }

    #[test]
    fn mag_thresholds_filter_completed_results_at_scaled_boundary() {
        let (_directory, path, group) = fixture();
        let flags = ModuleFlags::from_modules(&["Coverage".into(), "Misalignment".into()]);
        let mut bam = IndexedReader::from_path(&path).unwrap();
        let baseline = process_mag_group(&mut bam, &group, SequencingType::ShortSingle, flags,
            &config(0.0, 0.0), false, &[], &HashMap::new(), None).unwrap().unwrap();
        assert_eq!(baseline.mag_cov_stats.read_count, 2);
        assert_eq!(baseline.mag_cov_stats.aligned_fraction_pct_x10, 117); // 35/300 bp -> 11.7%
        assert_eq!(baseline.mag_cov_stats.coverage_mean_x10, 1); // 35/300 -> 0.1x after rounding
        assert!(!baseline.contig_blobs.is_empty());
        for (af, cov, rejected) in [(11.7, 0.1, false), (11.71, 0.0, true), (0.0, 0.11, true)] {
            let result = process_mag_group(&mut bam, &group, SequencingType::ShortSingle, flags,
                &config(af, cov), false, &[], &HashMap::new(), None).unwrap().unwrap();
            // Thresholds no longer bypass full processing or encoding.
            assert_eq!(result.mapped_reads, baseline.mapped_reads);
            assert_eq!(result.misassembly.len(), baseline.misassembly.len());
            assert_eq!(result.microdiversity.len(), baseline.microdiversity.len());
            assert_eq!(result.presences.len(), baseline.presences.len());
            assert_eq!(result.contig_blobs.len(), baseline.contig_blobs.len());
            for (actual, expected) in result.contig_blobs.iter().zip(&baseline.contig_blobs) {
                assert_eq!((&actual.0, actual.1, &actual.2.zoom, &actual.2.chunks),
                    (&expected.0, expected.1, &expected.2.zoom, &expected.2.chunks));
            }
            let (drop, count) = filter_failing_mags(&result.presences,
                &[("mag".into(), vec!["a".into(), "b".into(), "absent".into()])],
                &vec![(result.mag_name, result.mag_id, result.mag_cov_stats)], af, cov);
            assert_eq!(count, usize::from(rejected));
            assert_eq!(drop.len(), if rejected { 3 } else { 0 });
        }
    }


}

#[cfg(test)]
mod reader_budget_tests {
    use super::bam_reader_worker_limit;
    #[test]
    fn large_headers_and_low_headroom_bound_reader_count() {
        assert_eq!(bam_reader_worker_limit(20, 7_228_544, 3530, 273 * 1024, 62_589), 20);
        assert_eq!(bam_reader_worker_limit(20, 10_631_709, 4000, 273 * 1024, 72 * 1024), 13);
        assert_eq!(bam_reader_worker_limit(20, 10_631_709, 4000, 136 * 1024, 72 * 1024), 3);
        assert_eq!(bam_reader_worker_limit(20, 10_631_709, 4000, 136 * 1024, 100 * 1024), 1);
        assert_eq!(bam_reader_worker_limit(20, 100, 64, 136 * 1024, 1024), 20);
        assert_eq!(bam_reader_worker_limit(1, 10_631_709, 4000, 136 * 1024, 72 * 1024), 1);
        assert_eq!(bam_reader_worker_limit(20, 100, 64, 1024, 2048), 1);
        assert_eq!(bam_reader_worker_limit(0, 100, 64, 136 * 1024, 1024), 1);
    }
}
