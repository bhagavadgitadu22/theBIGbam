//! Parallel BAM file processing and feature calculation.
//!
//! Orchestrates parallel processing of BAM files:
//! 1. Each BAM file processed by separate thread (rayon)
//! 2. Results collected in memory
//! 3. Written to DuckDB sequentially (DuckDB is single-writer)
//!
//! BAM processing (95% of time) is fully parallelized.
//! Database writing (5% of time) runs sequentially after processing completes.

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

use crate::bam_reader::{detect_sequencing_type, process_contig_streaming, BamPhaseTimings};
use crate::compress::{
    Run,
    add_compressed_feature_with_stats,
};
use crate::db::{DbWriter, MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData, GCContentData, RepeatsData};
use crate::gc_content::{compute_gc_content, compute_gc_skew, GCParams};
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

/// Result of processing one MAG group (all member contigs + MAG-level blobs).
struct MagGroupResult {
    mag_name: String,
    mag_id: i64,
    contig_blobs: Vec<(String, i64, crate::blob::EncodedBlob)>,
    mag_blobs: Vec<(String, crate::blob::EncodedBlob)>,
    mag_cov_stats: crate::mag_blob::MagCoverageStats,
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    mapped_reads: u64,
}

/// Contig→MAG membership mapping, built once per run from mag_contig_map + DB IDs.
/// Maps contig_name → (mag_name, mag_id, offset_in_mag, mag_length).
pub type MagMemberLookup = HashMap<String, (String, i64, u32, u32)>;

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
    /// alongside the output database. Hidden CLI flag (--time).
    pub enable_timing: bool,
    /// When true, run BLAST-based repeat detection (autoblast in contig mode,
    /// per-MAG blast in MAG mode). Off by default.
    pub blast: bool,
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
    pub build_lookups: std::time::Duration,
    pub filter_tids: std::time::Duration,
    pub par_iter_wall: std::time::Duration,
    pub merge_results: std::time::Duration,
    // Per-contig sub-phase sums inside par_iter (nanoseconds, summed across threads).
    // bam_open_ns now counts map_init calls (~once per rayon worker chunk),
    // not once per contig.
    pub bam_open_ns: u64,
    pub cds_index_ns: u64,
    pub streaming_ns: u64,
    pub feature_calc_ns: u64,
    // Breakdown of streaming_ns (see BamPhaseTimings in bam_reader.rs).
    pub stream_fetch_seek_ns: u64,
    pub stream_records_ns: u64,
    pub stream_process_read_ns: u64,
    pub stream_finalize_ns: u64,
    // MAG blob encoding inside par_iter (nanoseconds, summed across threads).
    pub mag_encoding_ns: u64,
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
    cds_index_ns: AtomicU64,
    streaming_ns: AtomicU64,
    feature_calc_ns: AtomicU64,
    stream_fetch_seek_ns: AtomicU64,
    stream_records_ns: AtomicU64,
    stream_process_read_ns: AtomicU64,
    stream_finalize_ns: AtomicU64,
    mag_encoding_ns: AtomicU64,
}

impl TimingAccum {
    fn new() -> Self {
        Self {
            bam_open_ns: AtomicU64::new(0),
            cds_index_ns: AtomicU64::new(0),
            streaming_ns: AtomicU64::new(0),
            feature_calc_ns: AtomicU64::new(0),
            stream_fetch_seek_ns: AtomicU64::new(0),
            stream_records_ns: AtomicU64::new(0),
            stream_process_read_ns: AtomicU64::new(0),
            stream_finalize_ns: AtomicU64::new(0),
            mag_encoding_ns: AtomicU64::new(0),
        }
    }
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
    _contigs: &[ContigInfo],
    sample_name: &str,
) -> Result<bool> {
    let header = bam.header();

    // 1. Check @CO header for theBIGbam:circular=true/false
    let header_text = String::from_utf8_lossy(header.as_bytes());
    for line in header_text.lines() {
        if let Some(rest) = line.strip_prefix("@CO\t") {
            if let Some(val) = rest.strip_prefix("theBIGbam:circular=") {
                match val.trim() {
                    "true" => {
                        eprintln!("  Sample '{}': circular=true (from BAM header)", sample_name);
                        return Ok(true);
                    }
                    "false" => {
                        eprintln!("  Sample '{}': circular=false (from BAM header)", sample_name);
                        return Ok(false);
                    }
                    _ => {} // unrecognized value, fall through to length check
                }
            }
        }
    }

    // 2. No header tag found — assume linear
    eprintln!("  Sample '{}': circular=false (no circularity tag in BAM header)", sample_name);
    Ok(false)
}

/// Strict upfront validation of all inputs.
/// Collects all errors across all samples, prints them, and returns Err if any found.
/// This runs as the very first step before any parsing or database creation.
fn validate_inputs(
    bam_files: &[PathBuf],
    modules: &[String],
    _genbank_path: &Path,
    _assembly_path: &Path,
    _extend_db: &Path,
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

/// Quick scan of @CO headers only (no logging). Used for BAM-only contig extraction
/// before the full detection runs.
fn quick_co_header_scan(bam_files: &[PathBuf]) -> Result<HashMap<PathBuf, bool>> {
    let mut result = HashMap::new();
    for bam_path in bam_files {
        let bam = IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;
        let header_text = String::from_utf8_lossy(bam.header().as_bytes());
        let mut circular = false;
        for line in header_text.lines() {
            if let Some(rest) = line.strip_prefix("@CO\t") {
                if let Some(val) = rest.strip_prefix("theBIGbam:circular=") {
                    circular = val.trim() == "true";
                    break;
                }
            }
        }
        result.insert(bam_path.clone(), circular);
    }
    Ok(result)
}

/// Detect circularity for all samples upfront. Returns a map of BAM path → is_circular.
fn detect_all_sample_circularities(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
) -> Result<HashMap<PathBuf, bool>> {
    eprintln!("\n### Auto-detecting circularity per sample...");
    let mut result = HashMap::new();
    for bam_path in bam_files {
        let sample_name = bam_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .replace("_with_MD", "");

        let bam = IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM for circularity detection: {}", bam_path.display()))?;

        let is_circular = detect_sample_circularity(&bam, contigs, &sample_name)?;
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
    let mut matched = 0usize;
    let mut added = 0usize;
    for (name, seq) in &records {
        if let Some(contig) = contigs.iter_mut().find(|c| c.name == *name) {
            contig.sequence = Some(seq.clone());
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
fn run_autoblast(contigs: &[ContigInfo], _threads: usize) -> Result<Vec<RepeatsData>> {
    let contigs_with_seq: Vec<&ContigInfo> = contigs.iter().filter(|c| c.sequence.is_some()).collect();
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
    eprintln!("\n### Running repeat detection (BLAST) on {} contigs...", total_contigs);

    let done = std::sync::atomic::AtomicUsize::new(0);

    let all_repeats: Vec<Vec<RepeatsData>> = contigs_with_seq
        .par_iter()
        .filter_map(|contig| {
            let temp_file = match write_contigs_fasta(&[*contig]) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Warning: failed to write FASTA for {}: {}", contig.name, e);
                    return None;
                }
            };
            let temp_path = temp_file.path().to_path_buf();

            let output = match std::process::Command::new("blastn")
                .arg("-query").arg(&temp_path)
                .arg("-subject").arg(&temp_path)
                .arg("-evalue").arg("1e-10")
                .arg("-outfmt").arg("6")
                .output()
            {
                Ok(o) => o,
                Err(e) => {
                    eprintln!("Warning: blastn failed for {}: {}", contig.name, e);
                    return None;
                }
            };
            if !output.status.success() {
                return None;
            }

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
                eprintln!("  BLAST: {}/{} contigs", n, total_contigs);
            }
            Some(repeats)
        })
        .collect();

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
) -> Result<(Vec<RepeatsData>, usize)> {
    for tool in &["blastn", "makeblastdb"] {
        match std::process::Command::new(tool).arg("-version").output() {
            Ok(output) if output.status.success() => {}
            _ => {
                eprintln!("Warning: {} not found on PATH. Skipping MAG inter-contig BLAST.", tool);
                return Ok((Vec::new(), 0));
            }
        }
    }

    let name_to_contig: HashMap<&str, &ContigInfo> =
        contigs.iter().map(|c| (c.name.as_str(), c)).collect();

    eprintln!("\n### Running MAG inter-contig BLAST on {} MAGs...", mag_manifest.len());

    let total_repeat_count = std::sync::atomic::AtomicUsize::new(0);
    let total_inter_count = std::sync::atomic::AtomicUsize::new(0);
    let mags_done = std::sync::atomic::AtomicUsize::new(0);
    let total_mags = mag_manifest.len();

    let all_repeats: std::sync::Mutex<Vec<RepeatsData>> = std::sync::Mutex::new(Vec::new());

    mag_manifest.par_iter().for_each(|(mag_name, member_names)| {
        let members: Vec<&ContigInfo> = member_names
            .iter()
            .filter_map(|n| name_to_contig.get(n.as_str()).copied())
            .filter(|c| c.sequence.is_some())
            .collect();
        if members.is_empty() {
            mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return;
        }

        let mut temp_file = match tempfile::NamedTempFile::new() {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Warning: temp file for MAG '{}': {}", mag_name, e);
                mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return;
            }
        };
        use std::io::Write;
        for contig in &members {
            if writeln!(temp_file, ">{}", contig.name).is_err() { return; }
            if let Some(seq) = contig.sequence.as_ref() {
                for chunk in seq.chunks(80) {
                    if temp_file.write_all(chunk).is_err() { return; }
                    if writeln!(temp_file).is_err() { return; }
                }
            }
        }
        if temp_file.flush().is_err() { return; }
        let temp_path = temp_file.path().to_path_buf();

        let mkdb = std::process::Command::new("makeblastdb")
            .arg("-in").arg(&temp_path)
            .arg("-dbtype").arg("nucl")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status();
        match mkdb {
            Ok(s) if s.success() => {}
            Ok(s) => {
                eprintln!("Warning: makeblastdb failed for MAG '{}' (exit {})", mag_name, s);
                mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return;
            }
            Err(e) => {
                eprintln!("Warning: makeblastdb failed for MAG '{}': {}", mag_name, e);
                mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return;
            }
        }

        let output = match std::process::Command::new("blastn")
            .arg("-query").arg(&temp_path)
            .arg("-db").arg(&temp_path)
            .arg("-evalue").arg("1e-10")
            .arg("-outfmt").arg("6")
            .output()
        {
            Ok(o) => o,
            Err(e) => {
                eprintln!("Warning: blastn failed for MAG '{}': {}", mag_name, e);
                mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return;
            }
        };
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            eprintln!("Warning: blastn non-zero for MAG '{}': {}", mag_name, stderr);
            mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return;
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        let mut repeats: Vec<RepeatsData> = Vec::new();
        let mut inter: Vec<(String, String, i32, i32, i32, i32, f64)> = Vec::new();
        for line in stdout.lines() {
            let line = line.trim();
            if line.is_empty() { continue; }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 10 { continue; }
            let qid = f[0].to_string();
            let sid = f[1].to_string();
            let pident: f64 = f[2].parse().unwrap_or(0.0);
            let qstart: i32 = f[6].parse().unwrap_or(0);
            let qend:   i32 = f[7].parse().unwrap_or(0);
            let sstart: i32 = f[8].parse().unwrap_or(0);
            let send:   i32 = f[9].parse().unwrap_or(0);

            if qid == sid {
                if qstart == sstart && qend == send { continue; }
                let is_direct = (qstart < qend && sstart < send) || (qstart > qend && sstart > send);
                repeats.push(RepeatsData {
                    contig_name: qid,
                    position1: qstart, position2: qend,
                    position1prime: sstart, position2prime: send,
                    pident, is_direct,
                });
            } else {
                inter.push((qid, sid, qstart, qend, sstart, send, pident));
            }
        }

        let n_repeats = repeats.len();
        let n_inter = inter.len();

        if let Ok(conn) = db_writer.lock_conn() {
            use duckdb::params;
            for dup in &repeats {
                if let Some(contig_id) = db_writer.get_contig_id(&dup.contig_name) {
                    let pident_int = (dup.pident * 100.0).round() as i32;
                    let table = if dup.is_direct { "Contig_directRepeats" } else { "Contig_invertedRepeats" };
                    let _ = conn.execute(
                        &format!("INSERT INTO {} (Contig_id, Position1, Position2, Position1prime, Position2prime, Pident) VALUES (?, ?, ?, ?, ?, ?)", table),
                        params![contig_id, dup.position1, dup.position2, dup.position1prime, dup.position2prime, pident_int],
                    );
                }
            }
            if !inter.is_empty() {
                let mut stmt = conn.prepare(
                    "INSERT OR IGNORE INTO Contig_blast_hits
                       (Contig_id_1, Contig_id_2, Position1, Position2, Position1prime, Position2prime, Pident)
                     VALUES (?, ?, ?, ?, ?, ?, ?)"
                ).ok();
                if let Some(ref mut stmt) = stmt {
                    for (qid, sid, mut p1, mut p2, mut p1p, mut p2p, pident) in inter.iter().cloned() {
                        let cid_q = match db_writer.get_contig_id(&qid) { Some(v) => v, None => continue };
                        let cid_s = match db_writer.get_contig_id(&sid) { Some(v) => v, None => continue };
                        let (cid1, cid2) = if cid_q <= cid_s {
                            (cid_q, cid_s)
                        } else {
                            std::mem::swap(&mut p1, &mut p1p);
                            std::mem::swap(&mut p2, &mut p2p);
                            (cid_s, cid_q)
                        };
                        let pident_i = (pident * 100.0).round() as i32;
                        let _ = stmt.execute(params![cid1, cid2, p1, p2, p1p, p2p, pident_i]);
                    }
                }
            }
        }

        if let Some(p) = temp_path.to_str() {
            for ext in &[".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"] {
                let _ = std::fs::remove_file(format!("{}{}", p, ext));
            }
        }

        all_repeats.lock().unwrap().extend(repeats);

        total_repeat_count.fetch_add(n_repeats, std::sync::atomic::Ordering::Relaxed);
        total_inter_count.fetch_add(n_inter, std::sync::atomic::Ordering::Relaxed);
        let done = mags_done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        if done % 100 == 0 || done == total_mags {
            eprintln!("  MAG BLAST: {}/{} MAGs", done, total_mags);
        }
    });

    let repeat_count = total_repeat_count.load(std::sync::atomic::Ordering::Relaxed);
    let inter_count = total_inter_count.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!(
        "MAG BLAST: {} intra-contig repeat hits, {} inter-contig hits",
        repeat_count, inter_count,
    );
    let reps = all_repeats.into_inner().unwrap();
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
) -> (Option<PackagingData>, Option<(MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData)>) {
    use crate::blob::{encode_dense_blob, smooth_dense_values, encode_sparse_blob, EventMeta, MetadataFlags,
                       codon_category_to_id, codon_to_id, aa_to_id};
    use crate::types::{get_encoding, get_value_scale, Encoding};
    let pt_config = config.phagetermini_config;
    let contig_length = arrays.ref_length();

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
    if flags.coverage {
        // MAPQ - average mapping quality per position (needed for blob encoding below)
        let mapq_f64: Vec<f64> = arrays.sum_mapq.iter()
            .zip(&arrays.primary_reads)
            .map(|(&sum, &count)| if count > 0 { sum as f64 / count as f64 } else { 0.0 })
            .collect();

        // === BLOB encoding for dense coverage features ===
        let clen = contig_length as u32;

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

        // MAPQ: stored as ×100 integer
        let mut mapq_i32: Vec<i32> = mapq_f64.iter().map(|&x| (x * 100.0).round() as i32).collect();
        smooth_dense_values(&mut mapq_i32, config.variation_percentage);
        debug_assert_eq!(get_encoding("mapq"), Encoding::Dense);
        blob_output.push(("mapq".into(), cid, encode_dense_blob(&mapq_i32, get_value_scale("mapq"), clen)));
    }

    // Assemblycheck features
    let mut left_clip_runs: Vec<Run> = Vec::new();
    let mut right_clip_runs: Vec<Run> = Vec::new();
    if flags.mapping_metrics || flags.phagetermini {
        // Clippings and insertions with statistics
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
        left_clip_runs = add_compressed_feature_with_stats(&left_clip_counts, &left_clip_means, &left_clip_medians, &left_clip_stds,
            Some(&primary_reads_f64), "left_clippings", config);

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
        right_clip_runs = add_compressed_feature_with_stats(&right_clip_counts, &right_clip_means, &right_clip_medians, &right_clip_stds,
            Some(&primary_reads_f64), "right_clippings", config);
    }
            
    // Shared constants for sparse feature encoding (used by mapping_metrics + rna modules)
    let clen = contig_length as u32;
    let bar_threshold = config.bar_ratio * 0.01;
    let min_occ = config.min_occurrences;

    // Helper: filter positions by coverage threshold AND min_occurrences, produce (positions, values)
    let filter_sparse = |counts: &[f64], coverage: &[f64]| -> (Vec<u32>, Vec<i32>) {
        let mut positions = Vec::new();
        let mut values = Vec::new();
        let n = counts.len().min(coverage.len());
        for i in 0..n {
            let val = counts[i];
            let cov = coverage[i];
            if val >= cov * bar_threshold && val >= min_occ as f64 {
                positions.push(i as u32);
                values.push((val / cov * 1000.0).round() as i32);
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
        let threshold = config.bar_ratio * 0.01; // Convert percentage to fraction
        let dominant_mismatches = arrays.compute_dominant_mismatch_bases(&arrays.primary_reads, threshold);
        let dominant_insertions = FeatureArrays::compute_dominant_sequences(&arrays.insertion_sequences, &arrays.primary_reads, threshold);
        let dominant_left_clips = FeatureArrays::compute_dominant_sequences(&arrays.left_clip_sequences, &arrays.primary_reads, threshold);
        let dominant_right_clips = FeatureArrays::compute_dominant_sequences(&arrays.right_clip_sequences, &arrays.primary_reads, threshold);

        // Compute codon changes from position-level mismatch summaries
        let dominant_codons = compute_codon_changes_from_summaries(
            &arrays.mismatch_base_counts, &arrays.primary_reads,
            cds_index, contig_length, threshold,
        );

        // mismatches (with sequence + codons)
        {
            let (pos, vals) = filter_sparse(&mismatches_f64, &primary_reads_f64);
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
            let (pos, vals) = filter_sparse(&deletions_f64, &primary_reads_f64);
            let flags = MetadataFlags { sparse: true, ..Default::default() };
            debug_assert_eq!(get_encoding("deletions"), Encoding::Sparse);
            blob_output.push(("deletions".into(), cid,
                encode_sparse_blob(&pos, &vals, None, flags, get_value_scale("deletions"), clen)));
        }

        // insertions (with stats + sequence)
        {
            let (pos, vals) = filter_sparse(&insertion_counts, &primary_reads_f64);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((insertion_means.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    median: Some((insertion_medians.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    std: Some((insertion_stds.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
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

        // left_clippings (with stats + sequence)
        {
            let lc_counts: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| v.len() as f64).collect();
            let lc_means: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
                if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
            }).collect();
            let lc_medians: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
                if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
            }).collect();
            let lc_stds: Vec<f64> = arrays.left_clipping_lengths.iter().map(|v| {
                if v.len() <= 1 { 0.0 } else {
                    let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                    let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                    var.sqrt()
                }
            }).collect();
            let (pos, vals) = filter_sparse(&lc_counts, &primary_reads_f64);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((lc_means.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    median: Some((lc_medians.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    std: Some((lc_stds.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
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

        // right_clippings (with stats + sequence)
        {
            let rc_counts: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| v.len() as f64).collect();
            let rc_means: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
                if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
            }).collect();
            let rc_medians: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
                if v.is_empty() { 0.0 } else { let mut s = v.clone(); s.sort_unstable(); s[s.len()/2] as f64 }
            }).collect();
            let rc_stds: Vec<f64> = arrays.right_clipping_lengths.iter().map(|v| {
                if v.len() <= 1 { 0.0 } else {
                    let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                    let var = v.iter().map(|&x| { let d = x as f64 - mean; d*d }).sum::<f64>() / v.len() as f64;
                    var.sqrt()
                }
            }).collect();
            let (pos, vals) = filter_sparse(&rc_counts, &primary_reads_f64);
            let meta: Vec<EventMeta> = pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    mean: Some((rc_means.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    median: Some((rc_medians.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
                    std: Some((rc_stds.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
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
        let clen = contig_length as u32;
        let bar_threshold = config.bar_ratio * 0.01;
        let min_occ = config.min_occurrences;

        // non_inward_pairs, mate_not_mapped, mate_on_another_contig: sparse, value only
        let filter_sparse_pr = |counts: &[f64], coverage: &[f64]| -> (Vec<u32>, Vec<i32>) {
            let mut positions = Vec::new();
            let mut vals = Vec::new();
            let n = counts.len().min(coverage.len());
            for i in 0..n {
                let val = counts[i];
                let cov = coverage[i];
                if val >= cov * bar_threshold && val >= min_occ as f64 {
                    positions.push(i as u32);
                    vals.push((val / cov * 1000.0).round() as i32);
                }
            }
            (positions, vals)
        };

        let (pos, vals) = filter_sparse_pr(&non_inward_f64, &primary_reads_f64);
        let sp_flags = MetadataFlags { sparse: true, ..Default::default() };
        debug_assert_eq!(get_encoding("non_inward_pairs"), Encoding::Sparse);
        blob_output.push(("non_inward_pairs".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("non_inward_pairs"), clen)));

        let (pos, vals) = filter_sparse_pr(&mate_unmapped_f64, &primary_reads_f64);
        debug_assert_eq!(get_encoding("mate_not_mapped"), Encoding::Sparse);
        blob_output.push(("mate_not_mapped".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("mate_not_mapped"), clen)));

        let (pos, vals) = filter_sparse_pr(&mate_other_contig_f64, &primary_reads_f64);
        debug_assert_eq!(get_encoding("mate_on_another_contig"), Encoding::Sparse);
        blob_output.push(("mate_on_another_contig".into(), cid,
            encode_sparse_blob(&pos, &vals, None, sp_flags, get_value_scale("mate_on_another_contig"), clen)));

        // insert_sizes: dense curve
        let mut is_i32: Vec<i32> = values.iter().map(|&x| (x * 10.0).round() as i32).collect();
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
        let clen = contig_length as u32;
        let mut rl_i32: Vec<i32> = values.iter().map(|&x| (x * 10.0).round() as i32).collect();
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
        let threshold = config.bar_ratio * 0.01;
        let dominant_start_clips = FeatureArrays::compute_dominant_sequences(&arrays.start_clip_sequences, &arrays.coverage_reduced, threshold);
        let dominant_end_clips = FeatureArrays::compute_dominant_sequences(&arrays.end_clip_sequences, &arrays.coverage_reduced, threshold);

        // === BLOB encoding for phagetermini features ===
        {
            let clen = contig_length as u32;
            let bar_threshold = config.bar_ratio * 0.01;
            let min_occ = config.min_occurrences;

            // coverage_reduced: dense curve
            let mut cr_i32: Vec<i32> = arrays.coverage_reduced.iter().map(|&x| x as i32).collect();
            smooth_dense_values(&mut cr_i32, config.variation_percentage);
            debug_assert_eq!(get_encoding("coverage_reduced"), Encoding::Dense);
            blob_output.push(("coverage_reduced".into(), cid, encode_dense_blob(&cr_i32, get_value_scale("coverage_reduced"), clen)));

            // reads_starts: sparse with median + sequence
            let cov_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
            let mut rs_pos = Vec::new();
            let mut rs_vals = Vec::new();
            for i in 0..reads_starts_original.len().min(cov_reduced_f64.len()) {
                let val = reads_starts_original[i];
                let cov = cov_reduced_f64[i];
                if val >= cov * bar_threshold && val >= min_occ as f64 {
                    rs_pos.push(i as u32);
                    rs_vals.push((val / cov * 1000.0).round() as i32);
                }
            }
            let rs_meta: Vec<EventMeta> = rs_pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    median: Some((start_evt_medians.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
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
            let mut re_pos = Vec::new();
            let mut re_vals = Vec::new();
            for i in 0..reads_ends_original.len().min(cov_reduced_f64.len()) {
                let val = reads_ends_original[i];
                let cov = cov_reduced_f64[i];
                if val >= cov * bar_threshold && val >= min_occ as f64 {
                    re_pos.push(i as u32);
                    re_vals.push((val / cov * 1000.0).round() as i32);
                }
            }
            let re_meta: Vec<EventMeta> = re_pos.iter().map(|&p| {
                let idx = p as usize;
                let mut em = EventMeta {
                    median: Some((end_evt_medians.get(idx).copied().unwrap_or(0.0) * 100.0).round() as i32),
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
            &left_clip_runs,
            &right_clip_runs,
        );
        Some((misassembly, microdiversity, side_misassembly, topology))
    } else {
        None
    };

    (packaging_result, metrics_result)
}

/// Process one MAG group: all member contigs in sequence with one BAM reader.
/// Accumulates dense features into MAG-scale arrays and encodes MAG blobs
/// directly — no decode→reencode cycle.
fn process_mag_group(
    bam: &mut IndexedReader,
    group: &MagJobGroup,
    seq_type: SequencingType,
    flags: ModuleFlags,
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    annos_by_contig: &HashMap<i64, Vec<&FeatureAnnotation>>,
    accum_ref: Option<&TimingAccum>,
) -> Option<MagGroupResult> {
    use crate::blob::{encode_dense_blob, smooth_dense_values, EncodedBlob};
    use crate::types::{get_value_scale, Encoding, VARIABLES};

    let mag_len = group.mag_length as usize;
    if mag_len == 0 || group.members.is_empty() { return None; }

    // Dense accumulators — one Vec<i32> per dense sample feature.
    // Initialized to 0; each contig writes its slice at offset.
    struct DenseAccum {
        primary_reads: Vec<i32>,
        primary_reads_plus_only: Vec<i32>,
        primary_reads_minus_only: Vec<i32>,
        secondary_reads: Vec<i32>,
        supplementary_reads: Vec<i32>,
        mapq: Vec<i32>,
        coverage_reduced: Vec<i32>,
    }
    impl DenseAccum {
        fn new(len: usize) -> Self {
            Self {
                primary_reads: vec![0i32; len],
                primary_reads_plus_only: vec![0i32; len],
                primary_reads_minus_only: vec![0i32; len],
                secondary_reads: vec![0i32; len],
                supplementary_reads: vec![0i32; len],
                mapq: vec![0i32; len],
                coverage_reduced: vec![0i32; len],
            }
        }
    }

    let mut accum = DenseAccum::new(mag_len);
    let mut contig_blobs: Vec<(String, i64, EncodedBlob)> = Vec::new();
    let mut presences = Vec::new();
    let mut all_packaging = Vec::new();
    let mut all_misassembly = Vec::new();
    let mut all_microdiversity = Vec::new();
    let mut all_side_misassembly = Vec::new();
    let mut all_topology = Vec::new();
    let mut mapped_reads: u64 = 0;
    let mut contig_lengths_ordered: Vec<u32> = Vec::new();

    for member in &group.members {
        let off = member.offset_in_mag as usize;
        let clen = member.ref_length;
        contig_lengths_ordered.push(clen as u32);

        // CDS index
        let ts = accum_ref.map(|_| std::time::Instant::now());
        let cds_index = if flags.mapping_metrics {
            match annos_by_contig.get(&member.contig_id) {
                Some(slice) if !slice.is_empty() => {
                    let idx = CdsIndex::from_contig_annotations(slice);
                    if idx.is_empty() { None } else { Some(idx) }
                }
                _ => None,
            }
        } else {
            None
        };
        if let (Some(a), Some(t)) = (accum_ref, ts) {
            a.cds_index_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }

        // Process contig streaming (no per-contig AF filtering in MAG mode)
        let ts = accum_ref.map(|_| std::time::Instant::now());
        let mut phase_t: Option<BamPhaseTimings> = accum_ref.map(|_| BamPhaseTimings::default());
        let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming(
            bam, member.tid, &member.ref_name, member.ref_length, seq_type, flags,
            is_circular, 0.0, config.phagetermini_config.min_clipping_length, &mut phase_t,
        ) {
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
                eprintln!("Error processing contig {} in MAG {}: {}", member.ref_name, group.mag_name, e);
                continue;
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

        // Accumulate dense features into MAG arrays at offset (pre-smoothing raw values)
        let end = (off + clen).min(mag_len);
        let n = end - off;
        for i in 0..n.min(arrays.primary_reads.len()) {
            accum.primary_reads[off + i] = arrays.primary_reads[i] as i32;
            accum.primary_reads_plus_only[off + i] = arrays.primary_reads_plus_only[i] as i32;
            accum.primary_reads_minus_only[off + i] = arrays.primary_reads_minus_only[i] as i32;
            accum.secondary_reads[off + i] = arrays.secondary_reads[i] as i32;
            accum.supplementary_reads[off + i] = arrays.supplementary_reads[i] as i32;
            accum.coverage_reduced[off + i] = arrays.coverage_reduced[i] as i32;
        }
        // mapq: average per position
        for i in 0..n.min(arrays.sum_mapq.len()) {
            let pr = arrays.primary_reads[i];
            accum.mapq[off + i] = if pr > 0 { (arrays.sum_mapq[i] as f64 / pr as f64 * 100.0).round() as i32 } else { 0 };
        }

        // Build per-contig presence
        let coverage_mean = arrays.coverage_mean() as f32;
        let coverage_median = arrays.coverage_median() as f32;
        let coverage_trimmed_mean = arrays.coverage_trimmed_mean(0.05) as f32;
        let above_expected = coverage_pct >= (1.0 - (-0.883 * coverage_mean as f64).exp()) * 100.0;
        let coverage_relative_coverage_roughness = if arrays.primary_reads.len() > 1 && coverage_mean > 0.0 {
            let nn = arrays.primary_reads.len();
            let mean_cov = coverage_mean as f64;
            let sq_diff: f64 = arrays.primary_reads.windows(2)
                .map(|w| { let d = w[1] as f64 - w[0] as f64; d * d })
                .sum::<f64>() / (nn - 1) as f64;
            (sq_diff.sqrt() / mean_cov * 1000000.0) as f32
        } else { 0.0 };
        let coverage_coefficient_of_variation = if !arrays.primary_reads.is_empty() && coverage_mean > 0.0 {
            let mean_cov = coverage_mean as f64;
            let nn = arrays.primary_reads.len() as f64;
            let var: f64 = arrays.primary_reads.iter()
                .map(|&x| { let d = x as f64 - mean_cov; d * d }).sum::<f64>() / nn;
            (var.sqrt() / mean_cov * 1_000_000.0) as f32
        } else { 0.0 };

        presences.push(PresenceData {
            contig_id: member.contig_id,
            contig_name: member.ref_name.clone(),
            coverage_pct: coverage_pct as f32,
            above_expected_aligned_fraction: above_expected,
            read_count: primary_count,
            coverage_relative_coverage_roughness,
            coverage_coefficient_of_variation,
            coverage_mean,
            coverage_median,
            coverage_trimmed_mean,
        });
        mapped_reads += primary_count;

        // Encode per-contig blobs (existing logic, unchanged)
        let ts = accum_ref.map(|_| std::time::Instant::now());
        let mut blob_output: Vec<(String, i64, EncodedBlob)> = Vec::new();
        let (pkg, metrics) = add_features_from_arrays(
            &mut arrays, &member.ref_name, member.contig_id, config, is_circular, seq_type, flags,
            primary_count, repeats, cds_index.as_ref(), &mut blob_output,
        );
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
        return None;
    }

    // --- Encode dense MAG blobs from accumulators ---
    // Apply smoothing per-contig slice to match current behavior
    fn smooth_per_contig(values: &mut [i32], members: &[MagContigJob], variation_pct: f64) {
        let mut off = 0usize;
        for m in members {
            let end = (off + m.ref_length).min(values.len());
            smooth_dense_values(&mut values[off..end], variation_pct);
            off = end;
        }
    }

    let ts_mag_enc = accum_ref.map(|_| std::time::Instant::now());

    let vp = config.variation_percentage;
    let mut mag_blobs: Vec<(String, EncodedBlob)> = Vec::new();

    if flags.coverage {
        smooth_per_contig(&mut accum.primary_reads, &group.members, vp);
        mag_blobs.push(("primary_reads".into(), encode_dense_blob(&accum.primary_reads, get_value_scale("primary_reads"), group.mag_length)));
        smooth_per_contig(&mut accum.primary_reads_plus_only, &group.members, vp);
        mag_blobs.push(("primary_reads_plus_only".into(), encode_dense_blob(&accum.primary_reads_plus_only, get_value_scale("primary_reads_plus_only"), group.mag_length)));
        smooth_per_contig(&mut accum.primary_reads_minus_only, &group.members, vp);
        mag_blobs.push(("primary_reads_minus_only".into(), encode_dense_blob(&accum.primary_reads_minus_only, get_value_scale("primary_reads_minus_only"), group.mag_length)));
        smooth_per_contig(&mut accum.secondary_reads, &group.members, vp);
        mag_blobs.push(("secondary_reads".into(), encode_dense_blob(&accum.secondary_reads, get_value_scale("secondary_reads"), group.mag_length)));
        smooth_per_contig(&mut accum.supplementary_reads, &group.members, vp);
        mag_blobs.push(("supplementary_reads".into(), encode_dense_blob(&accum.supplementary_reads, get_value_scale("supplementary_reads"), group.mag_length)));
        smooth_per_contig(&mut accum.mapq, &group.members, vp);
        mag_blobs.push(("mapq".into(), encode_dense_blob(&accum.mapq, get_value_scale("mapq"), group.mag_length)));
    }
    if flags.phagetermini {
        smooth_per_contig(&mut accum.coverage_reduced, &group.members, vp);
        mag_blobs.push(("coverage_reduced".into(), encode_dense_blob(&accum.coverage_reduced, get_value_scale("coverage_reduced"), group.mag_length)));
    }

    // insert_sizes / read_lengths: accumulate from per-contig blobs since they
    // use derived values (running average, not raw counts). Use existing sparse
    // aggregation path below — they're dense but small and the decode is cheap.

    // --- Aggregate sparse MAG blobs from per-contig encoded blobs ---
    // Reuses existing aggregation: decode per-contig blob → shift positions → encode.
    // Sparse data is small so the decode cost is negligible.
    for var in VARIABLES.iter() {
        if crate::mag_blob::contig_blob_config_pub(var.name).is_some() { continue; }
        let is_dense_already_handled = matches!(var.name,
            "primary_reads" | "primary_reads_plus_only" | "primary_reads_minus_only" |
            "secondary_reads" | "supplementary_reads" | "mapq" | "coverage_reduced"
        );
        if is_dense_already_handled { continue; }

        let members_and_blobs: Vec<(crate::mag_blob::MagMember, &EncodedBlob)> = group.members.iter()
            .filter_map(|m| {
                let blob = contig_blobs.iter().find(|(fname, cid, _)| fname == var.name && *cid == m.contig_id);
                blob.map(|(_, _, b)| (crate::mag_blob::MagMember { contig_id: m.contig_id, offset: m.offset_in_mag, length: m.ref_length as u32 }, b))
            })
            .collect();
        if members_and_blobs.is_empty() { continue; }

        let encoded = match var.encoding {
            Encoding::Dense => crate::mag_blob::aggregate_feature_dense_inmem(&members_and_blobs, group.mag_length, var.value_scale),
            Encoding::Sparse => crate::mag_blob::aggregate_feature_sparse_inmem(&members_and_blobs, group.mag_length, var.value_scale),
        };
        if !encoded.zoom.is_empty() {
            mag_blobs.push((var.name.to_string(), encoded));
        }
    }

    // --- MAG coverage stats from raw primary_reads accumulator ---
    // Use the pre-smoothing values for stats. We already smoothed accum.primary_reads
    // above, so rebuild from presences read_count for the total.
    let total_reads_in_mag: u64 = presences.iter().map(|p| p.read_count).sum();
    // Recompute from the smoothed primary_reads (matches current behavior where
    // stats are computed from decoded blob values which are post-smoothing).
    let mag_cov_stats = crate::mag_blob::compute_mag_coverage_stats_from_raw(
        &accum.primary_reads, group.mag_length, &contig_lengths_ordered, total_reads_in_mag,
    );

    if let (Some(a), Some(t)) = (accum_ref, ts_mag_enc) {
        a.mag_encoding_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }

    Some(MagGroupResult {
        mag_name: group.mag_name.clone(),
        mag_id: group.mag_id,
        contig_blobs,
        mag_blobs,
        mag_cov_stats,
        presences,
        packaging: all_packaging,
        misassembly: all_misassembly,
        microdiversity: all_microdiversity,
        side_misassembly: all_side_misassembly,
        topology: all_topology,
        mapped_reads,
    })
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
    pub mapped_reads: u64,
    pub is_circular: bool,
    pub timings: Option<SampleTimings>,
    /// MAG-level blobs, ready to write. Empty in contig mode.
    /// Vec of (mag_name, mag_id, feature_blobs, coverage_stats).
    pub mag_results: Vec<(String, i64, Vec<(String, crate::blob::EncodedBlob)>, crate::mag_blob::MagCoverageStats)>,
}

pub fn process_sample(
    bam_path: &Path,
    _contigs: &[ContigInfo],
    flags: ModuleFlags,
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    _annotations: &[FeatureAnnotation],
    contig_by_name: &HashMap<&str, usize>,
    annos_by_contig: &HashMap<i64, Vec<&FeatureAnnotation>>,
    mag_member_lookup: Option<&MagMemberLookup>,
) -> Result<SampleProcessResult> {
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "");

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

    // --- Phase: unmapped_read_count + per-tid mapped counts ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Get per-tid stats from the BAM index in one call (avoids a second BAM open).
    // Each entry: (tid, length, mapped_records, unmapped_records).
    let idx_stats = temp_bam.index_stats()
        .with_context(|| format!("Failed to get index stats from: {}", bam_path.display()))?;
    let unmapped_reads: u64 = idx_stats.iter().map(|(_, _, _, unmapped)| *unmapped).sum();
    let mapped_per_tid: Vec<u64> = idx_stats.iter().map(|(_, _, mapped, _)| *mapped).collect();

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
            let contig_idx = contig_by_name.get(ref_name)?;
            let ref_length = temp_bam.header().target_len(tid).unwrap_or(0) as usize;
            Some((tid, (*contig_idx + 1) as i64, ref_name.to_string(), ref_length))
        })
        .collect();
    if jobs.is_empty() && n_refs > 0 {
        eprintln!("WARNING: Sample '{}': 0/{} BAM references matched known contigs. \
                   Check that contig names in your annotation/assembly files match the BAM reference names.",
                  sample_name, n_refs);
    }
    drop(temp_bam);
    if let (Some(s), Some(t)) = (st.as_mut(), t_filter) {
        s.filter_tids = t.elapsed();
    }

    // --- Phase: par_iter ---
    let accum = if time_on { Some(TimingAccum::new()) } else { None };
    let accum_ref = accum.as_ref();
    let t_par = if time_on { Some(std::time::Instant::now()) } else { None };

    // MAG mode: group jobs by MAG, process MAG-by-MAG (one BAM open per chunk).
    // Contig mode: par_chunks over individual jobs (unchanged).
    let (all_feature_blobs, all_presences, all_packaging, all_misassembly, all_microdiversity, all_side_misassembly, all_topology, mapped_reads, mag_results) = if let Some(lookup) = mag_member_lookup {
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
        // Sort members within each group by offset so smooth_per_contig works correctly
        let mut mag_groups: Vec<MagJobGroup> = mag_groups_map.into_values().collect();
        for g in &mut mag_groups {
            g.members.sort_by_key(|m| m.offset_in_mag);
        }

        let chunk_size = (mag_groups.len() + config.threads - 1).max(1) / config.threads.max(1);
        let mag_results_vec: Vec<MagGroupResult> = mag_groups
            .par_chunks(chunk_size.max(1))
            .flat_map_iter(|chunk| {
                let ts_open = accum_ref.map(|_| std::time::Instant::now());
                let mut bam = IndexedReader::from_path(bam_path)
                    .ok()
                    .and_then(|mut b| b.set_threads(1).ok().map(|_| b));
                if let (Some(a), Some(t)) = (accum_ref, ts_open) {
                    a.bam_open_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }

                chunk.iter().filter_map(move |group| {
                    let bam = bam.as_mut()?;
                    process_mag_group(bam, group, seq_type, flags, config, is_circular, repeats, annos_by_contig, accum_ref)
                })
            })
            .collect();

        if let (Some(s), Some(t)) = (st.as_mut(), t_par) {
            s.par_iter_wall = t.elapsed();
            s.n_contigs_processed = mag_results_vec.iter().map(|r| r.presences.len()).sum();
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
            mag_res.push((r.mag_name, r.mag_id, r.mag_blobs, r.mag_cov_stats));
        }

        if let (Some(s), Some(t)) = (st.as_mut(), t0) {
            s.merge_results = t.elapsed();
        }

        (all_blobs, all_pres, all_pkg, all_mis, all_micro, all_side, all_topo, mapped, mag_res)
    } else {
        // Contig mode: par_chunks over individual contigs
        let chunk_size = (jobs.len() + config.threads - 1).max(1) / config.threads.max(1);
        let results: Vec<_> = jobs
            .par_chunks(chunk_size.max(1))
            .flat_map_iter(|chunk| {
                let ts_open = accum_ref.map(|_| std::time::Instant::now());
                let mut bam = IndexedReader::from_path(bam_path)
                    .ok()
                    .and_then(|mut b| b.set_threads(1).ok().map(|_| b));
                if let (Some(a), Some(t)) = (accum_ref, ts_open) {
                    a.bam_open_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }

                chunk.iter().filter_map(move |(tid, contig_id, ref_name, ref_length)| {
                    let bam = bam.as_mut()?;
                    let tid = *tid;
                    let contig_id = *contig_id;
                    let ref_length = *ref_length;

                    let ts = accum_ref.map(|_| std::time::Instant::now());
                    let cds_index = if flags.mapping_metrics {
                        match annos_by_contig.get(&contig_id) {
                            Some(slice) if !slice.is_empty() => {
                                let idx = CdsIndex::from_contig_annotations(slice);
                                if idx.is_empty() { None } else { Some(idx) }
                            }
                            _ => None,
                        }
                    } else {
                        None
                    };
                    if let (Some(a), Some(t)) = (accum_ref, ts) {
                        a.cds_index_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    }

                    let ts = accum_ref.map(|_| std::time::Instant::now());
                    let mut phase_t: Option<BamPhaseTimings> = accum_ref.map(|_| BamPhaseTimings::default());
                    let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming(bam, tid, ref_name, ref_length, seq_type, flags, is_circular, config.min_aligned_fraction, config.phagetermini_config.min_clipping_length, &mut phase_t) {
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
                            eprintln!("Error processing contig {} in {}: {}", ref_name, bam_path.display(), e);
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
                        return None;
                    }

                    let ts = accum_ref.map(|_| std::time::Instant::now());
                    let mut feature_blobs: Vec<(String, i64, crate::blob::EncodedBlob)> = Vec::new();
                    let (packaging_info, metrics_info) = add_features_from_arrays(&mut arrays, ref_name, contig_id, config, is_circular, seq_type, flags, primary_count, repeats, cds_index.as_ref(), &mut feature_blobs);
                    let coverage_median = arrays.coverage_median() as f32;
                    let coverage_trimmed_mean = arrays.coverage_trimmed_mean(0.05) as f32;

                    let above_expected = coverage_pct >= (1.0 - (-0.883 * coverage_mean as f64).exp()) * 100.0;

                    let coverage_relative_coverage_roughness = if arrays.primary_reads.len() > 1 && coverage_mean > 0.0 {
                        let n = arrays.primary_reads.len();
                        let mean_cov = coverage_mean as f64;
                        let sum_squared_diff: f64 = arrays.primary_reads
                            .windows(2)
                            .map(|w| {
                                let diff = w[1] as f64 - w[0] as f64;
                                diff * diff
                            })
                            .sum::<f64>() / (n - 1) as f64;
                        let root_diff = sum_squared_diff.sqrt() / mean_cov;
                        (root_diff * 1000000.0) as f32
                    } else {
                        0.0
                    };

                    let coverage_coefficient_of_variation = if !arrays.primary_reads.is_empty() && coverage_mean > 0.0 {
                        let mean_cov = coverage_mean as f64;
                        let n = arrays.primary_reads.len() as f64;
                        let variance: f64 = arrays.primary_reads
                            .iter()
                            .map(|&x| {
                                let diff = x as f64 - mean_cov;
                                diff * diff
                            })
                            .sum::<f64>() / n;
                        let cv = variance.sqrt() / mean_cov;
                        (cv * 1_000_000.0) as f32
                    } else {
                        0.0
                    };

                    let presence = PresenceData {
                        contig_id,
                        contig_name: ref_name.clone(),
                        coverage_pct: coverage_pct as f32,
                        above_expected_aligned_fraction: above_expected,
                        read_count: primary_count,
                        coverage_relative_coverage_roughness,
                        coverage_coefficient_of_variation,
                        coverage_mean,
                        coverage_median,
                        coverage_trimmed_mean,
                    };

                    if let (Some(a), Some(t)) = (accum_ref, ts) {
                        a.feature_calc_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    }

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

    let total_reads = mapped_reads + unmapped_reads;

    if let Some(s) = st.as_mut() {
        s.total_reads = total_reads;
        if let Some(a) = accum.as_ref() {
            s.bam_open_ns = a.bam_open_ns.load(Ordering::Relaxed);
            s.cds_index_ns = a.cds_index_ns.load(Ordering::Relaxed);
            s.streaming_ns = a.streaming_ns.load(Ordering::Relaxed);
            s.feature_calc_ns = a.feature_calc_ns.load(Ordering::Relaxed);
            s.stream_fetch_seek_ns = a.stream_fetch_seek_ns.load(Ordering::Relaxed);
            s.stream_records_ns = a.stream_records_ns.load(Ordering::Relaxed);
            s.stream_process_read_ns = a.stream_process_read_ns.load(Ordering::Relaxed);
            s.stream_finalize_ns = a.stream_finalize_ns.load(Ordering::Relaxed);
            s.mag_encoding_ns = a.mag_encoding_ns.load(Ordering::Relaxed);
        }
    }

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
        is_circular,
        timings: st,
        mag_results,
    })
}

/// Extract contig information from BAM file headers.
/// Deduplicates contigs across all BAM files.
fn extract_contigs_from_bams(bam_files: &[PathBuf], _circularity_map: &HashMap<PathBuf, bool>) -> Result<Vec<ContigInfo>> {
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
    _create_indexes: bool, // Ignored - DuckDB uses zone maps instead of indexes
    extend_db: &Path,
) -> Result<ProcessResult> {
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    // 1. Strict upfront validation — fail early with actionable errors
    if !bam_files.is_empty() {
        validate_inputs(bam_files, modules, genbank_path, assembly_path, extend_db)?;
    }

    // 2. Parse annotations (GenBank/GFF3) or extract contigs from BAMs
    eprintln!("\n### Parsing input files...");

    let has_genbank = !genbank_path.as_os_str().is_empty();
    let is_extending = !extend_db.as_os_str().is_empty();

    // For BAM-only mode, do a quiet @CO-only scan for contig extraction (no log messages).
    // The full detection with logging happens once at step 5 below.
    let preliminary_map = if !bam_files.is_empty() && !has_genbank {
        quick_co_header_scan(bam_files)?
    } else {
        HashMap::new()
    };

    // In MAG mode, per-MAG membership (mag_name → member contig names) is collected during parsing.
    let mut mag_contig_map: Vec<(String, Vec<String>)> = Vec::new();

    let (mut contigs, mut annotations, contig_qualifiers) = if config.view_mode == ViewMode::Mag {
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
                for (local_idx, c) in mag_contigs.iter().enumerate() {
                    member_names.push(c.name.clone());
                    // Annotations/qualifiers reference parser-local ids 1..=N; shift by base_global.
                    let _ = local_idx;
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
                    // member_names already populated; just merge sequences.
                    let fasta_records = crate::parser::parse_fasta(asm)?;
                    for (name, seq) in fasta_records {
                        if let Some(existing) = all_contigs.iter_mut().find(|c| c.name == name) {
                            if existing.sequence.is_none() {
                                existing.sequence = Some(seq);
                            }
                        }
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
        let contigs = extract_contigs_from_bams(bam_files, &preliminary_map)?;
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
        annotations.len()
    );

    eprintln!("Found {} BAM files", bam_files.len());

    // 3. If assembly_path provided, parse FASTA and merge sequences into contigs.
    //    Skipped in MAG mode — sequences are already merged per-MAG above.
    if config.view_mode != ViewMode::Mag && !assembly_path.as_os_str().is_empty() {
        if assembly_path.is_file() {
            eprintln!("\n### Merging sequences from assembly FASTA...");
            merge_sequences_from_fasta(&mut contigs, assembly_path)?;
        } else if assembly_path.is_dir() {
            eprintln!("\n### Merging sequences from assembly FASTA directory...");
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
    if !annotations.is_empty() {
        compute_annotation_sequences(&contigs, &mut annotations);
    }

    if let Some(parent) = output_db.parent() {
        fs::create_dir_all(parent).context("Failed to create output directory")?;
    }

    std::thread::sleep(std::time::Duration::from_millis(50));

    // 4. Create or open database
    // For extend mode: read existing contigs, determine new contigs, open existing DB
    let (db_writer, new_contigs_for_blast, repeats_from_db) = if is_extending {
        // Read existing contig names from the database
        let existing_conn = duckdb::Connection::open(extend_db)
            .with_context(|| format!("Failed to open existing database: {}", extend_db.display()))?;

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
                "SELECT c.Contig_name, c.Contig_length, cs.Sequence \
                 FROM Contig c LEFT JOIN Contig_sequence cs ON c.Contig_id = cs.Contig_id"
            )?;
            let rows = stmt.query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?, row.get::<_, Option<String>>(2)?))
            })?;
            for row in rows {
                let (name, length, seq) = row?;
                existing_contig_names.insert(name.clone());
                existing_contigs.push(ContigInfo {
                    name,
                    aliases: vec![],
                    length: length as usize,
                    sequence: seq.map(|s| s.into_bytes()),
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

        // Load repeats from DB for DTR region detection on existing contigs
        let mut repeats_from_db: Vec<crate::db::RepeatsData> = Vec::new();
        {
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
                            pident: pident_int as f64 / 100.0,
                            is_direct: *is_direct,
                        });
                    }
                }
            }
        }

        // Load CDS annotations from DB for codon change computation (when no -g provided)
        let db_annotations: Vec<FeatureAnnotation> = if !has_genbank && db_has_annotations {
            let mut anns = Vec::new();
            let mut stmt = existing_conn.prepare(
                "SELECT ca.Contig_id, ca.\"Start\", ca.\"End\", ca.Strand, ca.\"Type\", \
                        aseq.Nucleotide_sequence \
                 FROM Contig_annotation_core ca \
                 LEFT JOIN Annotation_sequence aseq ON ca.Annotation_id = aseq.Annotation_id \
                 WHERE ca.\"Type\" = 'CDS'"
            )?;
            let rows = stmt.query_map([], |row| {
                Ok((
                    row.get::<_, i64>(0)?,
                    row.get::<_, i64>(1)?,
                    row.get::<_, i64>(2)?,
                    row.get::<_, i64>(3)?,
                    row.get::<_, String>(4)?,
                    row.get::<_, Option<String>>(5)?,
                ))
            })?;
            for row in rows {
                let (contig_id, start, end, strand, feature_type, nuc_seq) = row?;
                anns.push(FeatureAnnotation {
                    contig_id,
                    start,
                    end,
                    strand,
                    feature_type,
                    nucleotide_sequence: nuc_seq,
                    protein_sequence: None,
                    qualifiers: HashMap::new(),
                    segments: Vec::new(),
                    parent_key: None,
                    self_key: None,
                });
            }
            eprintln!("Loaded {} CDS annotations from existing database", anns.len());
            anns
        } else {
            Vec::new()
        };

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
            new_contig_qualifiers = contig_qualifiers.into_iter()
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
        annotations = if new_annotations.is_empty() && !db_annotations.is_empty() {
            db_annotations
        } else {
            new_annotations.clone()
        };

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
        writer.update_metadata_modification()?;
        writer.write_compression_metadata(
            config.min_aligned_fraction,
            config.min_coverage_depth,
            config.bar_ratio,
            config.variation_percentage,
            config.min_occurrences,
        )?;
        (writer, new_contigs, repeats_from_db)
    } else {
        let new_contigs_for_blast = contigs.clone();
        let writer = DbWriter::create(
            output_db,
            &contigs,
            &annotations,
            &contig_qualifiers,
            !bam_files.is_empty(),
            config.view_mode == ViewMode::Mag,
        )?;
        writer.write_metadata(
            modules,
            config.min_aligned_fraction,
            config.min_coverage_depth,
            config.bar_ratio,
            config.variation_percentage,
            config.min_occurrences,
            config.view_mode.as_str(),
            config.blast,
        )?;
        (writer, new_contigs_for_blast, Vec::new())
    };

    // 5. Auto-detect circularity per sample (now with contigs available for length comparison)
    let circularity_map = if !bam_files.is_empty() {
        detect_all_sample_circularities(bam_files, &contigs)?
    } else {
        HashMap::new()
    };

    // Open timing log early so pre-sample phases (autoblast, GC) are captured
    let mut timing_file: Option<fs::File> = if config.enable_timing {
        match open_timing_log(output_db, config.threads, bam_files.len()) {
            Ok(f) => Some(f),
            Err(e) => { eprintln!("Warning: could not open timing log: {}", e); None }
        }
    } else {
        None
    };

    // 6. If --blast enabled and sequences available, run repeat detection.
    // Contig mode: per-contig BLAST self-alignment (autoblast).
    // MAG mode: per-MAG BLAST handles both intra-contig repeats and inter-contig hits.
    let blast_contigs = &new_contigs_for_blast;
    let has_sequences = blast_contigs.iter().any(|c| c.sequence.is_some());
    let n_contigs_with_seq = blast_contigs.iter().filter(|c| c.sequence.is_some()).count();

    let t_autoblast = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let repeats = if config.blast && has_sequences {
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
            let (reps, inter_count) = run_mag_interblast(blast_contigs, &filtered_manifest, &db_writer)?;
            db_writer.compute_duplication_percentages(&reps)?;
            if inter_count > 0 {
                let t_hit_features = std::time::Instant::now();
                db_writer.write_mag_hit_features()?;
                eprintln!("  Per-position hit features: {:.1}s", t_hit_features.elapsed().as_secs_f64());
            }
            reps
        } else {
            let reps = run_autoblast(blast_contigs, config.threads)?;
            db_writer.write_repeats(&reps)?;
            db_writer.compute_duplication_percentages(&reps)?;
            reps
        }
    } else {
        if !config.blast && has_sequences {
            eprintln!("\n### BLAST repeat detection: skipped (use --blast to enable)");
        }
        Vec::new()
    };
    let autoblast_secs = t_autoblast.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);

    let repeats = if repeats.is_empty() && !repeats_from_db.is_empty() {
        repeats_from_db
    } else {
        repeats
    };

    // 7. Compute and write GC content and GC skew from sequence data (only for new contigs)
    let t_gc = if config.enable_timing { Some(std::time::Instant::now()) } else { None };
    let gc_data: Vec<GCContentData> = blast_contigs
        .iter()
        .filter_map(|contig| {
            contig.sequence.as_ref().map(|seq| {
                let (gc_values, stats) = compute_gc_content(seq, config.gc_params.gc_content_window_size);
                let (skew_values, skew_stats) = compute_gc_skew(seq, config.gc_params.gc_skew_window_size);
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

    if !gc_data.is_empty() {
        eprintln!("\n### Computing GC content for {} contigs...", gc_data.len());
        db_writer.write_gc_content(&gc_data)?;
        db_writer.write_gc_skew(&gc_data)?;
        db_writer.update_contig_gc_stats(&gc_data)?;
    }
    let gc_secs = t_gc.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);

    // Convert repeat detections into per-contig blobs (only for new contigs).
    let new_contig_ids: Vec<i64> = new_contigs_for_blast.iter()
        .filter_map(|c| db_writer.get_contig_id(&c.name))
        .collect();
    db_writer.convert_repeat_blobs(&new_contig_ids)?;

    // Write MAG rows after per-contig GC/duplication stats are finalized.
    if config.view_mode == ViewMode::Mag && !mag_contig_map.is_empty() {
        eprintln!("\n### Writing {} MAGs ({} associations)...",
            mag_contig_map.len(),
            mag_contig_map.iter().map(|(_, v)| v.len()).sum::<usize>(),
        );
        db_writer.write_mags(&mag_contig_map)?;
    }

    // Aggregate per-contig blobs into MAG-scale blobs once. Sample-independent.
    let mag_contig_blob_secs = if config.view_mode == ViewMode::Mag && !mag_contig_map.is_empty() {
        let t = std::time::Instant::now();
        crate::mag_blob::build_mag_contig_blobs_all(&db_writer)?;
        let secs = t.elapsed().as_secs_f64();
        if config.enable_timing {
            eprintln!("MAG contig-blob aggregation: {:.3} s", secs);
        }
        secs
    } else { 0.0 };

    // Write pre-sample timing block
    if let Some(ref mut f) = timing_file {
        use std::io::Write;
        let _ = writeln!(f, "=== Pre-sample phases ===");
        let _ = writeln!(f, "  Contigs with sequence: {}", n_contigs_with_seq);
        let _ = writeln!(f, "  Autoblast            : {:>10.3} s", autoblast_secs);
        let _ = writeln!(f, "  GC content + skew    : {:>10.3} s", gc_secs);
        let _ = writeln!(f, "  MAG contig blobs     : {:>10.3} s", mag_contig_blob_secs);
        let _ = writeln!(f);
    }

    // If no BAM files provided, skip sample processing (genbank-only mode)
    if bam_files.is_empty() {
        eprintln!("\n### No BAM files provided - skipping sample processing");
        eprintln!("Database populated with {} contigs and {} annotations", contigs.len(), annotations.len());

        db_writer.finalize()?;

        return Ok(ProcessResult {
            samples_processed: 0,
            samples_failed: 0,
            total_time_secs: 0.0,
            processing_time_secs: 0.0,
            writing_time_secs: 0.0,
        });
    }

    // 8. Process samples
    eprintln!("\n### Processing {} samples with {} threads", bam_files.len(), config.threads);
    eprintln!("Modules: {}\n", modules.join(", "));

    let flags = ModuleFlags::from_modules(modules);
    if flags.phagetermini && !config.blast {
        eprintln!("WARNING: Phage termini module is enabled without --blast. No terminal repeat (DTR/ITR) detection will be performed.\n");
    }
    let result = process_samples_parallel(&bam_files, &contigs, flags, config, &circularity_map, db_writer, &repeats, &annotations, output_db, timing_file, autoblast_secs, gc_secs, &mag_contig_map)?;

    print_summary(&result, output_db);

    Ok(result)
}

/// For a single sample, identify MAGs that fail the aligned-fraction or
/// coverage-depth thresholds and return the set of member contig names to drop.
/// All member contigs contribute their length to the denominator (even if they
/// have zero reads), so absent contigs correctly penalise the MAG aggregate.
fn filter_failing_mags(
    presences: &[PresenceData],
    mag_contig_map: &[(String, Vec<String>)],
    contig_lengths: &HashMap<String, usize>,
    min_aligned_fraction: f64,
    min_coverage_depth: f64,
) -> (HashSet<String>, usize) {
    let presence_map: HashMap<&str, &PresenceData> = presences
        .iter()
        .map(|p| (p.contig_name.as_str(), p))
        .collect();

    let mut drop_set: HashSet<String> = HashSet::new();
    let mut n_dropped: usize = 0;

    for (_mag_name, members) in mag_contig_map {
        let mut num_af: f64 = 0.0;
        let mut num_cov: f64 = 0.0;
        let mut total_len: f64 = 0.0;

        for cname in members {
            let clen = *contig_lengths.get(cname).unwrap_or(&0) as f64;
            total_len += clen;
            if let Some(p) = presence_map.get(cname.as_str()) {
                num_af += p.coverage_pct as f64 * clen;
                num_cov += p.coverage_mean as f64 * clen;
            }
        }

        if total_len == 0.0 {
            continue;
        }

        let af_mag = num_af / total_len;
        let cov_mag = num_cov / total_len;

        let af_fail = min_aligned_fraction > 0.0 && af_mag < min_aligned_fraction;
        let cov_fail = min_coverage_depth > 0.0 && cov_mag < min_coverage_depth;

        if af_fail || cov_fail {
            n_dropped += 1;
            for cname in members {
                drop_set.insert(cname.clone());
            }
        }
    }

    (drop_set, n_dropped)
}

type MagResultVec = Vec<(String, i64, Vec<(String, crate::blob::EncodedBlob)>, crate::mag_blob::MagCoverageStats)>;

/// Apply a contig drop set to all per-sample data vectors, removing entries
/// whose contig_name is in `drop`. Also filters `mag_results` by resolving
/// which MAGs had their contigs dropped via `mag_contig_map`.
fn apply_mag_filter(
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
    let drop_ids: HashSet<i64> = presences.iter()
        .filter(|p| drop.contains(&p.contig_name))
        .map(|p| p.contig_id)
        .collect();
    let dropped_mag_names: HashSet<&str> = mag_contig_map.iter()
        .filter(|(_, members)| members.iter().any(|c| drop.contains(c)))
        .map(|(name, _)| name.as_str())
        .collect();
    (
        presences.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        packaging.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        misassembly.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        microdiversity.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        side_misassembly.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        topology.into_iter().filter(|p| !drop.contains(&p.contig_name)).collect(),
        feature_blobs.into_iter().filter(|(_, cid, _)| !drop_ids.contains(cid)).collect(),
        mag_results.into_iter().filter(|(name, _, _, _)| !dropped_mag_names.contains(name.as_str())).collect(),
    )
}

fn write_mag_results(
    db_writer: &DbWriter,
    sample_name: &str,
    sample_id: i64,
    mag_results: &MagResultVec,
) -> u64 {
    if mag_results.is_empty() { return 0; }
    let mag_t = std::time::Instant::now();
    if let Ok(conn) = db_writer.lock_conn() {
        for (mag_name, mag_id, blobs, stats) in mag_results {
            if let Err(e) = db_writer.write_mag_blobs(&conn, *mag_id, sample_id, blobs) {
                eprintln!("\n### ERROR writing MAG blobs for {} MAG {}: {:?}", sample_name, mag_name, e);
            }
            if let Err(e) = db_writer.write_mag_coverage(&conn, *mag_id, sample_id, stats) {
                eprintln!("\n### ERROR writing MAG coverage for {} MAG {}: {:?}", sample_name, mag_name, e);
            }
        }
    }
    mag_t.elapsed().as_nanos() as u64
}

/// Holds processed sample data ready for database writing.
struct SampleResult {
    sample_name: String,
    sequencing_type: SequencingType,
    total_reads: u64,
    mapped_reads: u64,
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
    /// MAG-level blobs ready to write. Empty in contig mode.
    mag_results: Vec<(String, i64, Vec<(String, crate::blob::EncodedBlob)>, crate::mag_blob::MagCoverageStats)>,
}

fn process_samples_parallel(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    flags: ModuleFlags,
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    repeats: &[RepeatsData],
    annotations: &[FeatureAnnotation],
    output_db: &Path,
    timing_file: Option<fs::File>,
    autoblast_secs: f64,
    gc_secs: f64,
    mag_contig_map: &[(String, Vec<String>)],
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
        return process_samples_sequential(bam_files, contigs, flags, config, circularity_map, db_writer, is_tty, repeats, annotations, output_db, timing_file, autoblast_secs, gc_secs, mag_contig_map);
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
    let is_tty_writer = is_tty;
    let total_writer = total;

    let timing_threads = config.threads;

    // MAG threshold filtering data (moved into writer thread)
    let is_mag_mode = config.view_mode == ViewMode::Mag;
    let min_af = config.min_aligned_fraction;
    let min_cov = config.min_coverage_depth;
    let mag_map_owned: Vec<(String, Vec<String>)> = mag_contig_map.to_vec();
    let contig_length_map: HashMap<String, usize> = contigs.iter()
        .map(|c| (c.name.clone(), c.length))
        .collect();

    // Build MagMemberLookup before writer thread takes db_writer
    let mag_member_lookup: Option<MagMemberLookup> = if is_mag_mode {
        Some(crate::mag_blob::build_mag_member_lookup(&db_writer)?)
    } else {
        None
    };

    // Spawn dedicated writer thread
    let writer_handle = thread::spawn(move || -> Result<(usize, std::time::Duration, Vec<SampleTimings>, usize)> {
        let write_start = std::time::Instant::now();
        let mut written_count = 0usize;
        let mut mag_drops: usize = 0;
        let mut all_timings: Vec<SampleTimings> = Vec::new();
        let mut timing_file = timing_file;

        // Receive and write samples as they arrive
        for result in rx {
            written_count += 1;

            // MAG-level threshold pre-filter: drop contigs from failing MAGs
            let (presences, packaging, misassembly, microdiversity, side_misassembly, topology, feature_blobs, mag_results) =
                if is_mag_mode && !mag_map_owned.is_empty() {
                    let (drop_set, n_dropped) = filter_failing_mags(
                        &result.presences, &mag_map_owned, &contig_length_map, min_af, min_cov,
                    );
                    mag_drops += n_dropped;
                    if !drop_set.is_empty() {
                        apply_mag_filter(&drop_set, &mag_map_owned, result.presences, result.packaging, result.misassembly, result.microdiversity, result.side_misassembly, result.topology, result.feature_blobs, result.mag_results)
                    } else {
                        (result.presences, result.packaging, result.misassembly, result.microdiversity, result.side_misassembly, result.topology, result.feature_blobs, result.mag_results)
                    }
                } else {
                    (result.presences, result.packaging, result.misassembly, result.microdiversity, result.side_misassembly, result.topology, result.feature_blobs, result.mag_results)
                };

            if presences.is_empty() {
                if result.mapped_reads > 0 {
                    eprintln!(
                        "WARNING: Sample '{}': {} mapped reads but all contigs filtered out \
                         (min_aligned_fraction={:.1}%, min_coverage_depth={:.1}). \
                         Lowering thresholds may recover this sample.",
                        result.sample_name, result.mapped_reads,
                        min_af, min_cov,
                    );
                }
                write_pb_clone.inc(1);
                continue;
            }

            // Insert sample
            let sample_write_start = std::time::Instant::now();
            let t_insert = std::time::Instant::now();
            let sample_id = match db_writer.insert_sample(&result.sample_name, result.sequencing_type.as_str(), result.total_reads, result.mapped_reads, result.is_circular) {
                Ok(id) => id,
                Err(e) => {
                    eprintln!("\nError inserting sample {}: {}", result.sample_name, e);
                    let msg = format!("ERR: {}", result.sample_name);
                    write_pb_clone.set_message(msg.clone());
                    write_pb_clone.inc(1);
                    if !is_tty_writer {
                        eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                    }
                    continue;
                }
            };
            let insert_ns = t_insert.elapsed().as_nanos() as u64;

            // Write all data for this sample
            let mut contig_data_ns: u64 = 0;
            let mut mag_write_ns: u64 = 0;
            let t_contig = std::time::Instant::now();
            if let Err(e) = db_writer.write_sample_data(
                &result.sample_name,
                &presences,
                &packaging,
                &misassembly,
                &microdiversity,
                &side_misassembly,
                &topology,
                &feature_blobs,
                result.is_circular,
            ) {
                eprintln!("\n### ERROR writing data for {}: {:?}", result.sample_name, e);
                let msg = format!("ERR: {}", result.sample_name);
                write_pb_clone.set_message(msg.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                }
            } else {
                contig_data_ns = t_contig.elapsed().as_nanos() as u64;
                mag_write_ns = write_mag_results(&db_writer, &result.sample_name, sample_id, &mag_results);
                write_pb_clone.set_message(result.sample_name.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, result.sample_name);
                }
            }
            write_pb_clone.inc(1);

            // Write per-sample timings incrementally to the log file
            if let Some(mut t) = result.timings {
                t.db_write = sample_write_start.elapsed();
                t.insert_sample_ns = insert_ns;
                t.write_contig_data_ns = contig_data_ns;
                t.mag_write_ns = mag_write_ns;
                if let Some(ref mut f) = timing_file {
                    let _ = write_sample_timing(f, &t, timing_threads);
                }
                all_timings.push(t);
            }
        }

        // Finalize database
        db_writer.finalize()?;
        write_pb_clone.finish();

        Ok((written_count, write_start.elapsed(), all_timings, mag_drops))
    });

    // Build shared lookups once — avoids rebuilding per sample (critical for 2M+ contigs).
    let contig_by_name: HashMap<&str, usize> = {
        let mut map = HashMap::new();
        for (idx, c) in contigs.iter().enumerate() {
            map.insert(c.name.as_str(), idx);
            for alias in &c.aliases {
                map.entry(alias.as_str()).or_insert(idx);
            }
        }
        map
    };
    let annos_by_contig: HashMap<i64, Vec<&FeatureAnnotation>> = if flags.mapping_metrics
        && !annotations.is_empty()
    {
        let mut map: HashMap<i64, Vec<&FeatureAnnotation>> = HashMap::new();
        for a in annotations.iter() {
            map.entry(a.contig_id).or_default().push(a);
        }
        map
    } else {
        HashMap::new()
    };

    // Process samples in parallel, sending to channel immediately
    // send() blocks if channel is full (backpressure)
    let mag_member_lookup_ref = mag_member_lookup.as_ref();
    bam_files.par_iter().for_each(|bam_path| {
        let sample_start = std::time::Instant::now();

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        match process_sample(bam_path, contigs, flags, config, sample_circular, repeats, annotations, &contig_by_name, &annos_by_contig, mag_member_lookup_ref) {
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
                    eprintln!("Processing: {}", msg);
                }

                let _ = tx.send(SampleResult {
                    sample_name: r.sample_name,
                    sequencing_type: r.seq_type,
                    total_reads: r.total_reads,
                    mapped_reads: r.mapped_reads,
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
                });
            }
            Err(e) => {
                eprintln!("\nError processing {}: {}", bam_path.display(), e);
                let sample_time = sample_start.elapsed().as_secs_f64();
                completed_count.fetch_add(1, Ordering::SeqCst);
                failed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("ERR ({:.2}s)", sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("Processing: {}", msg);
                }
            }
        }
    });

    // Drop sender to signal writer thread that processing is complete
    drop(tx);

    process_pb.finish();
    let processing_time = start_time.elapsed();

    // Wait for writer thread to finish
    let (written_count, writing_time, all_timings, mag_drops) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    if mag_drops > 0 {
        eprintln!("### Dropped {} (MAG, Sample) groups that failed thresholds", mag_drops);
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
                if let Err(e) = finish_timing_log(&mut f, &all_timings, elapsed, autoblast_secs, gc_secs) {
                    eprintln!("Warning: failed to write timing totals: {}", e);
                }
            }
            Err(e) => eprintln!("Warning: could not reopen timing log: {}", e),
        }
    }

    Ok(ProcessResult {
        samples_processed: written_count,
        samples_failed: failed,
        total_time_secs: elapsed.as_secs_f64(),
        processing_time_secs: processing_time.as_secs_f64(),
        writing_time_secs: writing_time.as_secs_f64(),
    })
}

/// Sequential processing for single-threaded mode.
/// Process one sample → write it → repeat. No channel overhead.
fn process_samples_sequential(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    flags: ModuleFlags,
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    is_tty: bool,
    repeats: &[RepeatsData],
    annotations: &[FeatureAnnotation],
    _output_db: &Path,
    timing_file: Option<fs::File>,
    autoblast_secs: f64,
    gc_secs: f64,
    mag_contig_map: &[(String, Vec<String>)],
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let start_time = std::time::Instant::now();
    let is_mag_mode = config.view_mode == ViewMode::Mag;
    let contig_length_map: HashMap<String, usize> = contigs.iter()
        .map(|c| (c.name.clone(), c.length))
        .collect();
    let mut mag_drops: usize = 0;

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
    let contig_by_name: HashMap<&str, usize> = {
        let mut map = HashMap::new();
        for (idx, c) in contigs.iter().enumerate() {
            map.insert(c.name.as_str(), idx);
            for alias in &c.aliases {
                map.entry(alias.as_str()).or_insert(idx);
            }
        }
        map
    };
    let annos_by_contig: HashMap<i64, Vec<&FeatureAnnotation>> = if flags.mapping_metrics
        && !annotations.is_empty()
    {
        let mut map: HashMap<i64, Vec<&FeatureAnnotation>> = HashMap::new();
        for a in annotations.iter() {
            map.entry(a.contig_id).or_default().push(a);
        }
        map
    } else {
        HashMap::new()
    };

    for bam_path in bam_files {
        let sample_start = std::time::Instant::now();

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        match process_sample(bam_path, contigs, flags, config, sample_circular, repeats, annotations, &contig_by_name, &annos_by_contig, mag_member_lookup.as_ref()) {
            Ok(mut r) => {
                let process_elapsed = sample_start.elapsed();
                processing_time_total += process_elapsed;

                let write_start = std::time::Instant::now();

                // MAG-level threshold pre-filter
                let (presences, packaging, misassembly, microdiversity, side_misassembly, topology, feature_blobs, mag_results) =
                    if is_mag_mode && !mag_contig_map.is_empty() {
                        let (drop_set, n_dropped) = filter_failing_mags(
                            &r.presences, mag_contig_map, &contig_length_map,
                            config.min_aligned_fraction, config.min_coverage_depth,
                        );
                        mag_drops += n_dropped;
                        if !drop_set.is_empty() {
                            apply_mag_filter(&drop_set, mag_contig_map, r.presences, r.packaging, r.misassembly, r.microdiversity, r.side_misassembly, r.topology, r.feature_blobs, r.mag_results)
                        } else {
                            (r.presences, r.packaging, r.misassembly, r.microdiversity, r.side_misassembly, r.topology, r.feature_blobs, r.mag_results)
                        }
                    } else {
                        (r.presences, r.packaging, r.misassembly, r.microdiversity, r.side_misassembly, r.topology, r.feature_blobs, r.mag_results)
                    };

                if presences.is_empty() {
                    if r.mapped_reads > 0 {
                        eprintln!(
                            "WARNING: Sample '{}': {} mapped reads but all contigs filtered out \
                             (min_aligned_fraction={:.1}%, min_coverage_depth={:.1}). \
                             Lowering thresholds may recover this sample.",
                            r.sample_name, r.mapped_reads,
                            config.min_aligned_fraction, config.min_coverage_depth,
                        );
                    }
                    pb.inc(1);
                    continue;
                }

                // Insert sample
                let t_insert = std::time::Instant::now();
                let sample_id = match db_writer.insert_sample(&r.sample_name, r.seq_type.as_str(), r.total_reads, r.mapped_reads, r.is_circular) {
                    Ok(id) => id,
                    Err(e) => {
                        eprintln!("\nError inserting sample {}: {}", r.sample_name, e);
                        failed += 1;
                        pb.inc(1);
                        continue;
                    }
                };
                let insert_ns = t_insert.elapsed().as_nanos() as u64;

                // Write sample data
                let mut contig_data_ns: u64 = 0;
                let mut mag_write_ns: u64 = 0;
                let t_contig = std::time::Instant::now();
                if let Err(e) = db_writer.write_sample_data(
                    &r.sample_name,
                    &presences,
                    &packaging,
                    &misassembly,
                    &microdiversity,
                    &side_misassembly,
                    &topology,
                    &feature_blobs,
                    r.is_circular,
                ) {
                    eprintln!("\nError writing data for {}: {}", r.sample_name, e);
                    failed += 1;
                } else {
                    contig_data_ns = t_contig.elapsed().as_nanos() as u64;
                    mag_write_ns = write_mag_results(&db_writer, &r.sample_name, sample_id, &mag_results);
                    processed += 1;
                }

                writing_time_total += write_start.elapsed();

                // Write per-sample timings incrementally to the log file
                if let Some(ref mut t) = r.timings {
                    t.db_write = write_start.elapsed();
                    t.insert_sample_ns = insert_ns;
                    t.write_contig_data_ns = contig_data_ns;
                    t.mag_write_ns = mag_write_ns;
                    t.total = sample_start.elapsed();
                    if let Some(ref mut f) = timing_file {
                        let _ = write_sample_timing(f, t, config.threads);
                    }
                    all_timings.push(t.clone());
                }

                let total_sample_time = sample_start.elapsed().as_secs_f64();
                let msg = format!("{} ({:.2}s)", r.sample_name, total_sample_time);
                pb.set_message(msg.clone());
                pb.inc(1);
                if !is_tty {
                    eprintln!("{}", msg);
                }
            }
            Err(e) => {
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
    }

    if is_mag_mode && mag_drops > 0 {
        eprintln!("### Dropped {} (MAG, Sample) groups that failed thresholds", mag_drops);
    } else if !is_mag_mode {
        let dropped = config.contig_drops_counter.load(Ordering::Relaxed);
        if dropped > 0 {
            eprintln!("### Dropped {} (Contig, Sample) pairs that failed thresholds", dropped);
        }
    }

    // Finalize database
    db_writer.finalize()?;
    pb.finish_with_message("Done");
    if !is_tty {
        eprintln!("Done ({:.2}s)", start_time.elapsed().as_secs_f64());
    }

    // Append grand totals to the timing log
    if config.enable_timing && !all_timings.is_empty() {
        if let Some(ref mut f) = timing_file {
            if let Err(e) = finish_timing_log(f, &all_timings, start_time.elapsed(), autoblast_secs, gc_secs) {
                eprintln!("Warning: failed to write timing totals: {}", e);
            }
        }
    }

    Ok(ProcessResult {
        samples_processed: processed,
        samples_failed: failed,
        total_time_secs: start_time.elapsed().as_secs_f64(),
        processing_time_secs: processing_time_total.as_secs_f64(),
        writing_time_secs: writing_time_total.as_secs_f64(),
    })
}

fn ns_to_s(ns: u64) -> f64 {
    ns as f64 / 1e9
}

fn timing_log_path(output_db: &Path) -> PathBuf {
    let mut p = output_db.as_os_str().to_owned();
    p.push(".timings.log");
    PathBuf::from(p)
}

/// Open the timing log and write the header. Returns the open file handle
/// so callers can append per-sample sections as they finish.
fn open_timing_log(output_db: &Path, threads: usize, n_samples: usize) -> Result<fs::File> {
    use std::io::Write;
    let log_path = timing_log_path(output_db);
    let mut f = fs::File::create(&log_path)
        .with_context(|| format!("Failed to create timing log: {}", log_path.display()))?;

    writeln!(f, "=============================================================")?;
    writeln!(f, "theBIGbam calculate -- timing report")?;
    writeln!(f, "Output DB    : {}", output_db.display())?;
    writeln!(f, "Threads      : {}", threads)?;
    writeln!(f, "Samples      : {}", n_samples)?;
    writeln!(f, "=============================================================\n")?;

    eprintln!("Timing log: {}", log_path.display());
    Ok(f)
}

/// Append one sample's timing breakdown to an already-open timing log.
fn write_sample_timing(f: &mut fs::File, t: &SampleTimings, threads: usize) -> Result<()> {
    use std::io::Write;
    writeln!(f, "Sample: {}", t.sample_name)?;
    writeln!(f, "  BAM path             : {}", t.bam_path)?;
    writeln!(f, "  Contigs in BAM       : {}", t.n_contigs_in_bam)?;
    writeln!(f, "  Contigs processed    : {}", t.n_contigs_processed)?;
    writeln!(f, "  Total reads          : {}", t.total_reads)?;
    writeln!(f, "  ---- Wall-clock phases ----")?;
    writeln!(f, "  Header read          : {:>10.3} s", t.header_read.as_secs_f64())?;
    writeln!(f, "  Total read count     : {:>10.3} s", t.total_read_count.as_secs_f64())?;
    writeln!(f, "  Build lookups        : {:>10.3} s", t.build_lookups.as_secs_f64())?;
    writeln!(f, "  Filter tids (pre-par): {:>10.3} s", t.filter_tids.as_secs_f64())?;
    writeln!(f, "  Par-iter wall        : {:>10.3} s", t.par_iter_wall.as_secs_f64())?;
    writeln!(f, "  Merge contig results : {:>10.3} s", t.merge_results.as_secs_f64())?;
    writeln!(f, "  DB write             : {:>10.3} s", t.db_write.as_secs_f64())?;
    writeln!(f, "  TOTAL per sample     : {:>10.3} s", t.total.as_secs_f64())?;
    writeln!(f, "  ---- Par-iter sub-phase sums (CPU-seconds across {} threads) ----", threads)?;
    writeln!(f, "  BAM reader open      : {:>10.3} s", ns_to_s(t.bam_open_ns))?;
    writeln!(f, "  CDS index build      : {:>10.3} s", ns_to_s(t.cds_index_ns))?;
    writeln!(f, "  process_contig_stream: {:>10.3} s", ns_to_s(t.streaming_ns))?;
    writeln!(f, "    . fetch seek       : {:>10.3} s", ns_to_s(t.stream_fetch_seek_ns))?;
    writeln!(f, "    . record loop (htslib+prep): {:>10.3} s", ns_to_s(t.stream_records_ns))?;
    writeln!(f, "    . process_read     : {:>10.3} s", ns_to_s(t.stream_process_read_ns))?;
    writeln!(f, "    . finalize+cov     : {:>10.3} s", ns_to_s(t.stream_finalize_ns))?;
    writeln!(f, "  Feature calculation  : {:>10.3} s", ns_to_s(t.feature_calc_ns))?;
    writeln!(f, "  MAG blob encoding    : {:>10.3} s", ns_to_s(t.mag_encoding_ns))?;
    let sub_sum_ns = t.bam_open_ns + t.cds_index_ns
                   + t.streaming_ns + t.feature_calc_ns + t.mag_encoding_ns;
    writeln!(f, "  (sum of above)       : {:>10.3} s  ({}x threads -> wall ~ {:.1} s)",
             ns_to_s(sub_sum_ns), threads, ns_to_s(sub_sum_ns) / threads as f64)?;
    writeln!(f, "  ---- DB write sub-phases ----")?;
    writeln!(f, "  Insert sample        : {:>10.3} s", ns_to_s(t.insert_sample_ns))?;
    writeln!(f, "  Write contig data    : {:>10.3} s", ns_to_s(t.write_contig_data_ns))?;
    writeln!(f, "  Write MAG data       : {:>10.3} s", ns_to_s(t.mag_write_ns))?;
    writeln!(f)?;
    Ok(())
}

/// Append grand totals to the timing log and close it.
fn finish_timing_log(
    f: &mut fs::File,
    timings: &[SampleTimings],
    total_wall: std::time::Duration,
    autoblast_secs: f64,
    gc_secs: f64,
) -> Result<()> {
    use std::io::Write;
    writeln!(f, "=============================================================")?;
    writeln!(f, "GRAND TOTALS")?;
    writeln!(f, "=============================================================")?;
    writeln!(f, "  ---- Pre-sample phases ----")?;
    writeln!(f, "  Autoblast            : {:>10.3} s", autoblast_secs)?;
    writeln!(f, "  GC content + skew    : {:>10.3} s", gc_secs)?;
    writeln!(f, "  ---- Per-sample phases ----")?;
    let sum_dur = |g: fn(&SampleTimings) -> std::time::Duration| -> f64 {
        timings.iter().map(|t| g(t).as_secs_f64()).sum()
    };
    let sum_ns = |g: fn(&SampleTimings) -> u64| -> f64 {
        timings.iter().map(|t| ns_to_s(g(t))).sum()
    };
    writeln!(f, "  Header read          : {:>10.3} s", sum_dur(|t| t.header_read))?;
    writeln!(f, "  Total read count     : {:>10.3} s", sum_dur(|t| t.total_read_count))?;
    writeln!(f, "  Build lookups        : {:>10.3} s", sum_dur(|t| t.build_lookups))?;
    writeln!(f, "  Filter tids (pre-par): {:>10.3} s", sum_dur(|t| t.filter_tids))?;
    writeln!(f, "  Par-iter wall        : {:>10.3} s", sum_dur(|t| t.par_iter_wall))?;
    writeln!(f, "  Merge contig results : {:>10.3} s", sum_dur(|t| t.merge_results))?;
    writeln!(f, "  DB write             : {:>10.3} s", sum_dur(|t| t.db_write))?;
    writeln!(f, "  Total per-sample sum : {:>10.3} s", sum_dur(|t| t.total))?;
    writeln!(f, "  Wall-clock (calc)    : {:>10.3} s", total_wall.as_secs_f64())?;
    writeln!(f, "  ---- Par-iter sub-phase sums ----")?;
    writeln!(f, "  BAM reader open      : {:>10.3} s", sum_ns(|t| t.bam_open_ns))?;
    writeln!(f, "  CDS index build      : {:>10.3} s", sum_ns(|t| t.cds_index_ns))?;
    writeln!(f, "  process_contig_stream: {:>10.3} s", sum_ns(|t| t.streaming_ns))?;
    writeln!(f, "    . fetch seek       : {:>10.3} s", sum_ns(|t| t.stream_fetch_seek_ns))?;
    writeln!(f, "    . record loop (htslib+prep): {:>10.3} s", sum_ns(|t| t.stream_records_ns))?;
    writeln!(f, "    . process_read     : {:>10.3} s", sum_ns(|t| t.stream_process_read_ns))?;
    writeln!(f, "    . finalize+cov     : {:>10.3} s", sum_ns(|t| t.stream_finalize_ns))?;
    writeln!(f, "  Feature calculation  : {:>10.3} s", sum_ns(|t| t.feature_calc_ns))?;
    writeln!(f, "  MAG blob encoding    : {:>10.3} s", sum_ns(|t| t.mag_encoding_ns))?;
    writeln!(f, "  ---- DB write sub-phase sums ----")?;
    writeln!(f, "  Insert sample        : {:>10.3} s", sum_ns(|t| t.insert_sample_ns))?;
    writeln!(f, "  Write contig data    : {:>10.3} s", sum_ns(|t| t.write_contig_data_ns))?;
    writeln!(f, "  Write MAG data       : {:>10.3} s", sum_ns(|t| t.mag_write_ns))?;
    Ok(())
}

fn print_summary(result: &ProcessResult, output_db: &Path) {
    eprintln!();
    eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    eprintln!("### Complete");
    eprintln!();
    eprintln!("  Samples processed: {}/{}", result.samples_processed, result.samples_processed + result.samples_failed);
    eprintln!("  Total time:        {:.2}s", result.total_time_secs);
    eprintln!();
    eprintln!("  Output: {:?}", output_db);
}
