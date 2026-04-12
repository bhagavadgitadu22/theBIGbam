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
use std::sync::mpsc::sync_channel;
use std::thread;

use crate::bam_reader::{detect_sequencing_type, get_total_read_count, process_contig_streaming, BamPhaseTimings};
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

/// Configuration for processing.
#[derive(Clone)]
pub struct ProcessConfig {
    pub threads: usize,
    pub min_aligned_fraction: f64,
    pub min_coverage_depth: f64,
    /// Relative tolerance for RLE compression (e.g., 0.1 = 10% change threshold)
    pub curve_ratio: f64,
    pub bar_ratio: f64,
    /// Sequencing type: if None, auto-detect per sample; if Some, use for all samples
    pub sequencing_type: Option<SequencingType>,
    /// Phage termini detection configuration
    pub phagetermini_config: PhageTerminiConfig,
    /// GC content and GC skew parameters
    pub gc_params: GCParams,
    /// Minimum absolute event count for a sparse feature position to be kept.
    /// Applied in addition to bar_ratio filtering: position kept only if
    /// `value > coverage × bar_ratio AND value > min_occurrences`.
    /// Default: 2.
    pub min_occurrences: u32,
    /// When true, record per-sample phase timings and write a timing log
    /// alongside the output database. Hidden CLI flag (--time).
    pub enable_timing: bool,
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
    // Filled by the sample-level caller, not process_sample itself.
    pub db_write: std::time::Duration,
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
    genbank_path: &Path,
    assembly_path: &Path,
    extend_db: &Path,
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

    // 3. Sequence check (only when modules include Phage termini, skip in extend mode
    //    because sequences are already stored in the DB's Contig_sequence table)
    let is_extending = !extend_db.as_os_str().is_empty();
    if flags.phagetermini && !is_extending {
        let has_genbank = !genbank_path.as_os_str().is_empty() && genbank_path.exists();
        let has_assembly = !assembly_path.as_os_str().is_empty() && assembly_path.exists();
        if !has_genbank && !has_assembly {
            errors.push(
                "ERROR: Phage termini module requires sequence data for terminal repeat detection.\n\
                 Provide either:\n\
                 \x20 - A GenBank file with embedded sequences via -g/--genbank\n\
                 \x20 - An assembly FASTA file via -a/--assembly\n\
                 Or remove 'Phage termini' from -m/--modules."
                    .to_string(),
            );
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

/// Run autoblast (self-BLAST) on contigs that have sequence data.
/// Uses rayon for per-contig parallelization.
/// If blastn is not on PATH, prints a warning and returns empty vec.
fn run_autoblast(contigs: &[ContigInfo], _threads: usize) -> Result<Vec<RepeatsData>> {
    // Filter contigs that have sequence data
    let contigs_with_seq: Vec<&ContigInfo> = contigs.iter().filter(|c| c.sequence.is_some()).collect();
    if contigs_with_seq.is_empty() {
        return Ok(Vec::new());
    }

    // Check if blastn is on PATH
    match std::process::Command::new("blastn").arg("-version").output() {
        Ok(output) if output.status.success() => {}
        _ => {
            eprintln!("Warning: blastn not found on PATH. Skipping repeat detection.");
            eprintln!("Install BLAST+ to enable terminal repeat detection.");
            return Ok(Vec::new());
        }
    }

    eprintln!("\n### Running autoblast on {} contigs with sequences...", contigs_with_seq.len());

    // Use rayon to parallelize across contigs
    let all_repeats: Vec<Vec<RepeatsData>> = contigs_with_seq
        .par_iter()
        .filter_map(|contig| {
            let seq = contig.sequence.as_ref()?;

            // Write contig sequence to temp FASTA file
            let mut temp_file = match tempfile::NamedTempFile::new() {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Warning: Failed to create temp file for contig '{}': {}", contig.name, e);
                    return None;
                }
            };

            use std::io::Write;
            if let Err(e) = writeln!(temp_file, ">{}", contig.name) {
                eprintln!("Warning: Failed to write temp FASTA for '{}': {}", contig.name, e);
                return None;
            }
            // Write sequence in 80-char lines
            for chunk in seq.chunks(80) {
                if let Err(e) = temp_file.write_all(chunk) {
                    eprintln!("Warning: Failed to write temp FASTA for '{}': {}", contig.name, e);
                    return None;
                }
                if let Err(e) = writeln!(temp_file) {
                    eprintln!("Warning: Failed to write temp FASTA for '{}': {}", contig.name, e);
                    return None;
                }
            }

            // Flush to ensure all data is written before blastn reads it
            if let Err(e) = temp_file.flush() {
                eprintln!("Warning: Failed to flush temp FASTA for '{}': {}", contig.name, e);
                return None;
            }

            let temp_path = temp_file.path().to_path_buf();

            // Run blastn -query temp.fasta -subject temp.fasta -evalue 1e-10 -outfmt 6
            let output = match std::process::Command::new("blastn")
                .arg("-query").arg(&temp_path)
                .arg("-subject").arg(&temp_path)
                .arg("-evalue").arg("1e-10")
                .arg("-outfmt").arg("6")
                .output()
            {
                Ok(o) => o,
                Err(e) => {
                    eprintln!("Warning: blastn failed for contig '{}': {}", contig.name, e);
                    return None;
                }
            };

            if !output.status.success() {
                let stderr = String::from_utf8_lossy(&output.stderr);
                eprintln!("Warning: blastn returned non-zero for contig '{}': {}", contig.name, stderr);
                return None;
            }

            // Parse stdout (outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)
            let stdout = String::from_utf8_lossy(&output.stdout);
            let mut repeats = Vec::new();

            for line in stdout.lines() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 10 {
                    continue;
                }

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
                    contig_name: contig.name.clone(),
                    position1: qstart,
                    position2: qend,
                    position1prime: sstart,
                    position2prime: send,
                    pident,
                    is_direct,
                });
            }

            // Temp file auto-cleaned on drop
            Some(repeats)
        })
        .collect();

    let total: Vec<RepeatsData> = all_repeats.into_iter().flatten().collect();
    eprintln!("Found {} repeat hits across {} contigs", total.len(), contigs_with_seq.len());
    Ok(total)
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
    contig_length: usize,
    config: &ProcessConfig,
    is_circular: bool,
    seq_type: SequencingType,
    flags: ModuleFlags,
    primary_count: u64,
    repeats: &[RepeatsData],
    cds_index: Option<&CdsIndex>,
    blob_output: &mut Vec<(String, String, crate::blob::EncodedBlob)>,
) -> (Option<PackagingData>, Option<(MisassemblyData, MicrodiversityData, SideMisassemblyData, TopologyData)>) {
    use crate::blob::{encode_dense_blob, encode_sparse_blob, EventMeta, MetadataFlags, ValueScale,
                       codon_category_to_id, codon_to_id, aa_to_id};
    let pt_config = config.phagetermini_config;

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
    if flags.coverage {
        // MAPQ - average mapping quality per position (needed for blob encoding below)
        let mapq_f64: Vec<f64> = arrays.sum_mapq.iter()
            .zip(&arrays.primary_reads)
            .map(|(&sum, &count)| if count > 0 { sum as f64 / count as f64 } else { 0.0 })
            .collect();

        // === BLOB encoding for dense coverage features ===
        let clen = contig_length as u32;
        let cn = contig_name.to_string();

        // primary_reads: raw i32
        let pr_i32: Vec<i32> = arrays.primary_reads.iter().map(|&x| x as i32).collect();
        blob_output.push(("primary_reads".into(), cn.clone(), encode_dense_blob(&pr_i32, ValueScale::Raw, clen)));

        // plus/minus strand
        let pp_i32: Vec<i32> = arrays.primary_reads_plus_only.iter().map(|&x| x as i32).collect();
        blob_output.push(("primary_reads_plus_only".into(), cn.clone(), encode_dense_blob(&pp_i32, ValueScale::Raw, clen)));
        let pm_i32: Vec<i32> = arrays.primary_reads_minus_only.iter().map(|&x| x as i32).collect();
        blob_output.push(("primary_reads_minus_only".into(), cn.clone(), encode_dense_blob(&pm_i32, ValueScale::Raw, clen)));

        // secondary, supplementary
        let sec_i32: Vec<i32> = arrays.secondary_reads.iter().map(|&x| x as i32).collect();
        blob_output.push(("secondary_reads".into(), cn.clone(), encode_dense_blob(&sec_i32, ValueScale::Raw, clen)));
        let sup_i32: Vec<i32> = arrays.supplementary_reads.iter().map(|&x| x as i32).collect();
        blob_output.push(("supplementary_reads".into(), cn.clone(), encode_dense_blob(&sup_i32, ValueScale::Raw, clen)));

        // MAPQ: stored as ×100 integer
        let mapq_i32: Vec<i32> = mapq_f64.iter().map(|&x| (x * 100.0).round() as i32).collect();
        blob_output.push(("mapq".into(), cn.clone(), encode_dense_blob(&mapq_i32, ValueScale::Times100, clen)));
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
    let cn = contig_name.to_string();
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
            if val > cov * bar_threshold && val > min_occ as f64 {
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
            let flags = MetadataFlags { sparse: true, has_stats: false, has_sequence: true, has_codons: true };
            blob_output.push(("mismatches".into(), cn.clone(),
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, ValueScale::Times1000, clen)));
        }

        // deletions (value only)
        {
            let (pos, vals) = filter_sparse(&deletions_f64, &primary_reads_f64);
            let flags = MetadataFlags { sparse: true, ..Default::default() };
            blob_output.push(("deletions".into(), cn.clone(),
                encode_sparse_blob(&pos, &vals, None, flags, ValueScale::Times1000, clen)));
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
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false };
            blob_output.push(("insertions".into(), cn.clone(),
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, ValueScale::Times1000, clen)));
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
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false };
            blob_output.push(("left_clippings".into(), cn.clone(),
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, ValueScale::Times1000, clen)));
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
            let flags = MetadataFlags { sparse: true, has_stats: true, has_sequence: true, has_codons: false };
            blob_output.push(("right_clippings".into(), cn.clone(),
                encode_sparse_blob(&pos, &vals, Some(&meta), flags, ValueScale::Times1000, clen)));
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
        blob_output.push(("splicings".into(), cn.clone(),
            encode_sparse_blob(&spl_pos, &spl_vals, None, flags, ValueScale::Raw, clen)));
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
        let cn = contig_name.to_string();
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
                if val > cov * bar_threshold && val > min_occ as f64 {
                    positions.push(i as u32);
                    vals.push((val / cov * 1000.0).round() as i32);
                }
            }
            (positions, vals)
        };

        let (pos, vals) = filter_sparse_pr(&non_inward_f64, &primary_reads_f64);
        let sp_flags = MetadataFlags { sparse: true, ..Default::default() };
        blob_output.push(("non_inward_pairs".into(), cn.clone(),
            encode_sparse_blob(&pos, &vals, None, sp_flags, ValueScale::Times1000, clen)));

        let (pos, vals) = filter_sparse_pr(&mate_unmapped_f64, &primary_reads_f64);
        blob_output.push(("mate_not_mapped".into(), cn.clone(),
            encode_sparse_blob(&pos, &vals, None, sp_flags, ValueScale::Times1000, clen)));

        let (pos, vals) = filter_sparse_pr(&mate_other_contig_f64, &primary_reads_f64);
        blob_output.push(("mate_on_another_contig".into(), cn.clone(),
            encode_sparse_blob(&pos, &vals, None, sp_flags, ValueScale::Times1000, clen)));

        // insert_sizes: dense curve
        let is_i32: Vec<i32> = values.iter().map(|&x| (x * 10.0).round() as i32).collect();
        blob_output.push(("insert_sizes".into(), cn.clone(), encode_dense_blob(&is_i32, ValueScale::Times10, clen)));
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
        let cn = contig_name.to_string();
        let rl_i32: Vec<i32> = values.iter().map(|&x| (x * 10.0).round() as i32).collect();
        blob_output.push(("read_lengths".into(), cn, encode_dense_blob(&rl_i32, ValueScale::Times10, clen)));
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
            let cn = contig_name.to_string();
            let bar_threshold = config.bar_ratio * 0.01;
            let min_occ = config.min_occurrences;

            // coverage_reduced: dense curve
            let cr_i32: Vec<i32> = arrays.coverage_reduced.iter().map(|&x| x as i32).collect();
            blob_output.push(("coverage_reduced".into(), cn.clone(), encode_dense_blob(&cr_i32, ValueScale::Raw, clen)));

            // reads_starts: sparse with median + sequence
            let cov_reduced_f64: Vec<f64> = arrays.coverage_reduced.iter().map(|&x| x as f64).collect();
            let mut rs_pos = Vec::new();
            let mut rs_vals = Vec::new();
            for i in 0..reads_starts_original.len().min(cov_reduced_f64.len()) {
                let val = reads_starts_original[i];
                let cov = cov_reduced_f64[i];
                if val > cov * bar_threshold && val > min_occ as f64 {
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
            let pt_flags = MetadataFlags { sparse: true, has_stats: false, has_sequence: true, has_codons: false };
            blob_output.push(("reads_starts".into(), cn.clone(),
                encode_sparse_blob(&rs_pos, &rs_vals, Some(&rs_meta), pt_flags, ValueScale::Times1000, clen)));

            // reads_ends: sparse with median + sequence
            let mut re_pos = Vec::new();
            let mut re_vals = Vec::new();
            for i in 0..reads_ends_original.len().min(cov_reduced_f64.len()) {
                let val = reads_ends_original[i];
                let cov = cov_reduced_f64[i];
                if val > cov * bar_threshold && val > min_occ as f64 {
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
            blob_output.push(("reads_ends".into(), cn.clone(),
                encode_sparse_blob(&re_pos, &re_vals, Some(&re_meta), pt_flags, ValueScale::Times1000, clen)));
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
            );

            // Only return PackagingData if there's a detected mechanism (not "No_packaging")
            // All termini (both kept and discarded) are included with filtering metadata
            if mechanism != "No_packaging" {
                Some(PackagingData {
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

/// Process one BAM file using streaming (single-pass, optimized).
/// Parallelizes at the contig level for samples with many contigs.
/// Returns (feature_blobs, presences, packaging, misassembly, microdiversity, side_misassembly, topology, sample_name, seq_type, total_reads, mapped_reads)
pub fn process_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    is_circular: bool,
    repeats: &[RepeatsData],
    annotations: &[crate::types::FeatureAnnotation],
) -> Result<(Vec<(String, String, crate::blob::EncodedBlob)>, Vec<PresenceData>, Vec<PackagingData>, Vec<MisassemblyData>, Vec<MicrodiversityData>, Vec<SideMisassemblyData>, Vec<TopologyData>, String, SequencingType, u64, u64, bool, Option<SampleTimings>)> {
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "");

    let flags = ModuleFlags::from_modules(modules);

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
    let temp_bam = IndexedReader::from_path(bam_path)
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

    // --- Phase: total_read_count ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Get total read count from BAM index (fast, no iteration)
    let total_reads = get_total_read_count(bam_path)?;

    if let (Some(s), Some(t)) = (st.as_mut(), t0) {
        s.total_read_count = t.elapsed();
        s.total_reads = total_reads;
    }

    // --- Phase: build_lookups ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Build name -> index lookup once; avoids O(n_refs * n_contigs) linear scans
    // in the per-contig par_iter below (critical for projects with ~500k contigs).
    let contig_by_name: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(idx, c)| (c.name.as_str(), idx))
        .collect();

    // Pre-group annotations by contig_id once (only when mapping_metrics is on).
    // Avoids the O(n_annotations) linear filter inside CdsIndex::from_annotations
    // on every contig. For ~500k contigs x ~15M annotations this is a ~500k-fold
    // reduction in comparisons.
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

    if let (Some(s), Some(t)) = (st.as_mut(), t0) {
        s.build_lookups = t.elapsed();
    }

    // --- Phase: filter_tids (single-threaded) ---
    // Walk the BAM header once and materialise only the tids whose reference
    // name is also a known contig. This spares the par_iter from opening the
    // BAM for thousands of references we'd immediately discard.
    let t_filter = if time_on { Some(std::time::Instant::now()) } else { None };
    let jobs: Vec<(u32, i64, String, usize)> = (0..n_refs)
        .filter_map(|tid| {
            let ref_name = std::str::from_utf8(temp_bam.header().tid2name(tid)).ok()?;
            let contig_idx = contig_by_name.get(ref_name)?;
            let ref_length = temp_bam.header().target_len(tid).unwrap_or(0) as usize;
            Some((tid, (*contig_idx + 1) as i64, ref_name.to_string(), ref_length))
        })
        .collect();
    drop(temp_bam);
    if let (Some(s), Some(t)) = (st.as_mut(), t_filter) {
        s.filter_tids = t.elapsed();
    }

    // --- Phase: par_iter ---
    let accum = if time_on { Some(TimingAccum::new()) } else { None };
    let accum_ref = accum.as_ref();
    let t_par = if time_on { Some(std::time::Instant::now()) } else { None };

    // Each rayon worker chunk opens the BAM once via map_init and reuses it
    // across all of its jobs — the reader repositions cheaply on each fetch().
    // For a 362k-reference BAM this collapses ~362k opens to ~tens.
    let results: Vec<_> = jobs
        .par_iter()
        .map_init(
            || {
                let ts = accum_ref.map(|_| std::time::Instant::now());
                let reader = IndexedReader::from_path(bam_path)
                    .ok()
                    .and_then(|mut b| b.set_threads(1).ok().map(|_| b));
                if let (Some(a), Some(t)) = (accum_ref, ts) {
                    a.bam_open_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                reader
            },
            |bam_opt, (tid, contig_id, ref_name, ref_length)| {
                let bam = bam_opt.as_mut()?;
                let tid = *tid;
                let contig_id = *contig_id;
                let ref_length = *ref_length;

                // Sub-phase: CDS index build
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

                // Sub-phase: process_contig_streaming
                let ts = accum_ref.map(|_| std::time::Instant::now());
                let mut phase_t: Option<BamPhaseTimings> = accum_ref.map(|_| BamPhaseTimings::default());
                let (mut arrays, coverage_pct, primary_count) = match process_contig_streaming(bam, tid, ref_name, ref_length, seq_type, flags, is_circular, config.min_aligned_fraction, config.phagetermini_config.min_clipping_length, &mut phase_t) {
                    Ok(Some(result)) => result,
                    Ok(None) => {
                        // Flush phase_timings even on early-return so coverage-filtered
                        // contigs still show up in the log.
                        if let (Some(a), Some(pt)) = (accum_ref, phase_t.as_ref()) {
                            a.stream_fetch_seek_ns.fetch_add(pt.fetch_seek_ns, Ordering::Relaxed);
                            a.stream_records_ns.fetch_add(pt.records_ns, Ordering::Relaxed);
                            a.stream_process_read_ns.fetch_add(pt.process_read_ns, Ordering::Relaxed);
                            a.stream_finalize_ns.fetch_add(pt.finalize_ns, Ordering::Relaxed);
                        }
                        if let (Some(a), Some(t)) = (accum_ref, ts) {
                            a.streaming_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
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

                // Early-exit: skip contigs with low mean coverage depth (only when threshold > 0)
                let coverage_mean = arrays.coverage_mean() as f32;
                if config.min_coverage_depth > 0.0 && (coverage_mean as f64) < config.min_coverage_depth {
                    return None;
                }

                // Sub-phase: feature calculation
                let ts = accum_ref.map(|_| std::time::Instant::now());
                let mut feature_blobs: Vec<(String, String, crate::blob::EncodedBlob)> = Vec::new();
                let (packaging_info, metrics_info) = add_features_from_arrays(&mut arrays, ref_name, ref_length, config, is_circular, seq_type, flags, primary_count, repeats, cds_index.as_ref(), &mut feature_blobs);
                let coverage_median = arrays.coverage_median() as f32;
                let coverage_trimmed_mean = arrays.coverage_trimmed_mean(0.05) as f32;

                // Above expected aligned fraction: exact_af >= (1 - e^(-0.883 × coverage_mean)) × 100
                let above_expected = coverage_pct >= (1.0 - (-0.883 * coverage_mean as f64).exp()) * 100.0;

                // Calculate coverage variation using Fano factor style normalization
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

                // Coverage SD: Coefficient of Variation (CV) = std_dev / mean
                let coverage_coefficient_of_variation = if arrays.primary_reads.len() > 0 && coverage_mean > 0.0 {
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
            },
        )
        .filter_map(|x| x)
        .collect();

    if let (Some(s), Some(t)) = (st.as_mut(), t_par) {
        s.par_iter_wall = t.elapsed();
        s.n_contigs_processed = results.len();
    }

    // --- Phase: merge_results ---
    let t0 = if time_on { Some(std::time::Instant::now()) } else { None };

    // Merge results from all contigs
    let mut all_feature_blobs: Vec<(String, String, crate::blob::EncodedBlob)> = Vec::new();
    let mut all_presences = Vec::new();
    let mut all_packaging = Vec::new();
    let mut all_misassembly = Vec::new();
    let mut all_microdiversity = Vec::new();
    let mut all_side_misassembly = Vec::new();
    let mut all_topology = Vec::new();
    let mut mapped_reads: u64 = 0;

    for (blobs, presence, packaging, metrics, primary_count) in results {
        all_feature_blobs.extend(blobs);
        all_presences.push(presence);
        mapped_reads += primary_count;
        if let Some(pkg) = packaging {
            all_packaging.push(pkg);
        }
        if let Some((mis, micro, side, topo)) = metrics {
            all_misassembly.push(mis);
            all_microdiversity.push(micro);
            all_side_misassembly.push(side);
            all_topology.push(topo);
        }
    }

    if let (Some(s), Some(t)) = (st.as_mut(), t0) {
        s.merge_results = t.elapsed();
        if let Some(a) = accum.as_ref() {
            s.bam_open_ns = a.bam_open_ns.load(Ordering::Relaxed);
            s.cds_index_ns = a.cds_index_ns.load(Ordering::Relaxed);
            s.streaming_ns = a.streaming_ns.load(Ordering::Relaxed);
            s.feature_calc_ns = a.feature_calc_ns.load(Ordering::Relaxed);
            s.stream_fetch_seek_ns = a.stream_fetch_seek_ns.load(Ordering::Relaxed);
            s.stream_records_ns = a.stream_records_ns.load(Ordering::Relaxed);
            s.stream_process_read_ns = a.stream_process_read_ns.load(Ordering::Relaxed);
            s.stream_finalize_ns = a.stream_finalize_ns.load(Ordering::Relaxed);
        }
    }

    Ok((all_feature_blobs, all_presences, all_packaging, all_misassembly, all_microdiversity, all_side_misassembly, all_topology, sample_name, seq_type, total_reads, mapped_reads, is_circular, st))
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

    let (mut contigs, mut annotations, contig_qualifiers) = if !has_genbank {
        eprintln!("No GenBank file provided - extracting contigs from BAM headers");
        let contigs = extract_contigs_from_bams(bam_files, &preliminary_map)?;
        eprintln!("Found {} contigs from BAM files", contigs.len());
        (contigs, Vec::new(), Vec::new())
    } else {
        parse_annotations(genbank_path)?
    };
    eprintln!(
        "Found {} contigs with {} annotations",
        contigs.len(),
        annotations.len()
    );

    eprintln!("Found {} BAM files", bam_files.len());

    // 3. If assembly_path provided, parse FASTA and merge sequences into contigs
    if !assembly_path.as_os_str().is_empty() && assembly_path.exists() {
        eprintln!("\n### Merging sequences from assembly FASTA...");
        merge_sequences_from_fasta(&mut contigs, assembly_path)?;
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
    let (db_writer, new_contigs_for_blast) = if is_extending {
        // Read existing contig names from the database
        let existing_conn = duckdb::Connection::open(extend_db)
            .with_context(|| format!("Failed to open existing database: {}", extend_db.display()))?;

        let mut existing_contig_names: HashSet<String> = HashSet::new();
        let mut existing_contigs: Vec<ContigInfo> = Vec::new();
        {
            let mut stmt = existing_conn.prepare("SELECT Contig_name, Contig_length FROM Contig")?;
            let rows = stmt.query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?))
            })?;
            for row in rows {
                let (name, length) = row?;
                existing_contig_names.insert(name.clone());
                existing_contigs.push(ContigInfo {
                    name,
                    length: length as usize,
                    sequence: None,
                });
            }
        }

        // Detect if DB was created with annotations (genbank mode)
        let db_has_annotations: bool = existing_conn
            .query_row("SELECT COUNT(*) FROM Contig_annotation", [], |row| row.get::<_, i64>(0))
            .unwrap_or(0) > 0;

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
        annotations = new_annotations.clone();

        let writer = DbWriter::open(extend_db, &new_contigs, &new_annotations, &new_contig_qualifiers, !bam_files.is_empty())?;
        writer.update_metadata_modification()?;
        (writer, new_contigs)
    } else {
        let new_contigs_for_blast = contigs.clone();
        let writer = DbWriter::create(output_db, &contigs, &annotations, &contig_qualifiers, !bam_files.is_empty())?;
        writer.write_metadata(
            modules,
            config.min_aligned_fraction,
            config.min_coverage_depth,
            config.curve_ratio,
            config.bar_ratio,
        )?;
        (writer, new_contigs_for_blast)
    };

    // 5. Auto-detect circularity per sample (now with contigs available for length comparison)
    let circularity_map = if !bam_files.is_empty() {
        detect_all_sample_circularities(bam_files, &contigs)?
    } else {
        HashMap::new()
    };

    // 6. If sequences available, run autoblast with rayon (only for new contigs in extend mode)
    let blast_contigs = &new_contigs_for_blast;
    let has_sequences = blast_contigs.iter().any(|c| c.sequence.is_some());
    let repeats = if has_sequences {
        let reps = run_autoblast(blast_contigs, config.threads)?;
        db_writer.write_repeats(&reps)?;
        reps
    } else {
        Vec::new()
    };

    // 7. Compute and write GC content and GC skew from sequence data (only for new contigs)
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

    let result = process_samples_parallel(&bam_files, &contigs, modules, config, &circularity_map, db_writer, &repeats, &annotations, output_db)?;

    print_summary(&result, output_db);

    Ok(result)
}

/// Holds processed sample data ready for database writing.
struct SampleResult {
    sample_name: String,
    sequencing_type: SequencingType,
    total_reads: u64,
    mapped_reads: u64,
    is_circular: bool,
    /// Compressed BLOB data: Vec<(feature_name, contig_name, encoded_blob)>
    feature_blobs: Vec<(String, String, crate::blob::EncodedBlob)>,
    presences: Vec<PresenceData>,
    packaging: Vec<PackagingData>,
    misassembly: Vec<MisassemblyData>,
    microdiversity: Vec<MicrodiversityData>,
    side_misassembly: Vec<SideMisassemblyData>,
    topology: Vec<TopologyData>,
    timings: Option<SampleTimings>,
}

fn process_samples_parallel(
    bam_files: &[PathBuf],
    contigs: &[ContigInfo],
    modules: &[String],
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    repeats: &[RepeatsData],
    annotations: &[crate::types::FeatureAnnotation],
    output_db: &Path,
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
        return process_samples_sequential(bam_files, contigs, modules, config, circularity_map, db_writer, is_tty, repeats, annotations, output_db);
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

    // Open timing log incrementally so data is available even if calculate crashes
    let timing_file = if config.enable_timing {
        match open_timing_log(output_db, config.threads, total) {
            Ok(f) => Some(f),
            Err(e) => { eprintln!("Warning: could not open timing log: {}", e); None }
        }
    } else {
        None
    };
    let timing_threads = config.threads;

    // Spawn dedicated writer thread
    let writer_handle = thread::spawn(move || -> Result<(usize, std::time::Duration, Vec<SampleTimings>)> {
        let write_start = std::time::Instant::now();
        let mut written_count = 0usize;
        let mut all_timings: Vec<SampleTimings> = Vec::new();
        let mut timing_file = timing_file;

        // Receive and write samples as they arrive
        for result in rx {
            written_count += 1;

            // Insert sample
            if let Err(e) = db_writer.insert_sample(&result.sample_name, result.sequencing_type.as_str(), result.total_reads, result.mapped_reads, result.is_circular) {
                eprintln!("\nError inserting sample {}: {}", result.sample_name, e);
                let msg = format!("ERR: {}", result.sample_name);
                write_pb_clone.set_message(msg.clone());
                write_pb_clone.inc(1);
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                }
                continue;
            }

            // Write all data for this sample
            let sample_write_start = std::time::Instant::now();
            if let Err(e) = db_writer.write_sample_data(
                &result.sample_name,
                &result.presences,
                &result.packaging,
                &result.misassembly,
                &result.microdiversity,
                &result.side_misassembly,
                &result.topology,
                &result.feature_blobs,
                result.is_circular,
            ) {
                eprintln!("\nError writing data for {}: {}", result.sample_name, e);
                let msg = format!("ERR: {}", result.sample_name);
                write_pb_clone.set_message(msg.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, msg);
                }
            } else {
                write_pb_clone.set_message(result.sample_name.clone());
                if !is_tty_writer {
                    eprintln!("Writing:    [{}/{}] {}", written_count, total_writer, result.sample_name);
                }
            }
            write_pb_clone.inc(1);

            // Write per-sample timings incrementally to the log file
            if let Some(mut t) = result.timings {
                t.db_write = sample_write_start.elapsed();
                if let Some(ref mut f) = timing_file {
                    let _ = write_sample_timing(f, &t, timing_threads);
                }
                all_timings.push(t);
            }
        }

        // Finalize database
        db_writer.finalize()?;
        write_pb_clone.finish();

        Ok((written_count, write_start.elapsed(), all_timings))
    });

    // Process samples in parallel, sending to channel immediately
    // send() blocks if channel is full (backpressure)
    bam_files.par_iter().for_each(|bam_path| {
        let sample_start = std::time::Instant::now();

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        match process_sample(bam_path, contigs, modules, config, sample_circular, repeats, annotations) {
            Ok((feature_blobs, presences, packaging, misassembly, microdiversity, side_misassembly, topology, sample_name, sequencing_type, total_reads, mapped_reads, is_circular, mut timings)) => {
                let sample_time = sample_start.elapsed().as_secs_f64();
                // Fill the total processing time (DB write will be added by writer thread)
                if let Some(ref mut t) = timings {
                    t.total = sample_start.elapsed();
                }
                completed_count.fetch_add(1, Ordering::SeqCst);
                let msg = format!("{} ({:.2}s)", sample_name, sample_time);
                process_pb.set_message(msg.clone());
                process_pb.inc(1);
                if !is_tty {
                    eprintln!("Processing: {}", msg);
                }

                // Send to writer thread (blocks if channel full)
                let _ = tx.send(SampleResult {
                    sample_name,
                    sequencing_type,
                    total_reads,
                    mapped_reads,
                    is_circular,
                    feature_blobs,
                    presences,
                    packaging,
                    misassembly,
                    microdiversity,
                    side_misassembly,
                    topology,
                    timings,
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
    let (written_count, writing_time, all_timings) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    // Clear MultiProgress to avoid duplicate bar display
    mp.clear().ok();

    let elapsed = start_time.elapsed();
    let failed = failed_count.load(Ordering::SeqCst);

    // Append grand totals to the timing log
    if config.enable_timing && !all_timings.is_empty() {
        let log_path = timing_log_path(output_db);
        match fs::OpenOptions::new().append(true).open(&log_path) {
            Ok(mut f) => {
                if let Err(e) = finish_timing_log(&mut f, &all_timings, elapsed) {
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
    modules: &[String],
    config: &ProcessConfig,
    circularity_map: &HashMap<PathBuf, bool>,
    db_writer: DbWriter,
    is_tty: bool,
    repeats: &[RepeatsData],
    annotations: &[crate::types::FeatureAnnotation],
    output_db: &Path,
) -> Result<ProcessResult> {
    let total = bam_files.len();
    let start_time = std::time::Instant::now();

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

    // Open timing log incrementally so data is available even if calculate crashes
    let mut timing_file = if config.enable_timing {
        match open_timing_log(output_db, config.threads, total) {
            Ok(f) => Some(f),
            Err(e) => { eprintln!("Warning: could not open timing log: {}", e); None }
        }
    } else {
        None
    };

    for bam_path in bam_files {
        let sample_start = std::time::Instant::now();

        let sample_circular = circularity_map.get(bam_path).copied().unwrap_or(false);
        match process_sample(bam_path, contigs, modules, config, sample_circular, repeats, annotations) {
            Ok((feature_blobs, presences, packaging, misassembly, microdiversity, side_misassembly, topology, sample_name, sequencing_type, total_reads, mapped_reads, is_circular, mut timings)) => {
                let process_elapsed = sample_start.elapsed();
                processing_time_total += process_elapsed;

                let write_start = std::time::Instant::now();

                // Insert sample
                if let Err(e) = db_writer.insert_sample(&sample_name, sequencing_type.as_str(), total_reads, mapped_reads, is_circular) {
                    eprintln!("\nError inserting sample {}: {}", sample_name, e);
                    failed += 1;
                    pb.inc(1);
                    continue;
                }

                // Write sample data
                if let Err(e) = db_writer.write_sample_data(
                    &sample_name,
                    &presences,
                    &packaging,
                    &misassembly,
                    &microdiversity,
                    &side_misassembly,
                    &topology,
                    &feature_blobs,
                    is_circular,
                ) {
                    eprintln!("\nError writing data for {}: {}", sample_name, e);
                    failed += 1;
                } else {
                    processed += 1;
                }

                writing_time_total += write_start.elapsed();

                // Write per-sample timings incrementally to the log file
                if let Some(ref mut t) = timings {
                    t.db_write = write_start.elapsed();
                    t.total = sample_start.elapsed();
                    if let Some(ref mut f) = timing_file {
                        let _ = write_sample_timing(f, t, config.threads);
                    }
                    all_timings.push(t.clone());
                }

                let total_sample_time = sample_start.elapsed().as_secs_f64();
                let msg = format!("{} ({:.2}s)", sample_name, total_sample_time);
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

    // Finalize database
    db_writer.finalize()?;
    pb.finish_with_message("Done");
    if !is_tty {
        eprintln!("Done ({:.2}s)", start_time.elapsed().as_secs_f64());
    }

    // Append grand totals to the timing log
    if config.enable_timing && !all_timings.is_empty() {
        if let Some(ref mut f) = timing_file {
            if let Err(e) = finish_timing_log(f, &all_timings, start_time.elapsed()) {
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
    writeln!(f, "  ---- Par-iter sub-phase sums (CPU-seconds across threads) ----")?;
    writeln!(f, "  BAM reader open      : {:>10.3} s", ns_to_s(t.bam_open_ns))?;
    writeln!(f, "  CDS index build      : {:>10.3} s", ns_to_s(t.cds_index_ns))?;
    writeln!(f, "  process_contig_stream: {:>10.3} s", ns_to_s(t.streaming_ns))?;
    writeln!(f, "    . fetch seek       : {:>10.3} s", ns_to_s(t.stream_fetch_seek_ns))?;
    writeln!(f, "    . record loop (htslib+prep): {:>10.3} s", ns_to_s(t.stream_records_ns))?;
    writeln!(f, "    . process_read     : {:>10.3} s", ns_to_s(t.stream_process_read_ns))?;
    writeln!(f, "    . finalize+cov     : {:>10.3} s", ns_to_s(t.stream_finalize_ns))?;
    writeln!(f, "  Feature calculation  : {:>10.3} s", ns_to_s(t.feature_calc_ns))?;
    let sub_sum_ns = t.bam_open_ns + t.cds_index_ns
                   + t.streaming_ns + t.feature_calc_ns;
    writeln!(f, "  (sum of above)       : {:>10.3} s  ({}x threads -> wall ~ {:.1} s)",
             ns_to_s(sub_sum_ns), threads, ns_to_s(sub_sum_ns) / threads as f64)?;
    writeln!(f)?;
    Ok(())
}

/// Append grand totals to the timing log and close it.
fn finish_timing_log(
    f: &mut fs::File,
    timings: &[SampleTimings],
    total_wall: std::time::Duration,
) -> Result<()> {
    use std::io::Write;
    writeln!(f, "=============================================================")?;
    writeln!(f, "GRAND TOTALS")?;
    writeln!(f, "=============================================================")?;
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
