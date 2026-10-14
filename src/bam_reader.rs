//! BAM file reading and sequencing type detection.
//!
//! Processes BAM files using single-pass streaming to avoid loading all reads into memory.
//! Detects sequencing type (short paired, short single, or long reads) from file.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::path::Path;

use crate::features::{FeatureArrays, ModuleFlags};

mod accumulator;
pub(crate) use accumulator::ContigAccumulator;
use crate::types::SequencingType;

// ============================================================================
// Constants for Sequencing Type Detection
// ============================================================================

/// Reads longer than this are considered "long reads" (PacBio/Nanopore).
/// Illumina reads are typically 75-300bp, while PacBio/Nanopore are 1000-50000bp.
const LONG_READ_LENGTH_THRESHOLD: usize = 1000;

/// How many reads to examine when detecting sequencing type.
/// We only need to check a few reads - they should all be the same type.
const SEQUENCING_TYPE_SAMPLE_SIZE: usize = 100;

/// Distance from contig ends (in bp) within which a read is considered "near the edge".
/// Used for circularity confirmation and mate-pair topology.
const NEAR_CONTIG_END_DISTANCE: usize = 1000;

/// Maximum reasonable insert size (bp) for paired-end reads.
/// Pairs with template length exceeding this are treated as structural anomalies.
const MAX_REASONABLE_INSERT_SIZE: i32 = 10000;

/// Per-contig timing breakdown for `process_contig_streaming`.
/// Populated only when the caller passes `Some(&mut _)`.
#[derive(Default, Clone, Debug)]
pub struct BamPhaseTimings {
    /// `bam.fetch(...)` — BAM index seek cost.
    pub fetch_seek_ns: u64,
    /// Record-loop wall minus the time spent inside `process_read`. This is
    /// htslib record decoding + our per-record filtering (cigar/md/seq).
    pub records_ns: u64,
    /// Cumulative time inside `process_read(...)` — our feature-array work.
    pub process_read_ns: u64,
    /// Post-loop: `coverage_percentage()` + `finalize_strands()` (if phagetermini).
    pub finalize_ns: u64,
}

// ============================================================================
// Sequencing Type Detection
// ============================================================================

/// Detect the sequencing technology from a BAM file by examining the first few reads.
///
/// Returns Long if any read >1000bp, ShortPaired if any read is paired,
/// or ShortSingle after checking 100 reads.
pub fn detect_sequencing_type(bam_path: &Path) -> Result<SequencingType> {
    // Open the BAM file for reading
    // The `?` operator: if this fails, return the error immediately
    // .with_context() adds helpful information to error messages
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path.display()))?;

    let mut n_checked = 0;

    // Iterate through reads
    // In Rust, `for x in iterator` consumes the iterator one item at a time
    for result in bam.records() {
        // Each record could fail to parse, so we get Result<Record, Error>
        let record = result.context("Failed to read BAM record")?;

        // Skip unmapped reads - they don't tell us about sequencing type
        if record.is_unmapped() {
            continue;
        }

        // Check for long reads (PacBio, Nanopore)
        if record.seq_len() > LONG_READ_LENGTH_THRESHOLD {
            return Ok(SequencingType::Long);
        }

        // Check for paired-end reads (Illumina paired)
        if record.is_paired() {
            return Ok(SequencingType::ShortPaired);
        }

        // Count how many reads we've checked
        n_checked += 1;
        if n_checked >= SEQUENCING_TYPE_SAMPLE_SIZE {
            break;
        }
    }

    // If we get here, reads are short and not paired → single-end
    Ok(SequencingType::ShortSingle)
}

// ============================================================================
// Main Processing Function - Streaming Approach
// ============================================================================

/// Process all reads for a contig using single-pass streaming.
///
/// Returns (FeatureArrays, coverage_pct, primary_count) or None if contig has no reads.
pub fn process_contig_streaming(
    bam: &mut bam::IndexedReader,
    tid: u32,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_aligned_fraction: f64,
    min_clipping_length: u32,
    phase_timings: &mut Option<BamPhaseTimings>,
) -> Result<Option<(FeatureArrays, f64, u64)>> {
    process_contig_streaming_counted(bam,tid,contig_name,ref_length,seq_type,flags,circular,min_aligned_fraction,min_clipping_length,phase_timings).map(|(result,_)|result)
}

pub(crate) fn process_contig_streaming_counted(
    bam: &mut bam::IndexedReader,
    tid: u32,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_aligned_fraction: f64,
    min_clipping_length: u32,
    phase_timings: &mut Option<BamPhaseTimings>,
) -> Result<(Option<(FeatureArrays, f64, u64)>, u64)> {
    BAM_SCRATCH.with(|slot| {
        let mut buffers = slot.borrow_mut();
        let result = process_contig_with_buffers(bam, tid, contig_name, ref_length, seq_type, flags,
            circular, min_aligned_fraction, min_clipping_length, phase_timings, &mut buffers);
        // Reuse record storage without pinning the sample's huge BAM header.
        let empty_header = buffers.empty_header.clone();
        buffers.record.set_header(empty_header);
        result
    })
}

struct BamBuffers {
    record: bam::Record,
    cigar_buf: Vec<(u32, u32)>,
    seq_buf: Vec<u8>,
    empty_header: std::sync::Arc<bam::HeaderView>,
}

impl Default for BamBuffers {
    fn default() -> Self {
        Self {
            record: bam::Record::new(), cigar_buf: Vec::with_capacity(16), seq_buf: Vec::with_capacity(256),
            empty_header: std::sync::Arc::new(bam::HeaderView::from_header(&bam::Header::new())),
        }
    }
}

thread_local! {
    static BAM_SCRATCH: std::cell::RefCell<BamBuffers> =
        std::cell::RefCell::new(BamBuffers::default());
}

/// Decode BAM's packed CIGAR directly into reusable scratch, without the
/// intermediate CigarString allocation made by Record::cigar().
pub(crate) fn prepare_cigar(raw: &[u32], start: i64, scratch: &mut Vec<(u32, u32)>) -> i64 {
    const OPS: &[u8; 9] = b"MIDNSHP=X";
    scratch.clear();
    let mut end = start;
    for &packed in raw {
        let code = (packed & 15) as usize;
        let op = *OPS.get(code).expect("Unexpected cigar operation") as u32;
        let len = packed >> 4;
        scratch.push((op, len));
        if crate::cigar::raw_cigar_consumes_ref(op) { end += len as i64; }
    }
    end
}

fn process_contig_with_buffers(
    bam: &mut bam::IndexedReader,
    tid: u32,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_aligned_fraction: f64,
    min_clipping_length: u32,
    phase_timings: &mut Option<BamPhaseTimings>,
    buffers: &mut BamBuffers,
) -> Result<(Option<(FeatureArrays, f64, u64)>, u64)> {
    let time_on = phase_timings.is_some();

    // Fetch by tid directly — skips the name→tid hash lookup (and, crucially,
    // the lazy O(n_refs) hash-table build the first time any name lookup runs
    // on a freshly opened reader). Callers pre-validate the tid against the
    // BAM header, so no existence check is needed here.
    let t_fetch = if time_on { Some(std::time::Instant::now()) } else { None };
    bam.fetch(bam::FetchDefinition::Region(tid as i32, 0, ref_length as i64))
        .with_context(|| format!("Failed to fetch reads for contig: {}", contig_name))?;
    if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_fetch) {
        pt.fetch_seek_ns = t.elapsed().as_nanos() as u64;
    }

    let mut accumulator = ContigAccumulator::new(ref_length, seq_type, flags, circular,
        min_clipping_length, time_on);
    let BamBuffers { record, cigar_buf, seq_buf, .. } = buffers;
    let t_loop = time_on.then(std::time::Instant::now);
    while let Some(result) = bam.read(record) {
        result.with_context(||format!("Reading BAM records for {contig_name}"))?;
        accumulator.consume_record(record, cigar_buf, seq_buf);
    }
    if let Some(t) = t_loop {
        accumulator.record_loop_timing(t.elapsed().as_nanos() as u64, phase_timings);
    }
    let count=accumulator.primary_count();
    Ok((accumulator.finish(min_aligned_fraction, phase_timings),count))
}

// ============================================================================
// Read Count Statistics
// ============================================================================

/// Get unmapped read count from BAM index (fast, no iteration needed).
///
/// Returns the number of unmapped records (sum of the `unmapped` field across
/// all entries of `index_stats()`, including the tid=-1 entry for truly
/// unmapped reads). Unmapped reads have no secondary/supplementary alignments,
/// so one record == one distinct read.
///
/// The caller composes the true total as:
///     total_reads = mapped_reads (primary only) + unmapped_reads
/// which deduplicates multi-alignment mapped reads — unlike the BAM index's
/// `mapped` column, which counts alignment records (primary + secondary +
/// supplementary).
pub fn get_unmapped_read_count(bam_path: &Path) -> Result<u64> {
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;

    let stats = bam.index_stats()
        .with_context(|| format!("Failed to get index stats from: {}", bam_path.display()))?;

    let unmapped_total: u64 = stats.iter().map(|(_, _, _mapped, unmapped)| *unmapped).sum();

    Ok(unmapped_total)
}

#[cfg(test)]
mod buffer_tests {
    use super::*;

    #[test]
    fn packed_cigar_preparation_matches_htslib_view() {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let cigar = CigarString(vec![Cigar::HardClip(2), Cigar::SoftClip(5), Cigar::Match(10),
            Cigar::Ins(3), Cigar::Equal(4), Cigar::Diff(2), Cigar::Del(5), Cigar::RefSkip(20), Cigar::Pad(1)]);
        let mut record = bam::Record::new();
        record.set(b"read", Some(&cigar), &vec![b'A'; 24], &vec![30; 24]);
        record.set_pos(7);
        let mut scratch = vec![(999, 999)];
        let end = prepare_cigar(record.raw_cigar(), record.pos(), &mut scratch);
        let view = record.cigar();
        assert_eq!(end, view.end_pos());
        assert_eq!(scratch, view.iter().map(|c| (c.char() as u32, c.len())).collect::<Vec<_>>());
        assert_eq!(prepare_cigar(&[], 11, &mut scratch), 11);
        assert!(scratch.is_empty());
    }
}

pub(crate) fn count_reference_primary_reads(bam:&mut bam::IndexedReader,tid:u32,length:u64)->Result<u64> {
    bam.fetch((tid as i32,0,length as i64))?;
    let mut record=bam::Record::new(); let mut count=0;
    while let Some(result)=bam.read(&mut record) { result?; if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {count+=1;} }
    Ok(count)
}
