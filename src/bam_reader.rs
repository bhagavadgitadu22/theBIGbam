//! BAM file reading and sequencing type detection.
//!
//! Processes BAM files using single-pass streaming to avoid loading all reads into memory.
//! Detects sequencing type (short paired, short single, or long reads) from file.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::path::Path;

use crate::features::{process_read, FeatureArrays, ModuleFlags};
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
/// Returns FeatureArrays with calculated features, or None if contig has no reads.
pub fn process_contig_streaming(
    bam: &mut bam::IndexedReader,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_coverage: f64,
) -> Result<Option<(FeatureArrays, f64)>> {
    // -------------------------------------------------------------------------
    // Step 1: Check if this contig exists in the BAM file
    // -------------------------------------------------------------------------
    // BAM files have a header listing all reference sequences. If our contig
    // isn't in the header, there are no reads for it.
    if bam.header().tid(contig_name.as_bytes()).is_none() {
        return Ok(None);
    }

    bam.fetch((contig_name, 0, ref_length as i64 * 2))
        .with_context(|| format!("Failed to fetch reads for contig: {}", contig_name))?;

    let mut arrays = FeatureArrays::new(ref_length);
    let need_md = flags.needs_md();
    let mut has_reads = false;
    let mut cigar_buf: Vec<(u32, u32)> = Vec::with_capacity(16);

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped() {
            continue;
        }

        has_reads = true;

        // Parse CIGAR string and MD tag (if needed)
        let cigar_view = record.cigar();
        cigar_buf.clear();
        cigar_buf.extend(cigar_view.iter().map(|c| (c.char() as u32, c.len())));

        let md_tag: Option<&[u8]> = if need_md {
            record.aux(b"MD").ok().and_then(|aux| match aux {
                rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                _ => None,
            })
        } else {
            None
        };

        // Compute template length and proper_pair with circular correction if needed
        let (template_length, is_proper_pair) = if circular && seq_type.is_short_paired() {
            // In circular mode, insert_size and is_proper_pair from BAM are based on doubled assembly
            // We need to recompute based on the actual circular genome coordinates
            let pos1 = record.pos() as usize;
            let pos2 = record.mpos() as usize;
            
            // Both positions modulo actual genome length
            let p1 = pos1 % ref_length;
            let p2 = pos2 % ref_length;
            
            // Distance could wrap around - take the minimum of direct and wrapped distance
            let direct = (p1 as i32 - p2 as i32).abs();
            let wrapped = ref_length as i32 - direct;
            let corrected_tlen = direct.min(wrapped);
            
            // Recompute proper_pair: pairs are proper if they're on same contig, opposite strands,
            // facing each other, and within reasonable insert size
            let same_ref = record.tid() == record.mtid();
            let opposite_strands = record.is_reverse() != record.is_mate_reverse();
            let reasonable_distance = corrected_tlen > 0 && corrected_tlen < 10000; // typical paired-end insert size range
            let corrected_proper = same_ref && opposite_strands && reasonable_distance;
            
            (corrected_tlen, corrected_proper)
        } else {
            (record.insert_size().abs() as i32, record.is_proper_pair())
        };

        process_read(
            &mut arrays,
            record.pos(),
            cigar_view.end_pos(),
            record.seq_len() as i32,
            template_length,
            record.is_first_in_template(),
            is_proper_pair,
            record.is_reverse(),
            record.is_secondary(),
            record.is_supplementary(),
            record.is_mate_unmapped(),
            record.tid() != record.mtid(),
            &cigar_buf,
            md_tag,
            record.mapq(),
            seq_type,
            flags,
            circular,
        );
    }

    if !has_reads {
        return Ok(None);
    }

    // Check coverage percentage before post-processing and feature calculation
    // This saves time on low-coverage contigs by skipping:
    // 1. finalize_strands() (phagetermini-specific, expensive)
    // 2. All feature compression in the caller (very expensive)
    // 3. Database writes
    let coverage_pct = arrays.coverage_percentage();
    if coverage_pct < min_coverage {
        return Ok(None);
    }

    // Finalize strand-specific reads_starts/reads_ends for phagetermini
    // This merges forward/reverse strand tracking into final arrays
    // Only needed for phagetermini; other features are already finalized
    if flags.phagetermini {
        arrays.finalize_strands(seq_type);
    }

    Ok(Some((arrays, coverage_pct)))
}
