//! BAM file processing functions.
//!
//! This module handles reading BAM files and extracting read data.
//! Supports both batch processing (legacy) and streaming processing (optimized).

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::path::Path;

use crate::features::{process_read, FeatureArrays, ModuleFlags};
use crate::types::SequencingType;

/// Thresholds for sequencing type detection.
const LONG_READ_LENGTH_THRESHOLD: usize = 1000;
const SEQUENCING_TYPE_SAMPLE_SIZE: usize = 100;

/// Detect sequencing type from first N reads.
pub fn detect_sequencing_type(bam_path: &Path) -> Result<SequencingType> {
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path.display()))?;

    let mut n_checked = 0;

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;

        if record.is_unmapped() {
            continue;
        }

        if record.seq_len() > LONG_READ_LENGTH_THRESHOLD {
            return Ok(SequencingType::Long);
        }

        if record.is_paired() {
            return Ok(SequencingType::ShortPaired);
        }

        n_checked += 1;
        if n_checked >= SEQUENCING_TYPE_SAMPLE_SIZE {
            break;
        }
    }

    Ok(SequencingType::ShortSingle)
}

/// Process reads for a contig using streaming (single-pass, no intermediate storage).
///
/// This is the optimized version that processes reads directly from BAM
/// and updates feature arrays in a single pass.
pub fn process_contig_streaming(
    bam: &mut bam::IndexedReader,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
) -> Result<Option<FeatureArrays>> {
    // Check if contig exists in BAM
    if bam.header().tid(contig_name.as_bytes()).is_none() {
        return Ok(None);
    }

    // Fetch reads for contig (doubled length for circular genomes)
    bam.fetch((contig_name, 0, ref_length as i64 * 2))
        .with_context(|| format!("Failed to fetch reads for contig: {}", contig_name))?;

    let mut arrays = FeatureArrays::new(ref_length);
    let need_md = flags.needs_md();
    let mut has_reads = false;

    // Reusable buffer to avoid per-read allocations
    let mut cigar_buf: Vec<(u32, u32)> = Vec::with_capacity(16);

    // Process reads directly from BAM - no intermediate storage
    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped() {
            continue;
        }

        has_reads = true;

        // Extract CIGAR into reusable buffer
        let cigar_view = record.cigar();
        cigar_buf.clear();
        cigar_buf.extend(cigar_view.iter().map(|c| (c.char() as u32, c.len())));

        // Extract MD tag - borrow directly without allocation
        let md_tag: Option<&[u8]> = if need_md {
            record.aux(b"MD").ok().and_then(|aux| match aux {
                rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                _ => None,
            })
        } else {
            None
        };

        // Process this read immediately - single pass
        process_read(
            &mut arrays,
            record.pos(),
            cigar_view.end_pos(),
            record.seq_len() as i32,
            record.insert_size().abs() as i32,
            record.is_first_in_template(),
            record.is_proper_pair(),
            record.is_reverse(),
            &cigar_buf,
            md_tag,
            seq_type,
            flags,
        );
    }

    if !has_reads {
        return Ok(None);
    }

    // Finalize strand arrays for phagetermini
    if flags.phagetermini {
        arrays.finalize_strands(seq_type);
    }

    Ok(Some(arrays))
}
