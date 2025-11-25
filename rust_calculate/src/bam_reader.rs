//! BAM file processing functions.
//!
//! This module handles reading BAM files and extracting read data.
//! Python equivalent: Functions in calculating_data.py related to BAM processing.

use anyhow::Result;
use rust_htslib::bam::{self, Read as BamRead};
use std::path::Path;

use crate::types::{ReadData, SequencingType};

/// Detect sequencing type from first N reads.
/// Python equivalent: `find_sequencing_type_from_bam()` in calculating_data.py:55-78
pub fn detect_sequencing_type(bam_path: &Path) -> Result<SequencingType> {
    // calculating_data.py:55-78 - find_sequencing_type_from_bam()
    let mut bam = bam::Reader::from_path(bam_path)?;  // py:63
    let mut n_checked = 0;

    for result in bam.records() {  // py:65
        let record = result?;
        if record.is_unmapped() {  // py:66
            continue;
        }

        // py:68 - If any read > 1000 bp, reads are "long"
        if record.seq_len() > 1000 {
            return Ok(SequencingType::Long);  // py:70
        }

        // py:71 - Elif any read.is_paired, reads are "short-paired"
        if record.is_paired() {
            return Ok(SequencingType::ShortPaired);  // py:72
        }

        n_checked += 1;
        if n_checked >= 100 {  // py:74 - n_reads_check=100
            break;
        }
    }

    Ok(SequencingType::ShortSingle)  // py:78 - default to "short-single"
}

/// Extract read data from BAM for one contig.
/// Python equivalent: `preprocess_reads()` in calculating_data.py:497-602
pub fn process_reads_for_contig(
    bam: &mut bam::IndexedReader,
    contig_name: &str,
    contig_length: usize,
    modules: &[String],
    _seq_type: SequencingType,
) -> Result<Vec<ReadData>> {
    // calculating_data.py:497-602 - preprocess_reads()
    let tid = bam.header().tid(contig_name.as_bytes());
    if tid.is_none() {
        return Ok(Vec::new());
    }

    // py:609 - reads_mapped = bam_file.fetch(ref_name)
    bam.fetch((contig_name, 0, contig_length as i64 * 2))?;

    // py:500-501 - Determine which attributes we need based on modules
    let need_md = modules.contains(&"phagetermini".to_string())
        || modules.contains(&"assemblycheck".to_string());

    let mut reads = Vec::new();

    // py:529-578 - Single pass through reads
    for result in bam.records() {
        let record = result?;
        if record.is_unmapped() {  // py:530
            continue;
        }

        // py:565-566 - Extract CIGAR tuples
        let cigar: Vec<(u32, u32)> = record
            .cigar()
            .iter()
            .map(|c| (c.char() as u32, c.len()))
            .collect();

        // py:567-576 - Extract MD tag if needed
        let md_tag = if need_md {
            record.aux(b"MD").ok().and_then(|aux| {
                match aux {
                    rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes().to_vec()),
                    _ => None,
                }
            })
        } else {
            None
        };

        // py:551-564 - Extract read attributes
        reads.push(ReadData {
            ref_start: record.pos(),             // py:552 - read.reference_start
            ref_end: record.cigar().end_pos(),   // py:553 - read.reference_end
            query_length: record.seq_len() as i32,        // py:557 - read.query_length
            template_length: record.insert_size().abs() as i32,  // py:559 - abs(read.template_length)
            is_read1: record.is_first_in_template(),      // py:560 - read.is_read1
            is_proper_pair: record.is_proper_pair(),      // py:561 - read.is_proper_pair
            is_reverse: record.is_reverse(),              // py:564 - read.is_reverse
            cigar,
            md_tag,
        });
    }

    Ok(reads)
}
