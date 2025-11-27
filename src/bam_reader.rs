//! BAM file processing functions.
//!
//! This module handles reading BAM files and extracting read data for feature calculation.
//!
//! # Bioinformatics context:
//!
//! BAM (Binary Alignment Map) files contain sequencing reads aligned to a reference genome.
//! Each read has:
//! - Position on the reference (contig name, start position)
//! - CIGAR string describing how the read aligns (matches, insertions, deletions, clips)
//! - MD tag describing mismatches
//! - Flags (paired, proper pair, reverse strand, etc.)
//!
//! # Architecture - Streaming Processing:
//!
//! Instead of loading all reads into memory (which could be millions of reads),
//! we use "streaming" - processing each read immediately as we read it from the BAM file.
//! This is much more memory-efficient:
//!
//! ```text
//! OLD APPROACH (batch):
//!   Read all BAM records → Store in Vec<ReadData> → Process all reads
//!   Memory usage: O(number of reads) - could be gigabytes!
//!
//! NEW APPROACH (streaming):
//!   For each BAM record: Read → Process immediately → Discard
//!   Memory usage: O(1) - constant, regardless of read count!
//! ```
//!
//! # For Python developers:
//!
//! ## Key Rust concepts in this file:
//!
//! - **`Result<T, E>`**: Rust's way of handling errors. Instead of raising exceptions,
//!   functions return `Ok(value)` on success or `Err(error)` on failure.
//!   The `?` operator propagates errors automatically (like Python's `raise`).
//!
//! - **`&mut`**: A mutable reference. Allows modifying the borrowed value.
//!   In Python everything is mutable by default; Rust requires explicit marking.
//!
//! - **`anyhow::Result`**: A convenient error type that can hold any error.
//!   Similar to Python's generic `Exception`.

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
/// # Algorithm:
/// 1. Open the BAM file
/// 2. Look at each read until we find evidence of the sequencing type:
///    - If any read is >1000bp → Long reads (PacBio/Nanopore)
///    - If any read is paired → Short paired (Illumina paired-end)
///    - After 100 reads with neither → Short single (Illumina single-end)
///
/// # Arguments
/// * `bam_path` - Path to the BAM file (must be sorted and indexed)
///
/// # Returns
/// * `Ok(SequencingType)` - The detected sequencing type
/// * `Err(...)` - If the file couldn't be opened or read
///
/// # Python equivalent:
/// ```python
/// def find_sequencing_type_from_bam(bam_path):
///     with pysam.AlignmentFile(bam_path, "rb") as bam:
///         for i, read in enumerate(bam):
///             if read.is_unmapped:
///                 continue
///             if read.query_length > 1000:
///                 return "long"
///             if read.is_paired:
///                 return "short_paired"
///             if i >= 100:
///                 break
///     return "short_single"
/// ```
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

/// Process all reads for a contig using streaming (single-pass, no intermediate storage).
///
/// This is the core optimized function. For each read in the BAM file:
/// 1. Read the BAM record from disk
/// 2. Extract needed fields (position, CIGAR, MD tag, flags)
/// 3. Update feature arrays immediately
/// 4. Discard the record (don't store it)
///
/// # Why streaming is better:
///
/// For a contig with 1 million reads, the old approach would:
/// - Allocate memory for 1 million ReadData structs (~200 bytes each = 200MB)
/// - Store all reads, then iterate through them multiple times
///
/// The streaming approach:
/// - Only keeps one read in memory at a time (~200 bytes)
/// - Processes everything in a single pass through the BAM file
///
/// # Arguments
/// * `bam` - The indexed BAM reader (must be positioned at the contig)
/// * `contig_name` - Name of the contig to process (e.g., "phage_lambda")
/// * `ref_length` - Length of the reference sequence in base pairs
/// * `seq_type` - Sequencing type (affects which features to calculate)
/// * `flags` - Which modules are enabled (coverage, phagetermini, assemblycheck)
///
/// # Returns
/// * `Ok(Some(FeatureArrays))` - Arrays containing all calculated features
/// * `Ok(None)` - If the contig doesn't exist in the BAM or has no reads
/// * `Err(...)` - If there was an error reading the BAM file
///
/// # Performance notes:
/// - We reuse a single buffer for CIGAR data to avoid allocations per read
/// - MD tag is borrowed directly from the record (zero-copy)
/// - All feature calculations happen in `process_read()` in one pass
pub fn process_contig_streaming(
    bam: &mut bam::IndexedReader,
    contig_name: &str,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
) -> Result<Option<FeatureArrays>> {
    // -------------------------------------------------------------------------
    // Step 1: Check if this contig exists in the BAM file
    // -------------------------------------------------------------------------
    // BAM files have a header listing all reference sequences. If our contig
    // isn't in the header, there are no reads for it.
    if bam.header().tid(contig_name.as_bytes()).is_none() {
        return Ok(None);
    }

    // -------------------------------------------------------------------------
    // Step 2: Fetch reads for this contig
    // -------------------------------------------------------------------------
    // The BAM index (.bai file) allows us to jump directly to reads for this contig.
    // We fetch reads from position 0 to 2x the reference length because the reference
    // was doubled to handle circular genomes (reads can wrap around).
    bam.fetch((contig_name, 0, ref_length as i64 * 2))
        .with_context(|| format!("Failed to fetch reads for contig: {}", contig_name))?;

    // -------------------------------------------------------------------------
    // Step 3: Initialize feature arrays
    // -------------------------------------------------------------------------
    // FeatureArrays holds all the per-position counters we're calculating:
    // coverage, reads_starts, reads_ends, insertions, deletions, etc.
    // Each array has one slot per base pair in the reference.
    let mut arrays = FeatureArrays::new(ref_length);

    // Check if any module needs the MD tag (mismatches)
    // The MD tag is optional in BAM files, and extracting it has some cost
    let need_md = flags.needs_md();

    // Track if we found any reads (empty contigs return None)
    let mut has_reads = false;

    // -------------------------------------------------------------------------
    // Step 4: Pre-allocate reusable buffer
    // -------------------------------------------------------------------------
    // OPTIMIZATION: Instead of allocating a new Vec for each read's CIGAR,
    // we reuse one buffer. This avoids millions of small allocations.
    //
    // Rust concept - Vec::with_capacity():
    // Pre-allocates space for 16 elements. Most CIGARs have <16 operations.
    // If we need more, the Vec will grow automatically, but that's rare.
    let mut cigar_buf: Vec<(u32, u32)> = Vec::with_capacity(16);

    // -------------------------------------------------------------------------
    // Step 5: Process reads one at a time (streaming)
    // -------------------------------------------------------------------------
    // This loop reads each BAM record, processes it, and lets it be discarded.
    // Memory usage stays constant regardless of how many reads there are.
    for result in bam.records() {
        // Handle potential read errors (corrupted BAM, etc.)
        let record = match result {
            Ok(r) => r,
            Err(_) => continue, // Skip bad records
        };

        // Skip unmapped reads - they don't contribute to coverage
        if record.is_unmapped() {
            continue;
        }

        has_reads = true;

        // ---------------------------------------------------------------------
        // Extract CIGAR operations into our reusable buffer
        // ---------------------------------------------------------------------
        // CIGAR describes how the read aligns: matches (M), insertions (I),
        // deletions (D), soft clips (S), etc.
        //
        // We convert to (operation_char, length) tuples:
        // "10M2I30M" → [(M, 10), (I, 2), (M, 30)]
        let cigar_view = record.cigar();
        cigar_buf.clear(); // Reuse buffer - much faster than allocating new Vec
        cigar_buf.extend(cigar_view.iter().map(|c| (c.char() as u32, c.len())));

        // ---------------------------------------------------------------------
        // Extract MD tag (describes mismatches)
        // ---------------------------------------------------------------------
        // The MD tag is optional and only needed for phagetermini and assemblycheck.
        // We borrow directly from the record to avoid copying.
        //
        // Rust concept - Option<&[u8]>:
        // This is an optional reference to a byte slice. If the MD tag exists,
        // we get Some(&bytes); if not, we get None.
        let md_tag: Option<&[u8]> = if need_md {
            record.aux(b"MD").ok().and_then(|aux| match aux {
                rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                _ => None,
            })
        } else {
            None
        };

        // ---------------------------------------------------------------------
        // Process this read immediately - update all feature arrays
        // ---------------------------------------------------------------------
        // This is where the magic happens! All features are calculated in one
        // function call, updating the arrays in place.
        process_read(
            &mut arrays,           // Arrays to update (mutable reference)
            record.pos(),          // Reference start position (0-based)
            cigar_view.end_pos(),  // Reference end position
            record.seq_len() as i32,           // Query (read) length
            record.insert_size().abs() as i32, // Insert size for paired reads
            record.is_first_in_template(),     // Is this read1 of a pair?
            record.is_proper_pair(),           // Are mate pairs properly oriented?
            record.is_reverse(),               // Is read on reverse strand?
            &cigar_buf,            // CIGAR operations
            md_tag,                // MD tag (optional)
            seq_type,              // Sequencing type
            flags,                 // Which modules to calculate
        );
    }

    // -------------------------------------------------------------------------
    // Step 6: Handle empty contigs
    // -------------------------------------------------------------------------
    if !has_reads {
        return Ok(None);
    }

    // -------------------------------------------------------------------------
    // Step 7: Finalize strand-specific arrays for phagetermini
    // -------------------------------------------------------------------------
    // During processing, we track read starts/ends separately by strand.
    // Now we combine them according to sequencing type:
    // - Short reads: Only count starts from + strand, ends from - strand
    // - Long reads: Sum both strands
    if flags.phagetermini {
        arrays.finalize_strands(seq_type);
    }

    Ok(Some(arrays))
}
