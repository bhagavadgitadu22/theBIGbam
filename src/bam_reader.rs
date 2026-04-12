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

    let mut arrays = FeatureArrays::new(ref_length);
    let need_md = flags.needs_md();
    let mut has_reads = false;
    let mut primary_count: u64 = 0;
    let mut cigar_buf: Vec<(u32, u32)> = Vec::with_capacity(16);
    let mut seq_buf: Vec<u8> = Vec::with_capacity(256);
    let need_seq = flags.mapping_metrics;

    // Sum of time spent strictly inside process_read(). Subtracted from total
    // loop wall below to yield records_ns (decode + per-record filtering).
    let mut process_read_ns: u64 = 0;
    let t_loop = if time_on { Some(std::time::Instant::now()) } else { None };

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped() {
            continue;
        }

        has_reads = true;

        // Count primary alignments (not secondary, not supplementary)
        if !record.is_secondary() && !record.is_supplementary() {
            primary_count += 1;
        }

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

        // Extract SEQ bytes for sequence tracking (reuses buffer to avoid per-read allocation)
        if need_seq {
            seq_buf.clear();
            let seq = record.seq();
            let seq_len = seq.len();
            seq_buf.reserve(seq_len);
            for i in 0..seq_len {
                seq_buf.push(seq[i]);
            }
        }

        // Compute template length and proper_pair
        // BAM proper pair flag is unreliable (mapper-dependent), so recompute for all short-paired reads
        let (template_length, is_proper_pair, is_non_inward) = if seq_type.is_short_paired() {
            let pos1 = record.pos() as usize;
            let pos2 = record.mpos() as usize;

            let (corrected_tlen, corrected_proper) = if circular {
                // Circular: shortest path around the genome using actual (non-duplicated) length
                let p1 = pos1 % ref_length;
                let p2 = pos2 % ref_length;
                let direct = (p1 as i32 - p2 as i32).abs();
                let wrapped = ref_length as i32 - direct;
                let tlen = direct.min(wrapped);

                let same_ref = record.tid() == record.mtid();
                let opposite_strands = record.is_reverse() != record.is_mate_reverse();
                let reasonable = tlen > 0 && tlen < 10000;
                (tlen, same_ref && opposite_strands && reasonable)
            } else {
                // Linear: use absolute insert size, recompute proper pair
                let tlen = record.insert_size().abs() as i32;
                let same_ref = record.tid() == record.mtid();
                let opposite_strands = record.is_reverse() != record.is_mate_reverse();
                let reasonable = tlen > 0 && tlen < 10000;
                (tlen, same_ref && opposite_strands && reasonable)
            };

            let non_inward = record.tid() == record.mtid()
                && !record.is_mate_unmapped()
                && !corrected_proper;

            (corrected_tlen, corrected_proper, non_inward)
        } else {
            (record.insert_size().abs() as i32, record.is_proper_pair(), false)
        };

        // Track circularising reads (primary alignments only, circular mode only)
        // Origin-crossing reads: alignment end exceeds contig length
        let _is_circularising_read = if circular && !record.is_secondary() && !record.is_supplementary() {
            let raw_start = record.pos() as usize;
            let raw_end = cigar_view.end_pos() as usize;

            if raw_end > ref_length {
                arrays.circularising_reads_count += 1;
                // Gate: confirm circularity if ≥20bp mapped on both sides
                let left_overlap = ref_length - raw_start;
                let right_overlap = raw_end - ref_length;
                if left_overlap >= 20 && right_overlap >= 20 {
                    arrays.circularising_confirmed = true;
                }
                arrays.circularising_min_overlaps.push(left_overlap.min(right_overlap));
                true
            } else {
                false
            }
        } else {
            false
        };

        // Track circularising inserts & contig-end anomalies (primary paired-end only)
        // Skip reads already counted as circularising reads to avoid double-counting
        if seq_type.is_short_paired() && !record.is_secondary() && !record.is_supplementary() {
            // Detect reads near contig ends for circularisation analysis
            let pos = record.pos() as usize;
            let near_left = pos < 1000;
            let near_right = pos >= ref_length.saturating_sub(1000);

            if near_left || near_right {
                // Circularising inserts: read1 near one end, mate near opposite end
                if record.is_first_in_template() {
                    if circular {
                        // In SAM-spec circular mode, both mates have POS < LN.
                        // Junction-spanning pairs have one mate near each end.
                        let mpos = record.mpos() as usize;
                        let mate_near_left = mpos < 1000;
                        let mate_near_right = mpos >= ref_length.saturating_sub(1000);
                        if (near_left && mate_near_right) || (near_right && mate_near_left) {
                            let bam_proper = record.is_proper_pair();
                            let bam_non_inward = !record.is_proper_pair()
                                && record.tid() == record.mtid()
                                && !record.is_mate_unmapped();
                            if bam_proper || bam_non_inward {
                                arrays.circularising_inserts_count += 1;
                                arrays.circularising_insert_sizes.push(template_length);
                            }
                        }
                    } else if is_non_inward {
                        // Linear mode: non-inward pairs with mate on opposite end
                        let mpos = record.mpos() as usize;
                        let mate_near_left = mpos < 1000;
                        let mate_near_right = mpos >= ref_length.saturating_sub(1000);
                        if (near_left && mate_near_right) || (near_right && mate_near_left) {
                            arrays.circularising_inserts_count += 1;
                            // Use wrapped distance: these pairs span the contig boundary,
                            // so the true insert size is ref_length - raw_tlen
                            let wrapped_tlen = ref_length as i32 - template_length;
                            arrays.circularising_insert_sizes.push(wrapped_tlen.abs());
                        }
                    }
                }

                // Unmapped mate with orientation suggesting mate beyond contig end
                // Forward on right end → mate expected further right → beyond right end
                // Reverse on left end  → mate expected further left  → beyond left end
                if record.is_mate_unmapped() {
                    if (near_right && !record.is_reverse()) || (near_left && record.is_reverse()) {
                        arrays.contig_end_unmapped_mates += 1;
                    }
                }

                // Mate on another contig near contig ends
                if record.tid() != record.mtid() && !record.is_mate_unmapped() {
                    arrays.contig_end_mates_mapped_on_another_contig += 1;
                }
            }

            // Collect all proper-pair insert sizes for overall mean/median baseline
            if record.is_first_in_template() && template_length > 0 && is_proper_pair {
                arrays.all_proper_insert_sizes.push(template_length);
            }
        }

        let t_pr = if time_on { Some(std::time::Instant::now()) } else { None };
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
            if need_seq { &seq_buf } else { &[] },
            record.mapq(),
            seq_type,
            flags,
            circular,
            min_clipping_length,
        );
        if let Some(t) = t_pr {
            process_read_ns += t.elapsed().as_nanos() as u64;
        }
    }

    // Flush loop timings before any early return — records_ns = (total loop
    // wall) − (time inside process_read).
    if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_loop) {
        let total_ns = t.elapsed().as_nanos() as u64;
        pt.records_ns = total_ns.saturating_sub(process_read_ns);
        pt.process_read_ns = process_read_ns;
    }

    if !has_reads {
        return Ok(None);
    }

    // Check coverage percentage before post-processing and feature calculation
    // This saves time on low-coverage contigs by skipping:
    // 1. finalize_strands() (phagetermini-specific, expensive)
    // 2. All feature compression in the caller (very expensive)
    // 3. Database writes
    let t_fin = if time_on { Some(std::time::Instant::now()) } else { None };
    let coverage_pct = arrays.coverage_percentage();
    if coverage_pct < min_aligned_fraction {
        if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_fin) {
            pt.finalize_ns = t.elapsed().as_nanos() as u64;
        }
        return Ok(None);
    }

    // Finalize strand-specific reads_starts/reads_ends for phagetermini
    // This merges forward/reverse strand tracking into final arrays
    // Only needed for phagetermini; other features are already finalized
    if flags.phagetermini {
        arrays.finalize_strands(seq_type);
    }
    if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_fin) {
        pt.finalize_ns = t.elapsed().as_nanos() as u64;
    }

    Ok(Some((arrays, coverage_pct, primary_count)))
}

// ============================================================================
// Read Count Statistics
// ============================================================================

/// Get total read count from BAM index (fast, no iteration needed).
/// Returns total reads (mapped + unmapped) from index_stats().
pub fn get_total_read_count(bam_path: &Path) -> Result<u64> {
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed BAM: {}", bam_path.display()))?;

    let stats = bam.index_stats()
        .with_context(|| format!("Failed to get index stats from: {}", bam_path.display()))?;

    // Sum mapped + unmapped across all entries (including tid=-1 for unmapped)
    let total: u64 = stats.iter().map(|(_, _, mapped, unmapped)| mapped + unmapped).sum();

    Ok(total)
}
