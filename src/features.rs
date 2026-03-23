//! Feature calculation from aligned sequencing reads.
//!
//! Calculates coverage, phage termini features (reads_starts, reads_ends, tau),
//! and assembly check features (clippings, indels, mismatches, orientations).
//!
//! All features are calculated in a single pass through the BAM reads for efficiency.

use crate::cigar::{raw_cigar_consumes_ref, raw_cigar_is_clipping, raw_boundary_event_length};
use crate::circular::{increment_circular, increment_circular_long, increment_range};
use crate::types::{FeatureAnnotation, FeatureMap, SequencingType};
use std::collections::HashMap;

/// Maximum number of unique sequence variants tracked per position.
/// 10 is enough to reliably capture the dominant biological variant
/// while bounding memory usage.
const MAX_SEQS_PER_POS: usize = 10;

/// Track a sequence variant at a given position, capping the number of
/// unique variants stored per position to [`MAX_SEQS_PER_POS`].
fn track_sequence(
    map: &mut HashMap<usize, HashMap<Vec<u8>, u32>>,
    pos: usize,
    seq: &[u8],
) {
    let inner = map.entry(pos).or_default();
    if let Some(count) = inner.get_mut(seq) {
        *count += 1;
    } else if inner.len() < MAX_SEQS_PER_POS {
        inner.insert(seq.to_vec(), 1);
    }
}

// ============================================================================
// Feature Arrays - Central Data Structure
// ============================================================================

/// Container for all feature arrays calculated during processing.
///
/// Each array has one element per base pair in the reference sequence.
/// For a 50kb genome, each array has 50,000 elements.
///
/// # Memory Usage
/// With 18 arrays × 50,000 positions × 8 bytes = ~7MB per contig.
/// This is allocated once and reused, not per-read.
///
/// # Why separate arrays instead of one struct per position?
/// This "struct of arrays" layout is more cache-friendly than an
/// "array of structs" layout. When incrementing coverage for 150
/// consecutive positions, all the data we need is contiguous in memory.
pub struct FeatureArrays {
    // -------------------------------------------------------------------------
    // Coverage Module
    // -------------------------------------------------------------------------
    /// Number of primary reads overlapping each position (standard read depth).
    /// Only counts reads that are neither secondary nor supplementary.
    pub primary_reads: Vec<u64>,

    /// Number of primary reads on the plus strand at each position.
    pub primary_reads_plus_only: Vec<u64>,

    /// Number of primary reads on the minus strand at each position.
    pub primary_reads_minus_only: Vec<u64>,

    /// Number of secondary alignments at each position.
    /// Secondary reads have flag 0x100 set.
    pub secondary_reads: Vec<u64>,

    /// Number of supplementary alignments at each position.
    /// Supplementary reads have flag 0x800 set.
    pub supplementary_reads: Vec<u64>,

    /// Sum of MAPQ values at each position (for primary reads only).
    /// Average MAPQ = sum_mapq / primary_reads.
    pub sum_mapq: Vec<u64>,

    // -------------------------------------------------------------------------
    // Phagetermini Module
    // -------------------------------------------------------------------------
    /// Coverage counting only "clean" reads (no clipping/mismatch at ends).
    /// Used for phage terminus detection - noisy reads would obscure the signal.
    pub coverage_reduced: Vec<u64>,

    /// Count of read 5' ends at each position.
    /// High peaks indicate potential DNA packaging/cutting sites.
    pub reads_starts: Vec<u64>,

    /// Count of read 3' ends at each position.
    /// Combined with reads_starts to identify terminus types.
    pub reads_ends: Vec<u64>,

    /// Event lengths at each read start position (0 for exact match, clip/insertion length for near-match).
    /// Used for per-position Mean/Median/Std statistics on reads_starts.
    pub start_event_lengths: Vec<Vec<u32>>,

    /// Event lengths at each read end position (0 for exact match, clip/insertion length for near-match).
    /// Used for per-position Mean/Median/Std statistics on reads_ends.
    pub end_event_lengths: Vec<Vec<u32>>,

    // -------------------------------------------------------------------------
    // Assemblycheck Module
    // -------------------------------------------------------------------------
    /// Lengths of soft/hard clippings at their 5' end (one entry per clipped read).
    /// High values may indicate assembly breakpoints or structural variants.
    pub left_clipping_lengths: Vec<Vec<u32>>,

    /// Lengths of soft/hard clippings at their 3' end (one entry per clipped read).
    pub right_clipping_lengths: Vec<Vec<u32>>,

    /// Lengths of insertion events from CIGAR (one entry per insertion).
    pub insertion_lengths: Vec<Vec<u32>>,

    /// Deletion events from CIGAR (sequence in reference but not in read).
    pub deletions: Vec<u64>,

    /// Lengths of deletion events from CIGAR (one entry per deletion at each position).
    pub deletion_lengths: Vec<Vec<u32>>,

    /// Base mismatches (from MD tag).
    /// High mismatch rates may indicate errors or strain variation.
    pub mismatches: Vec<u64>,

    /// Per-position counts of mismatched read bases: [A, C, G, T] counts.
    /// Dense array for O(1) access; 16 bytes per position.
    pub mismatch_base_counts: Vec<[u32; 4]>,

    /// Sparse map: position → (insertion_sequence → count).
    /// Full insertion sequences (no truncation). Only positions with events are allocated.
    pub insertion_sequences: HashMap<usize, HashMap<Vec<u8>, u32>>,

    /// Sparse map: position → (clip_sequence → count).
    /// Clip sequences truncated to 20bp. Only positions with events are allocated.
    pub left_clip_sequences: HashMap<usize, HashMap<Vec<u8>, u32>>,

    /// Sparse map: position → (clip_sequence → count).
    /// Clip sequences truncated to 20bp. Only positions with events are allocated.
    pub right_clip_sequences: HashMap<usize, HashMap<Vec<u8>, u32>>,

    /// Sparse map: position → (clip_sequence → count).
    /// Only clips < min_clipping_length (used by reads_starts). Truncated to 20bp.
    pub start_clip_sequences: HashMap<usize, HashMap<Vec<u8>, u32>>,

    /// Sparse map: position → (clip_sequence → count).
    /// Only clips < min_clipping_length (used by reads_ends). Truncated to 20bp.
    pub end_clip_sequences: HashMap<usize, HashMap<Vec<u8>, u32>>,

    /// Sum of read lengths at each position (for computing average).
    /// Only used for long-read data.
    pub sum_read_lengths: Vec<u64>,

    /// Count of reads at each position (denominator for read length average).
    pub count_read_lengths: Vec<u64>,

    /// Sum of insert sizes at each position (for computing average).
    /// Only used for paired-end data. Insert size = distance between mates.
    pub sum_insert_sizes: Vec<u64>,

    /// Count of proper pairs at each position (denominator for insert size average).
    pub count_insert_sizes: Vec<u64>,

    /// Reads where mates are on the same contig but not in proper inward orientation.
    /// High values may indicate structural rearrangements.
    pub non_inward_pairs: Vec<u64>,

    /// Reads where the mate is not mapped at all.
    pub mate_not_mapped: Vec<u64>,

    /// Reads where the mate is mapped to a different contig.
    pub mate_on_another_contig: Vec<u64>,

    // -------------------------------------------------------------------------
    // Circularising reads tracking
    // -------------------------------------------------------------------------
    /// Total count of reads that support genome circularity.
    /// In circular mode: reads whose alignment extends past LN (SAM-spec circular).
    pub circularising_reads_count: u64,
    /// At least one read with ≥20bp mapped on both sides of the junction.
    pub circularising_confirmed: bool,
    /// For each circularising read, min(left_overlap, right_overlap).
    pub circularising_min_overlaps: Vec<usize>,

    /// Count of non-inward read pairs spanning both contig ends (paired-end only).
    pub circularising_inserts_count: u64,
    /// Template lengths of circularising insert pairs (for mean/median).
    pub circularising_insert_sizes: Vec<i32>,
    /// Template lengths of all proper pairs (read1 only, for overall mean/median baseline).
    pub all_proper_insert_sizes: Vec<i32>,
    /// Reads near contig ends with unmapped mate and outward-facing orientation.
    pub contig_end_unmapped_mates: u64,
    /// Reads near contig ends with mate on a different contig.
    pub contig_end_mates_mapped_on_another_contig: u64,

    /// Number of primary reads with at least one clean terminus (no clip/mismatch).
    /// For long reads: each read is split in two, so a read with both termini clean
    /// counts as 2, and a read with only one clean terminus counts as 1.
    /// For short reads: counts reads where the start is clean.
    pub clean_reads_count: u64,

    // -------------------------------------------------------------------------
    // Internal: Strand-specific tracking for phagetermini
    // -------------------------------------------------------------------------
    // We track starts/ends separately by strand during processing,
    // then combine them at the end based on sequencing type.

    /// Read starts on the forward (+) strand
    start_plus: Vec<u64>,
    /// Read starts on the reverse (-) strand
    start_minus: Vec<u64>,
    /// Read ends on the forward (+) strand
    end_plus: Vec<u64>,
    /// Read ends on the reverse (-) strand
    end_minus: Vec<u64>,
}

impl FeatureArrays {
    /// Create new feature arrays for a given reference length.
    ///
    /// # Arguments
    /// * `ref_length` - Length of the reference sequence in base pairs
    ///
    /// # Rust Concept - `vec![value; count]` macro:
    /// Creates a Vec with `count` copies of `value`.
    /// `vec![0u64; 50000]` creates a Vec with 50,000 zeros.
    #[inline]
    pub fn new(ref_length: usize) -> Self {
        Self {
            // Initialize all arrays to zero
            primary_reads: vec![0u64; ref_length],
            primary_reads_plus_only: vec![0u64; ref_length],
            primary_reads_minus_only: vec![0u64; ref_length],
            secondary_reads: vec![0u64; ref_length],
            supplementary_reads: vec![0u64; ref_length],
            sum_mapq: vec![0u64; ref_length],
            coverage_reduced: vec![0u64; ref_length],
            reads_starts: vec![0u64; ref_length],
            reads_ends: vec![0u64; ref_length],
            start_event_lengths: vec![Vec::new(); ref_length],
            end_event_lengths: vec![Vec::new(); ref_length],
            left_clipping_lengths: vec![Vec::new(); ref_length],
            right_clipping_lengths: vec![Vec::new(); ref_length],
            insertion_lengths: vec![Vec::new(); ref_length],
            deletions: vec![0u64; ref_length],
            deletion_lengths: vec![Vec::new(); ref_length],
            mismatches: vec![0u64; ref_length],
            mismatch_base_counts: vec![[0u32; 4]; ref_length],
            insertion_sequences: HashMap::new(),
            left_clip_sequences: HashMap::new(),
            right_clip_sequences: HashMap::new(),
            start_clip_sequences: HashMap::new(),
            end_clip_sequences: HashMap::new(),
            sum_read_lengths: vec![0u64; ref_length],
            count_read_lengths: vec![0u64; ref_length],
            sum_insert_sizes: vec![0u64; ref_length],
            count_insert_sizes: vec![0u64; ref_length],
            non_inward_pairs: vec![0u64; ref_length],
            mate_not_mapped: vec![0u64; ref_length],
            mate_on_another_contig: vec![0u64; ref_length],
            circularising_reads_count: 0,
            circularising_confirmed: false,
            circularising_min_overlaps: Vec::new(),
            circularising_inserts_count: 0,
            circularising_insert_sizes: Vec::new(),
            all_proper_insert_sizes: Vec::new(),
            contig_end_unmapped_mates: 0,
            contig_end_mates_mapped_on_another_contig: 0,
            clean_reads_count: 0,
            start_plus: vec![0u64; ref_length],
            start_minus: vec![0u64; ref_length],
            end_plus: vec![0u64; ref_length],
            end_minus: vec![0u64; ref_length],
        }
    }

    /// Get the reference length (number of positions in arrays).
    #[inline]
    pub fn ref_length(&self) -> usize {
        self.primary_reads.len()
    }

    /// Calculate what percentage of the reference has at least 1x coverage.
    ///
    /// # Returns
    /// Percentage as a float (0.0 to 100.0)
    ///
    /// # Rust Concept - Iterator chain:
    /// `.iter().filter(...).count()` is like Python's:
    /// `sum(1 for x in self.primary_reads if x > 0)`
    #[inline]
    pub fn coverage_percentage(&self) -> f64 {
        // Count how many positions have at least 1 primary read
        let covered_bp = self.primary_reads.iter().filter(|&&x| x > 0).count();
        // Convert to percentage
        (covered_bp as f64 / self.ref_length() as f64) * 100.0
    }

    /// Finalize phagetermini strand arrays based on sequencing type.
    ///
    /// During processing, we track read starts/ends separately by strand.
    /// At the end, we combine them differently based on the sequencing type:
    ///
    /// - **Short reads**: For paired-end short reads, the 5' end of read1 (+ strand)
    ///   represents the true fragment start, and the 5' end of read2 (- strand)
    ///   represents the true fragment end.
    ///
    /// - **Long reads**: Both strands contribute equally, so we sum them.
    ///
    /// # Rust Concept - `std::mem::swap`:
    /// Efficiently exchanges the contents of two variables without copying.
    /// This is O(1) - just swaps pointers, doesn't copy 50,000 elements!
    pub fn finalize_strands(&mut self, seq_type: SequencingType) {
        match seq_type {
            SequencingType::ShortPaired | SequencingType::ShortSingle => {
                // Short reads: starts from + strand, ends from - strand
                // swap() is O(1) - just exchanges the Vec pointers
                std::mem::swap(&mut self.reads_starts, &mut self.start_plus);
                std::mem::swap(&mut self.reads_ends, &mut self.end_minus);
            }
            SequencingType::Long => {
                for i in 0..self.ref_length() {
                    self.reads_starts[i] = self.start_plus[i] + self.start_minus[i];
                    self.reads_ends[i] = self.end_plus[i] + self.end_minus[i];
                }
            }
        }
    }

    /// Convert to FeatureMap for phagetermini features.
    ///
    /// # Rust Concept - `.clone()`:
    /// Creates a deep copy of the Vec. This is necessary because we're
    /// moving data out of the struct into a HashMap.
    pub fn to_phagetermini_map(&self) -> FeatureMap {
        let mut results = FeatureMap::with_capacity(3);
        results.insert("coverage_reduced".to_string(), self.coverage_reduced.clone());
        results.insert("reads_starts".to_string(), self.reads_starts.clone());
        results.insert("reads_ends".to_string(), self.reads_ends.clone());
        results
    }

    /// Compute count of clipping/insertion events at each position.
    fn compute_counts(lengths: &[Vec<u32>]) -> Vec<u64> {
        lengths.iter().map(|v| v.len() as u64).collect()
    }

    /// Compute mean length at each position (0 if no events).
    fn compute_means(lengths: &[Vec<u32>]) -> Vec<u64> {
        lengths.iter().map(|v| {
            if v.is_empty() {
                0
            } else {
                let sum: u64 = v.iter().map(|&x| x as u64).sum();
                sum / v.len() as u64
            }
        }).collect()
    }

    /// Compute median length at each position (0 if no events).
    fn compute_medians(lengths: &[Vec<u32>]) -> Vec<u64> {
        lengths.iter().map(|v| {
            if v.is_empty() {
                0
            } else {
                let mut sorted = v.clone();
                sorted.sort_unstable();
                sorted[sorted.len() / 2] as u64
            }
        }).collect()
    }

    /// Compute standard deviation at each position (0 if no events).
    fn compute_std_devs(lengths: &[Vec<u32>]) -> Vec<u64> {
        lengths.iter().map(|v| {
            if v.len() <= 1 {
                0
            } else {
                let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
                let variance = v.iter()
                    .map(|&x| {
                        let diff = x as f64 - mean;
                        diff * diff
                    })
                    .sum::<f64>() / v.len() as f64;
                variance.sqrt() as u64
            }
        }).collect()
    }

    /// Compute mean coverage across all positions.
    /// Returns the average of primary_reads as f64.
    pub fn coverage_mean(&self) -> f64 {
        if self.primary_reads.is_empty() {
            0.0
        } else {
            self.primary_reads.iter().sum::<u64>() as f64 / self.primary_reads.len() as f64
        }
    }

    /// Compute median coverage across all positions.
    /// Returns the middle value of sorted primary_reads as f64.
    pub fn coverage_median(&self) -> f64 {
        if self.primary_reads.is_empty() {
            0.0
        } else {
            let mut sorted: Vec<u64> = self.primary_reads.clone();
            sorted.sort_unstable();
            sorted[sorted.len() / 2] as f64
        }
    }

    /// Compute trimmed mean coverage (trim `trim_fraction` from both sides).
    /// For example, trim_fraction=0.05 trims the bottom and top 5%.
    pub fn coverage_trimmed_mean(&self, trim_fraction: f64) -> f64 {
        if self.primary_reads.is_empty() {
            return 0.0;
        }
        let mut sorted: Vec<u64> = self.primary_reads.clone();
        sorted.sort_unstable();
        let n = sorted.len();
        let trim = (n as f64 * trim_fraction) as usize;
        let trimmed = &sorted[trim..n.saturating_sub(trim)];
        if trimmed.is_empty() {
            0.0
        } else {
            trimmed.iter().sum::<u64>() as f64 / trimmed.len() as f64
        }
    }

    /// Compute dominant mismatch base per position.
    /// Returns Vec of (base_char, percentage_x10) for each position.
    /// base_char: b'A', b'C', b'G', b'T' or 0 if no mismatches.
    /// percentage_x10: prevalence × 1000 relative to primary_reads (e.g., 35 = 3.5%)
    /// Only includes positions where prevalence >= threshold.
    pub fn compute_dominant_mismatch_bases(
        &self,
        primary_reads: &[u64],
        threshold: f64,
    ) -> Vec<(u8, i32)> {
        self.mismatch_base_counts.iter().enumerate().map(|(pos, counts)| {
            let total: u32 = counts.iter().sum();
            if total == 0 {
                return (0, 0);
            }
            let (max_idx, &max_count) = counts.iter().enumerate()
                .max_by_key(|(_, &c)| c).unwrap();
            let total_reads = if pos < primary_reads.len() { primary_reads[pos] } else { 0 };
            if total_reads == 0 {
                return (0, 0);
            }
            let prevalence = max_count as f64 / total_reads as f64;
            if prevalence < threshold {
                return (0, 0);
            }
            let base = match max_idx {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => 0,
            };
            let percentage_x10 = (prevalence * 1000.0) as i32;
            (base, percentage_x10)
        }).collect()
    }

    /// Compute dominant sequence per position from a sparse HashMap.
    /// Returns HashMap<usize, (String, i32)> mapping position -> (sequence, percentage_x10)
    /// Prevalence is relative to primary_reads (total reads at that position).
    /// Only includes positions where prevalence >= threshold.
    pub fn compute_dominant_sequences(
        seqs: &HashMap<usize, HashMap<Vec<u8>, u32>>,
        primary_reads: &[u64],
        threshold: f64,
    ) -> HashMap<usize, (String, i32)> {
        seqs.iter().filter_map(|(&pos, seq_counts)| {
            let (dom_seq, &max_count) = seq_counts.iter()
                .max_by_key(|(_, &c)| c).unwrap();
            let total_reads = if pos < primary_reads.len() { primary_reads[pos] } else { 0 };
            if total_reads == 0 { return None; }
            let prevalence = max_count as f64 / total_reads as f64;
            if prevalence < threshold { return None; }
            let percentage_x10 = (prevalence * 1000.0) as i32;
            let seq_str = String::from_utf8_lossy(dom_seq).into_owned();
            Some((pos, (seq_str, percentage_x10)))
        }).collect()
    }

    /// Convert to FeatureMap for assemblycheck features.
    ///
    /// Only includes features relevant to the sequencing type:
    /// - Long reads: includes read_lengths
    /// - Paired reads: includes insert_sizes, mate orientation data
    /// 
    /// Returns a tuple: (FeatureMap with counts, statistics for clippings/insertions)
    pub fn to_assemblycheck_data(&self, seq_type: SequencingType) -> (
        FeatureMap,
        HashMap<String, (Vec<u64>, Vec<u64>, Vec<u64>)> // mean, median, std
    ) {
        let mut results = FeatureMap::with_capacity(10);
        let mut stats = HashMap::new();
        
        // Clipping and insertion with statistics
        results.insert("left_clippings".to_string(), Self::compute_counts(&self.left_clipping_lengths));
        stats.insert("left_clippings".to_string(), (
            Self::compute_means(&self.left_clipping_lengths),
            Self::compute_medians(&self.left_clipping_lengths),
            Self::compute_std_devs(&self.left_clipping_lengths),
        ));
        
        results.insert("right_clippings".to_string(), Self::compute_counts(&self.right_clipping_lengths));
        stats.insert("right_clippings".to_string(), (
            Self::compute_means(&self.right_clipping_lengths),
            Self::compute_medians(&self.right_clipping_lengths),
            Self::compute_std_devs(&self.right_clipping_lengths),
        ));
        
        results.insert("insertions".to_string(), Self::compute_counts(&self.insertion_lengths));
        stats.insert("insertions".to_string(), (
            Self::compute_means(&self.insertion_lengths),
            Self::compute_medians(&self.insertion_lengths),
            Self::compute_std_devs(&self.insertion_lengths),
        ));
        
        results.insert("deletions".to_string(), self.deletions.clone());
        results.insert("mismatches".to_string(), self.mismatches.clone());

        // Long reads: include read length data
        if seq_type.is_long() {
            results.insert("sum_read_lengths".to_string(), self.sum_read_lengths.clone());
            results.insert("count_read_lengths".to_string(), self.count_read_lengths.clone());
        }

        // Paired reads: include insert size and orientation data
        if seq_type.is_short_paired() {
            results.insert("sum_insert_sizes".to_string(), self.sum_insert_sizes.clone());
            results.insert("count_insert_sizes".to_string(), self.count_insert_sizes.clone());
            results.insert("non_inward_pairs".to_string(), self.non_inward_pairs.clone());
            results.insert("mate_not_mapped".to_string(), self.mate_not_mapped.clone());
            results.insert("mate_on_another_contig".to_string(), self.mate_on_another_contig.clone());
        }

        (results, stats)
    }
}

// ============================================================================
// CDS Index for Codon Tracking
// ============================================================================

/// A single CDS interval for binary search lookup.
#[derive(Clone, Debug)]
pub struct CdsInterval {
    /// 1-based start position on contig
    pub start: i64,
    /// 1-based end position on contig (inclusive)
    pub end: i64,
    /// Strand: 1 or -1
    pub strand: i64,
    /// Nucleotide sequence of the CDS (already strand-corrected)
    pub nucleotide_sequence: String,
}

/// Sorted CDS index for fast position-to-CDS lookup via binary search.
pub struct CdsIndex {
    /// CDS intervals sorted by start position
    intervals: Vec<CdsInterval>,
}

impl CdsIndex {
    /// Build a CDS index from annotations for a specific contig.
    pub fn from_annotations(annotations: &[FeatureAnnotation], contig_id: i64) -> Self {
        let mut intervals: Vec<CdsInterval> = annotations
            .iter()
            .filter(|a| {
                a.contig_id == contig_id
                    && a.feature_type == "CDS"
                    && a.nucleotide_sequence.is_some()
            })
            .map(|a| CdsInterval {
                start: a.start,
                end: a.end,
                strand: a.strand,
                nucleotide_sequence: a.nucleotide_sequence.clone().unwrap(),
            })
            .collect();
        intervals.sort_by_key(|c| c.start);
        Self { intervals }
    }

    /// Find the CDS covering a given 1-based position, if any.
    /// Uses binary search on start positions, then scans backwards for overlapping CDS.
    pub fn find_cds(&self, pos: i64) -> Option<&CdsInterval> {
        if self.intervals.is_empty() {
            return None;
        }
        // Find rightmost CDS with start <= pos
        let idx = self.intervals.partition_point(|c| c.start <= pos);
        if idx == 0 {
            return None;
        }
        // Scan backwards to find a CDS that covers pos
        for i in (0..idx).rev() {
            if self.intervals[i].end >= pos {
                return Some(&self.intervals[i]);
            }
        }
        None
    }

    /// Check if this index has any CDS intervals.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }
}

/// Codon change info: (codon_category, codon_change, aa_change)
type CodonInfo = (String, String, String);

/// Standard genetic code lookup for codon translation.
pub(crate) fn translate_codon(codon: &[u8; 3]) -> Option<(char, &'static str)> {
    let c = [codon[0].to_ascii_uppercase(), codon[1].to_ascii_uppercase(), codon[2].to_ascii_uppercase()];
    Some(match &c {
        b"TTT" | b"TTC" => ('F', "Phenylalanine"),
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => ('L', "Leucine"),
        b"ATT" | b"ATC" | b"ATA" => ('I', "Isoleucine"),
        b"ATG" => ('M', "Methionine"),
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => ('V', "Valine"),
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => ('S', "Serine"),
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => ('P', "Proline"),
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => ('T', "Threonine"),
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => ('A', "Alanine"),
        b"TAT" | b"TAC" => ('Y', "Tyrosine"),
        b"TAA" | b"TAG" | b"TGA" => ('*', "Stop"),
        b"CAT" | b"CAC" => ('H', "Histidine"),
        b"CAA" | b"CAG" => ('Q', "Glutamine"),
        b"AAT" | b"AAC" => ('N', "Asparagine"),
        b"AAA" | b"AAG" => ('K', "Lysine"),
        b"GAT" | b"GAC" => ('D', "Aspartate"),
        b"GAA" | b"GAG" => ('E', "Glutamate"),
        b"TGT" | b"TGC" => ('C', "Cysteine"),
        b"TGG" => ('W', "Tryptophan"),
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => ('R', "Arginine"),
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => ('G', "Glycine"),
        _ => return None,
    })
}

/// Compute codon changes from per-position mismatch summaries (post-processing).
///
/// Instead of analyzing codons per-read (O(reads_with_mismatches)), this computes
/// codon changes once from the accumulated mismatch_base_counts (O(positions_with_mismatches)).
/// For each position with mismatches above the threshold, finds the dominant mismatch base,
/// looks up the covering CDS, and classifies the codon change.
///
/// # Trade-off
/// Each mismatch position is evaluated independently. In the rare case (~0.01%) where two
/// mismatches in the same codon originate from the same read, they are not combined into
/// a single mutant codon. This has negligible impact on classification accuracy.
///
/// # Arguments
/// * `mismatch_base_counts` - Per-position [A, C, G, T] mismatch counts (dense array)
/// * `primary_reads` - Per-position primary read depth (for prevalence threshold)
/// * `cds_index` - CDS lookup index for this contig (None if no annotations)
/// * `ref_length` - Reference length for circular modulo
/// * `threshold` - Minimum prevalence fraction to report a dominant mismatch
pub fn compute_codon_changes_from_summaries(
    mismatch_base_counts: &[[u32; 4]],
    primary_reads: &[u64],
    cds_index: Option<&CdsIndex>,
    ref_length: usize,
    threshold: f64,
) -> HashMap<usize, CodonInfo> {
    let cds_idx = match cds_index {
        Some(idx) if !idx.is_empty() => idx,
        _ => return HashMap::new(),
    };

    let complement = |b: u8| -> u8 {
        match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        }
    };

    let bases = [b'A', b'C', b'G', b'T'];
    let mut result = HashMap::new();

    for (pos, counts) in mismatch_base_counts.iter().enumerate() {
        let total: u32 = counts.iter().sum();
        if total == 0 {
            continue;
        }

        // Find dominant mismatch base
        let (max_idx, &max_count) = counts.iter().enumerate()
            .max_by_key(|(_, &c)| c).unwrap();
        let total_reads = if pos < primary_reads.len() { primary_reads[pos] } else { 0 };
        if total_reads == 0 {
            continue;
        }
        let prevalence = max_count as f64 / total_reads as f64;
        if prevalence < threshold {
            continue;
        }

        let dominant_base = bases[max_idx];
        let pos_1based = (pos % ref_length) as i64 + 1;

        if let Some(cds) = cds_idx.find_cds(pos_1based) {
            let nuc_bytes = cds.nucleotide_sequence.as_bytes();

            let offset = if cds.strand == -1 {
                (cds.end - pos_1based) as usize
            } else {
                (pos_1based - cds.start) as usize
            };

            if offset >= nuc_bytes.len() {
                continue;
            }

            let codon_idx = offset / 3;
            let pos_in_codon = offset % 3;

            let codon_start = codon_idx * 3;
            let codon_end = codon_start + 3;
            if codon_end > nuc_bytes.len() {
                continue;
            }

            let ref_codon = &nuc_bytes[codon_start..codon_end];

            // Build mutant codon with dominant mismatch at this position
            let mut mut_codon = [
                ref_codon[0].to_ascii_uppercase(),
                ref_codon[1].to_ascii_uppercase(),
                ref_codon[2].to_ascii_uppercase(),
            ];
            let alt_base = if cds.strand == -1 {
                complement(dominant_base)
            } else {
                dominant_base
            };
            mut_codon[pos_in_codon] = alt_base;

            // Translate both codons
            let ref_codon_upper = [
                ref_codon[0].to_ascii_uppercase(),
                ref_codon[1].to_ascii_uppercase(),
                ref_codon[2].to_ascii_uppercase(),
            ];
            let Some((ref_aa, _)) = translate_codon(&ref_codon_upper) else { continue };
            let Some((mut_aa, mut_name)) = translate_codon(&mut_codon) else { continue };

            let category = if ref_aa == mut_aa { "Synonymous" } else { "Non-synonymous" };
            let codon_change_str = String::from_utf8_lossy(&mut_codon).into_owned();
            let aa_change_str = format!("{} ({})", mut_aa, mut_name);

            result.insert(pos, (category.to_string(), codon_change_str, aa_change_str));
        } else {
            result.insert(pos, ("Intergenic".to_string(), String::new(), String::new()));
        }
    }

    result
}

// ============================================================================
// Module Flags - Configuration
// ============================================================================

/// Configuration flags indicating which feature modules to calculate.
///
/// This allows skipping unnecessary calculations. For example, if the user
/// only wants coverage, we don't need to process CIGAR strings for indels.
///
/// # Rust Concept - `#[derive(Clone, Copy)]`:
/// `Copy` means this struct can be copied implicitly (like integers).
/// It's only 3 booleans (3 bytes), so copying is cheap.
#[derive(Clone, Copy)]
pub struct ModuleFlags {
    /// Calculate basic coverage (read depth)
    pub coverage: bool,
    /// Calculate misalignment metrics per position (clippings, indels, mismatches)
    pub mapping_metrics: bool,
    /// Calculate long-read specific metrics (read_lengths)
    pub long_read_metrics: bool,
    /// Calculate paired-read specific metrics (insert_sizes, non_inward_pairs, etc.)
    pub paired_read_metrics: bool,
    /// Calculate phage termini features (coverage_reduced, reads_starts, reads_ends, tau)
    pub phagetermini: bool,
}

impl ModuleFlags {
    /// Create flags from a list of module names.
    ///
    /// # Arguments
    /// * `modules` - List of module names matching types.rs Variable.module values:
    ///   "Coverage", "Misalignment", "Long-reads",
    ///   "Paired-reads", "Phage termini"
    ///
    /// # Note
    /// Coverage is automatically enabled when dependent modules are requested
    /// (mapping_metrics, long_read_metrics, paired_read_metrics need coverage for normalization).
    pub fn from_modules(modules: &[String]) -> Self {
        let has = |name: &str| modules.iter().any(|m| m.eq_ignore_ascii_case(name));

        let mapping_metrics = has("Misalignment");
        let long_read_metrics = has("Long-reads");
        let paired_read_metrics = has("Paired-reads");

        Self {
            // Coverage needed if explicitly requested or if dependent modules are enabled
            coverage: has("Coverage") || mapping_metrics || long_read_metrics || paired_read_metrics,
            mapping_metrics,
            long_read_metrics,
            paired_read_metrics,
            phagetermini: has("Phage termini"),
        }
    }

    /// Check if we need to extract the MD tag from reads.
    ///
    /// MD tag extraction has some overhead, so we only do it when needed.
    /// It's required for:
    /// - phagetermini: to check for mismatches at read ends
    /// - mapping_metrics: to count mismatches
    #[inline]
    pub fn needs_md(&self) -> bool {
        self.phagetermini || self.mapping_metrics
    }
}

// ============================================================================
// Core Processing Function
// ============================================================================

/// Process a single read and update all feature arrays in one pass.
///
/// This is the **heart of the feature calculation**. It's called once for each
/// read in the BAM file and updates all relevant arrays based on the read's
/// properties (position, CIGAR, MD tag, flags).
///
/// # Arguments
///
/// * `arrays` - Mutable reference to the feature arrays to update
/// * `ref_start` - Start position on reference (0-based)
/// * `ref_end` - End position on reference (0-based, exclusive)
/// * `query_length` - Length of the read sequence
/// * `template_length` - Insert size for paired reads (distance between mates)
/// * `is_read1` - True if this is read1 of a pair
/// * `is_proper_pair` - True if mates are properly oriented
/// * `is_reverse` - True if read is on reverse strand
/// * `is_secondary` - True if this is a secondary alignment
/// * `is_supplementary` - True if this is a supplementary alignment
/// * `mate_unmapped` - True if the mate is not mapped
/// * `mate_other_contig` - True if the mate is mapped to a different contig
/// * `cigar_raw` - CIGAR operations as (operation_char, length) tuples
/// * `md_tag` - Optional MD tag bytes for mismatch detection
/// * `mapq` - Mapping quality (0-255)
/// * `seq_type` - Sequencing type (affects which features to calculate)
/// * `flags` - Which modules are enabled
/// * `circular` - True if genome is circular (affects position calculations)
///
/// # Performance Optimizations
///
/// 1. **Branch hoisting**: Conditions like `if flags.coverage` are checked once,
///    not inside tight loops.
///
/// 2. **Non-wrapping fast path**: Most reads don't wrap around circular genomes,
///    so we have optimized code paths that skip modulo operations.
///
/// 3. **Combined loops**: When updating multiple arrays over the same range,
///    we combine them into one loop to improve cache locality.
///
/// 4. **Zero-allocation**: We work directly with raw CIGAR slices instead of
///    allocating new data structures.
#[inline]
pub fn process_read(
    arrays: &mut FeatureArrays,
    ref_start: i64,
    ref_end: i64,
    query_length: i32,
    template_length: i32,
    is_read1: bool,
    is_proper_pair: bool,
    is_reverse: bool,
    is_secondary: bool,
    is_supplementary: bool,
    mate_unmapped: bool,
    mate_other_contig: bool,
    cigar_raw: &[(u32, u32)],
    md_tag: Option<&[u8]>,
    seq: &[u8],
    mapq: u8,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_clipping_length: u32,
) {
    // -------------------------------------------------------------------------
    // Calculate positions
    // -------------------------------------------------------------------------
    let ref_length = arrays.ref_length();
    let raw_start = ref_start as usize;
    let raw_end = ref_end as usize;

    // For circular genomes (SAM-spec), apply modulo for array bounds.
    // Keep raw_end for increment_circular to detect wrapping correctly.
    let start = raw_start % ref_length;
    let end = raw_end;

    // -------------------------------------------------------------------------
    // Coverage module: primary mappings only
    // -------------------------------------------------------------------------
    if flags.coverage {
        // Track secondary and supplementary reads separately
        if is_secondary {
            if circular {
                if seq_type.is_long() {
                    increment_circular_long(&mut arrays.secondary_reads, raw_start, raw_end, 1);
                } else {
                    increment_circular(&mut arrays.secondary_reads, start, end, 1);
                }
            } else {
                increment_range(&mut arrays.secondary_reads, start, end, 1);
            }
        } else if is_supplementary {
            if circular {
                if seq_type.is_long() {
                    increment_circular_long(&mut arrays.supplementary_reads, raw_start, raw_end, 1);
                } else {
                    increment_circular(&mut arrays.supplementary_reads, start, end, 1);
                }
            } else {
                increment_range(&mut arrays.supplementary_reads, start, end, 1);
            }
        } else {
            // Only count primary mappings in main coverage
            let mapq_val = mapq as u64;
            if circular {
                if seq_type.is_long() {
                    increment_circular_long(&mut arrays.primary_reads, raw_start, raw_end, 1);
                    increment_circular_long(&mut arrays.sum_mapq, raw_start, raw_end, mapq_val);
                    // Also track strand-specific coverage
                    if is_reverse {
                        increment_circular_long(&mut arrays.primary_reads_minus_only, raw_start, raw_end, 1);
                    } else {
                        increment_circular_long(&mut arrays.primary_reads_plus_only, raw_start, raw_end, 1);
                    }
                } else {
                    increment_circular(&mut arrays.primary_reads, start, end, 1);
                    increment_circular(&mut arrays.sum_mapq, start, end, mapq_val);
                    // Also track strand-specific coverage
                    if is_reverse {
                        increment_circular(&mut arrays.primary_reads_minus_only, start, end, 1);
                    } else {
                        increment_circular(&mut arrays.primary_reads_plus_only, start, end, 1);
                    }
                }
            } else {
                increment_range(&mut arrays.primary_reads, start, end, 1);
                increment_range(&mut arrays.sum_mapq, start, end, mapq_val);
                // Also track strand-specific coverage
                if is_reverse {
                    increment_range(&mut arrays.primary_reads_minus_only, start, end, 1);
                } else {
                    increment_range(&mut arrays.primary_reads_plus_only, start, end, 1);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Long-reads module - primary mappings only
    // -------------------------------------------------------------------------
    if flags.long_read_metrics && seq_type.is_long() && !is_secondary && !is_supplementary {
        // --- Read lengths ---
        // Track sum and count to compute average read length at each position
        let ql = query_length as u64;

        // OPTIMIZATION: Non-wrapping case (most common)
        // Avoids modulo operation in tight loop
        if raw_end <= ref_length {
            for p in raw_start..raw_end {
                arrays.sum_read_lengths[p] += ql;
                arrays.count_read_lengths[p] += 1;
            }
        } else {
            // Wrapping case (read spans the origin of circular genome)
            for pos in raw_start..raw_end {
                let p = pos % ref_length;
                arrays.sum_read_lengths[p] += ql;
                arrays.count_read_lengths[p] += 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Paired-reads module - primary mappings only
    // -------------------------------------------------------------------------
    if flags.paired_read_metrics && seq_type.is_short_paired() && !is_secondary && !is_supplementary {
        // --- Insert sizes and mate orientation issues ---
        // Only track insert size for read1 with valid template length
        // (read2 would double-count the same fragment)
        let track_insert = is_read1 && template_length > 0;
        let tl = template_length as u64;

        // Classify mate issues
        let non_inward = !is_proper_pair && !mate_unmapped && !mate_other_contig;
        let mate_unmapped_flag = mate_unmapped;
        let mate_other_contig_flag = mate_other_contig;

        // OPTIMIZATION: Separate non-wrapping and wrapping paths
        if raw_end <= ref_length {
            // Non-wrapping fast path
            for p in raw_start..raw_end {
                if track_insert {
                    arrays.sum_insert_sizes[p] += tl;
                    arrays.count_insert_sizes[p] += 1;
                }
                if non_inward {
                    arrays.non_inward_pairs[p] += 1;
                }
                if mate_unmapped_flag {
                    arrays.mate_not_mapped[p] += 1;
                }
                if mate_other_contig_flag {
                    arrays.mate_on_another_contig[p] += 1;
                }
            }
        } else {
            // Wrapping case
            for pos in raw_start..raw_end {
                let p = pos % ref_length;
                if track_insert {
                    arrays.sum_insert_sizes[p] += tl;
                    arrays.count_insert_sizes[p] += 1;
                }
                if non_inward {
                    arrays.non_inward_pairs[p] += 1;
                }
                if mate_unmapped_flag {
                    arrays.mate_not_mapped[p] += 1;
                }
                if mate_other_contig_flag {
                    arrays.mate_on_another_contig[p] += 1;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Mapping metrics module - clippings, indels, mismatches (primary only)
    // Also needed for phagetermini (clippings used for completeness detection)
    // -------------------------------------------------------------------------
    if (flags.mapping_metrics || flags.phagetermini) && !is_secondary && !is_supplementary {
        // --- Clippings from CIGAR ---
        // Check first and last CIGAR operations for soft/hard clips
        // Left clippings recorded at first aligned position (start)
        // Right clippings recorded at last aligned position (end - 1)
        if !cigar_raw.is_empty() {
            // Check first operation for left clipping
            if let Some(&(op, len)) = cigar_raw.first() {
                if raw_cigar_is_clipping(op) {
                    // Record at first aligned position
                    arrays.left_clipping_lengths[start].push(len);
                    // Extract clip sequence (truncated to 20bp) for soft clips
                    if flags.mapping_metrics && op as u8 as char == 'S' && !seq.is_empty() {
                        let clip_len = (len as usize).min(20).min(seq.len());
                        track_sequence(&mut arrays.left_clip_sequences, start, &seq[..clip_len]);
                    }
                }
            }
            // Check last operation for right clipping
            if let Some(&(op, len)) = cigar_raw.last() {
                if raw_cigar_is_clipping(op) {
                    // Record at last aligned position
                    let clip_pos = if end > 0 {
                        let pos = end - 1;
                        if circular { pos % ref_length } else { pos }
                    } else {
                        ref_length - 1
                    };
                    arrays.right_clipping_lengths[clip_pos].push(len);
                    // Extract clip sequence (truncated to 20bp) for soft clips
                    if flags.mapping_metrics && op as u8 as char == 'S' && !seq.is_empty() {
                        let len_usize = len as usize;
                        let clip_start = seq.len().saturating_sub(len_usize);
                        let clip_len = len_usize.min(20);
                        let clip_end = (clip_start + clip_len).min(seq.len());
                        track_sequence(&mut arrays.right_clip_sequences, clip_pos, &seq[clip_start..clip_end]);
                    }
                }
            }
        }
    }

    // --- Unified CIGAR + MD + SEQ walk (mapping_metrics only) ---
    // Single walk handles indels, mismatches, and sequence extraction
    if flags.mapping_metrics && !is_secondary && !is_supplementary {
        let mut ref_pos = raw_start;
        let mut query_pos: usize = 0;

        // MD tag parsing state
        let md_bytes = md_tag.unwrap_or(&[]);
        let mut md_pos: usize = 0;
        let mut md_match_remaining: usize = 0;
        let has_md = !md_bytes.is_empty();
        let has_seq = !seq.is_empty();

        for &(op, len) in cigar_raw {
            let c = op as u8 as char;
            let len_usize = len as usize;

            match c {
                'M' | '=' | 'X' => {
                    // Match/mismatch region - walk MD tag in parallel
                    if has_md {
                        let mut remaining = len_usize;
                        while remaining > 0 {
                            if md_match_remaining > 0 {
                                // Consume matching bases
                                let consume = md_match_remaining.min(remaining);
                                ref_pos += consume;
                                query_pos += consume;
                                remaining -= consume;
                                md_match_remaining -= consume;
                            } else if md_pos < md_bytes.len() {
                                let b = md_bytes[md_pos];
                                if b.is_ascii_digit() {
                                    // Parse match count
                                    let mut num = 0usize;
                                    while md_pos < md_bytes.len() && md_bytes[md_pos].is_ascii_digit() {
                                        num = num * 10 + (md_bytes[md_pos] - b'0') as usize;
                                        md_pos += 1;
                                    }
                                    md_match_remaining = num;
                                } else if b.is_ascii_uppercase() {
                                    // Mismatch: record count and read base
                                    let normalized_pos = ref_pos % ref_length;
                                    arrays.mismatches[normalized_pos] += 1;

                                    if has_seq && query_pos < seq.len() {
                                        let read_base = seq[query_pos];
                                        let base_idx = match read_base {
                                            b'A' | b'a' => 0usize,
                                            b'C' | b'c' => 1,
                                            b'G' | b'g' => 2,
                                            b'T' | b't' => 3,
                                            _ => 4, // N or other - skip
                                        };
                                        if base_idx < 4 {
                                            arrays.mismatch_base_counts[normalized_pos][base_idx] += 1;
                                        }
                                    }

                                    ref_pos += 1;
                                    query_pos += 1;
                                    remaining -= 1;
                                    md_pos += 1;
                                } else if b == b'^' {
                                    // Deletion marker in MD - shouldn't happen inside M block
                                    // but handle gracefully
                                    break;
                                } else {
                                    md_pos += 1;
                                }
                            } else {
                                // MD exhausted, just advance
                                ref_pos += remaining;
                                query_pos += remaining;
                                remaining = 0;
                            }
                        }
                    } else {
                        // No MD tag, just advance positions
                        ref_pos += len_usize;
                        query_pos += len_usize;
                    }
                }
                'I' => {
                    // Insertion: extra sequence in read, not in reference
                    let normalized_pos = ref_pos % ref_length;
                    arrays.insertion_lengths[normalized_pos].push(len);

                    // Extract insertion sequence (capped to 40bp: first 20 + "..." + last 20)
                    if has_seq && query_pos + len_usize <= seq.len() {
                        if len_usize <= 40 {
                            track_sequence(&mut arrays.insertion_sequences, normalized_pos, &seq[query_pos..query_pos + len_usize]);
                        } else {
                            let mut v = seq[query_pos..query_pos + 20].to_vec();
                            v.extend_from_slice(b"...");
                            v.extend_from_slice(&seq[query_pos + len_usize - 20..query_pos + len_usize]);
                            track_sequence(&mut arrays.insertion_sequences, normalized_pos, &v);
                        }
                    }

                    query_pos += len_usize;
                }
                'D' => {
                    // Deletion: sequence in reference, not in read
                    let normalized_pos = ref_pos % ref_length;
                    arrays.deletion_lengths[normalized_pos].push(len);
                    for j in 0..len_usize {
                        arrays.deletions[(ref_pos + j) % ref_length] += 1;
                    }
                    ref_pos += len_usize;

                    // Skip deletion bases in MD tag (^ACGT...)
                    if md_pos < md_bytes.len() && md_bytes[md_pos] == b'^' {
                        md_pos += 1;
                        while md_pos < md_bytes.len() && md_bytes[md_pos].is_ascii_uppercase() {
                            md_pos += 1;
                        }
                    }
                }
                'S' => {
                    // Soft clip - only advances query position
                    query_pos += len_usize;
                }
                'N' => {
                    // Reference skip
                    ref_pos += len_usize;
                }
                'H' => {
                    // Hard clip - no sequence in BAM, skip
                }
                _ => {
                    // Other ops
                    if raw_cigar_consumes_ref(op) {
                        ref_pos += len_usize;
                    }
                }
            }
        }

    }


    // -------------------------------------------------------------------------
    // Phagetermini module - primary mappings only
    // -------------------------------------------------------------------------
    if flags.phagetermini && !is_secondary && !is_supplementary {
        // Compute boundary event lengths: Some(0) for exact match, Some(len) for near-match, None if neither
        let start_event = raw_boundary_event_length(cigar_raw, md_tag, !is_reverse, min_clipping_length);
        let start_matches = start_event.is_some();

        if seq_type.is_long() {
            // Long reads: split read in half and check each terminus independently
            // This avoids losing ~80% of reads that have clipping at one end
            let end_event = raw_boundary_event_length(cigar_raw, md_tag, is_reverse, min_clipping_length);
            let end_matches = end_event.is_some();
            let midpoint = (raw_start as usize + raw_end as usize) / 2;

            // Count clean termini for clipped_ratio calculation
            // Each clean terminus counts as 1 (a read with both clean = 2)
            if start_matches { arrays.clean_reads_count += 1; }
            if end_matches { arrays.clean_reads_count += 1; }

            // Coverage: split based on which termini match
            // - Both match: full read coverage
            // - Only start: coverage from start to midpoint
            // - Only end: coverage from midpoint to end
            match (start_matches, end_matches) {
                (true, true) => {
                    // Both termini match: count full read coverage
                    if circular {
                        increment_circular_long(&mut arrays.coverage_reduced, raw_start, raw_end, 1);
                    } else {
                        increment_range(&mut arrays.coverage_reduced, start, end, 1);
                    }
                }
                (true, false) => {
                    // Only start matches: count first half coverage
                    if circular {
                        increment_circular_long(&mut arrays.coverage_reduced, raw_start, midpoint, 1);
                    } else {
                        increment_range(&mut arrays.coverage_reduced, start, midpoint, 1);
                    }
                }
                (false, true) => {
                    // Only end matches: count second half coverage
                    if circular {
                        increment_circular_long(&mut arrays.coverage_reduced, midpoint, raw_end, 1);
                    } else {
                        increment_range(&mut arrays.coverage_reduced, midpoint % ref_length, end, 1);
                    }
                }
                (false, false) => {
                    // Neither matches: no coverage counted
                }
            }

            // Count start position and collect event length only if start matches
            if let Some(evt_len) = start_event {
                arrays.start_event_lengths[start].push(evt_len);
                if is_reverse {
                    arrays.start_minus[start] += 1;
                } else {
                    arrays.start_plus[start] += 1;
                }
                // Collect short clip sequence for reads_starts tooltip
                if evt_len > 0 {
                    if let Some(&(op, _)) = cigar_raw.first() {
                        if op as u8 as char == 'S' && !seq.is_empty() {
                            let clip_len = (evt_len as usize).min(20).min(seq.len());
                            track_sequence(&mut arrays.start_clip_sequences, start, &seq[..clip_len]);
                        }
                    }
                }
            }

            // Count end position and collect event length only if end matches
            if let Some(evt_len) = end_event {
                let end_pos = if end > 0 {
                    let pos = end - 1;
                    if circular { pos % ref_length } else { pos }
                } else {
                    0
                };

                arrays.end_event_lengths[end_pos].push(evt_len);
                if is_reverse {
                    arrays.end_minus[end_pos] += 1;
                } else {
                    arrays.end_plus[end_pos] += 1;
                }
                // Collect short clip sequence for reads_ends tooltip
                if evt_len > 0 {
                    if let Some(&(op, _)) = cigar_raw.last() {
                        if op as u8 as char == 'S' && !seq.is_empty() {
                            let clip_start = seq.len().saturating_sub(evt_len as usize);
                            let clip_len = (evt_len as usize).min(20);
                            let clip_end = (clip_start + clip_len).min(seq.len());
                            track_sequence(&mut arrays.end_clip_sequences, end_pos, &seq[clip_start..clip_end]);
                        }
                    }
                }
            }
        } else {
            // Short reads: only check start, count both positions if valid
            if start_matches {
                arrays.clean_reads_count += 1;

                if circular {
                    increment_circular(&mut arrays.coverage_reduced, start, end, 1);
                } else {
                    increment_range(&mut arrays.coverage_reduced, start, end, 1);
                }

                let end_pos = if end > 0 {
                    let pos = end - 1;
                    if circular { pos % ref_length } else { pos }
                } else {
                    0
                };

                // Collect start event length
                let start_evt_len = start_event.unwrap();
                arrays.start_event_lengths[start].push(start_evt_len);

                // Collect short clip sequence for reads_starts tooltip
                if start_evt_len > 0 {
                    if let Some(&(op, _)) = cigar_raw.first() {
                        if op as u8 as char == 'S' && !seq.is_empty() {
                            let clip_len = (start_evt_len as usize).min(20).min(seq.len());
                            track_sequence(&mut arrays.start_clip_sequences, start, &seq[..clip_len]);
                        }
                    }
                }

                // Compute and collect end event length separately
                let end_event = raw_boundary_event_length(cigar_raw, md_tag, is_reverse, min_clipping_length);
                let end_evt_len = end_event.unwrap_or(0);
                arrays.end_event_lengths[end_pos].push(end_evt_len);

                // Collect short clip sequence for reads_ends tooltip
                if end_evt_len > 0 {
                    if let Some(&(op, _)) = cigar_raw.last() {
                        if op as u8 as char == 'S' && !seq.is_empty() {
                            let clip_start = seq.len().saturating_sub(end_evt_len as usize);
                            let clip_len = (end_evt_len as usize).min(20);
                            let clip_end = (clip_start + clip_len).min(seq.len());
                            track_sequence(&mut arrays.end_clip_sequences, end_pos, &seq[clip_start..clip_end]);
                        }
                    }
                }

                if is_reverse {
                    arrays.start_minus[start] += 1;
                    arrays.end_minus[end_pos] += 1;
                } else {
                    arrays.start_plus[start] += 1;
                    arrays.end_plus[end_pos] += 1;
                }
            }
        }
    }
}
