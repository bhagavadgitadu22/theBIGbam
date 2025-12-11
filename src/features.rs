//! Feature calculation from aligned sequencing reads.
//!
//! Calculates coverage, phage termini features (reads_starts, reads_ends, tau),
//! and assembly check features (clippings, indels, mismatches, orientations).
//!
//! All features are calculated in a single pass through the BAM reads for efficiency.

use crate::cigar::{raw_cigar_consumes_ref, raw_cigar_is_clipping, raw_has_match_at_position, MdTag};
use crate::circular::{increment_circular, increment_range};
use crate::types::{FeatureMap, SequencingType};

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

    /// Number of secondary alignments at each position.
    /// Secondary reads have flag 0x100 set.
    pub secondary_reads: Vec<u64>,

    /// Number of supplementary alignments at each position.
    /// Supplementary reads have flag 0x800 set.
    pub supplementary_reads: Vec<u64>,

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

    // -------------------------------------------------------------------------
    // Assemblycheck Module
    // -------------------------------------------------------------------------
    /// Reads with soft/hard clipping at their 5' end.
    /// High values may indicate assembly breakpoints or structural variants.
    pub left_clippings: Vec<u64>,

    /// Reads with soft/hard clipping at their 3' end.
    pub right_clippings: Vec<u64>,

    /// Insertion events from CIGAR (sequence in read but not in reference).
    pub insertions: Vec<u64>,

    /// Deletion events from CIGAR (sequence in reference but not in read).
    pub deletions: Vec<u64>,

    /// Base mismatches (from MD tag).
    /// High mismatch rates may indicate errors or strain variation.
    pub mismatches: Vec<u64>,

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

    /// Reads where mates have unexpected orientation (not proper pairs).
    /// High values may indicate structural rearrangements.
    pub bad_orientations: Vec<u64>,

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
            secondary_reads: vec![0u64; ref_length],
            supplementary_reads: vec![0u64; ref_length],
            coverage_reduced: vec![0u64; ref_length],
            reads_starts: vec![0u64; ref_length],
            reads_ends: vec![0u64; ref_length],
            left_clippings: vec![0u64; ref_length],
            right_clippings: vec![0u64; ref_length],
            insertions: vec![0u64; ref_length],
            deletions: vec![0u64; ref_length],
            mismatches: vec![0u64; ref_length],
            sum_read_lengths: vec![0u64; ref_length],
            count_read_lengths: vec![0u64; ref_length],
            sum_insert_sizes: vec![0u64; ref_length],
            count_insert_sizes: vec![0u64; ref_length],
            bad_orientations: vec![0u64; ref_length],
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
                // Long reads: sum both strands
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

    /// Convert to FeatureMap for assemblycheck features.
    ///
    /// Only includes features relevant to the sequencing type:
    /// - Long reads: includes read_lengths
    /// - Paired reads: includes insert_sizes and bad_orientations
    pub fn to_assemblycheck_map(&self, seq_type: SequencingType) -> FeatureMap {
        let mut results = FeatureMap::with_capacity(10);
        results.insert("left_clippings".to_string(), self.left_clippings.clone());
        results.insert("right_clippings".to_string(), self.right_clippings.clone());
        results.insert("insertions".to_string(), self.insertions.clone());
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
            results.insert("bad_orientations".to_string(), self.bad_orientations.clone());
        }

        results
    }
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
    /// Calculate phage termini features (coverage_reduced, reads_starts, reads_ends)
    pub phagetermini: bool,
    /// Calculate assembly check features (clippings, indels, mismatches, etc.)
    pub assemblycheck: bool,
}

impl ModuleFlags {
    /// Create flags from a list of module names.
    ///
    /// # Arguments
    /// * `modules` - List of module names: "coverage", "phagetermini", "assemblycheck"
    ///
    /// # Note
    /// If assemblycheck is enabled, coverage is automatically enabled too
    /// (assemblycheck needs coverage data for normalization).
    pub fn from_modules(modules: &[String]) -> Self {
        let assemblycheck = modules.iter().any(|m| m == "assemblycheck");
        Self {
            // Coverage is needed if explicitly requested OR if assemblycheck is enabled
            coverage: modules.iter().any(|m| m == "coverage") || assemblycheck,
            phagetermini: modules.iter().any(|m| m == "phagetermini"),
            assemblycheck,
        }
    }

    /// Check if we need to extract the MD tag from reads.
    ///
    /// MD tag extraction has some overhead, so we only do it when needed.
    /// It's required for:
    /// - phagetermini: to check for mismatches at read ends
    /// - assemblycheck: to count mismatches
    #[inline]
    pub fn needs_md(&self) -> bool {
        self.phagetermini || self.assemblycheck
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
/// * `cigar_raw` - CIGAR operations as (operation_char, length) tuples
/// * `md_tag` - Optional MD tag bytes for mismatch detection
/// * `seq_type` - Sequencing type (affects which features to calculate)
/// * `flags` - Which modules are enabled
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
    cigar_raw: &[(u32, u32)],
    md_tag: Option<&[u8]>,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
) {
    // -------------------------------------------------------------------------
    // Calculate positions
    // -------------------------------------------------------------------------
    let ref_length = arrays.ref_length();
    let raw_start = ref_start as usize;
    let raw_end = ref_end as usize;

    // For circular genomes with doubled references, always apply modulo for array bounds.
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
                increment_circular(&mut arrays.secondary_reads, start, end, 1);
            } else {
                increment_range(&mut arrays.secondary_reads, start, end, 1);
            }
        } else if is_supplementary {
            if circular {
                increment_circular(&mut arrays.supplementary_reads, start, end, 1);
            } else {
                increment_range(&mut arrays.supplementary_reads, start, end, 1);
            }
        } else {
            // Only count primary mappings in main coverage
            if circular {
                increment_circular(&mut arrays.primary_reads, start, end, 1);
            } else {
                increment_range(&mut arrays.primary_reads, start, end, 1);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Phagetermini module - primary mappings only
    // -------------------------------------------------------------------------
    if flags.phagetermini && !is_secondary && !is_supplementary {
        // For long reads, we also check the end for clipping/mismatches
        let check_end = seq_type.is_long();

        // Check if read has clean alignment at start (no clip/mismatch)
        // raw_has_match_at_position checks both CIGAR and MD tag
        let start_matches = raw_has_match_at_position(cigar_raw, md_tag, true);
        let end_matches = !check_end || raw_has_match_at_position(cigar_raw, md_tag, false);

        // Only count "clean" primary reads for phage termini detection
        if start_matches && end_matches {
            // Update coverage_reduced (clean coverage from primary mappings only)
            // Note: end is exclusive (one past last position), so use non-inclusive increment
            if circular {
                increment_circular(&mut arrays.coverage_reduced, start, end, 1);
            } else {
                increment_range(&mut arrays.coverage_reduced, start, end, 1);
            }

            // Track start/end positions by strand
            // We separate strands here; they're combined in finalize_strands()
            // Note: 'end' is exclusive, so the actual last position is end-1
            let end_pos = if end > 0 { 
                let pos = end - 1;
                if circular { pos % ref_length } else { pos }
            } else { 
                0 
            };
            if is_reverse {
                arrays.start_minus[start] += 1;
                arrays.end_minus[end_pos] += 1;
            } else {
                arrays.start_plus[start] += 1;
                arrays.end_plus[end_pos] += 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Assemblycheck module - primary mappings only
    // -------------------------------------------------------------------------
    if flags.assemblycheck && !is_secondary && !is_supplementary {
        // --- Read lengths (long reads only) ---
        // Track sum and count to compute average read length at each position
        if seq_type.is_long() {
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

        // --- Insert sizes and bad orientations (short paired only) ---
        if seq_type.is_short_paired() {
            // Only track insert size for read1 with valid template length
            // (read2 would double-count the same fragment)
            let track_insert = is_read1 && template_length > 0;
            let tl = template_length as u64;
            let bad_orient = !is_proper_pair;

            // OPTIMIZATION: Hoist conditions outside loop and combine updates
            // This creates 4 versions of the loop, but each one is tight
            if raw_end <= ref_length {
                // Non-wrapping fast path
                if track_insert && bad_orient {
                    for p in raw_start..raw_end {
                        arrays.sum_insert_sizes[p] += tl;
                        arrays.count_insert_sizes[p] += 1;
                        arrays.bad_orientations[p] += 1;
                    }
                } else if track_insert {
                    for p in raw_start..raw_end {
                        arrays.sum_insert_sizes[p] += tl;
                        arrays.count_insert_sizes[p] += 1;
                    }
                } else if bad_orient {
                    for p in raw_start..raw_end {
                        arrays.bad_orientations[p] += 1;
                    }
                }
                // If neither track_insert nor bad_orient, nothing to do
            } else {
                // Wrapping case
                for pos in raw_start..raw_end {
                    let p = pos % ref_length;
                    if track_insert {
                        arrays.sum_insert_sizes[p] += tl;
                        arrays.count_insert_sizes[p] += 1;
                    }
                    if bad_orient {
                        arrays.bad_orientations[p] += 1;
                    }
                }
            }
        }

        // --- Clippings from CIGAR ---
        // Check first and last CIGAR operations for soft/hard clips
        if !cigar_raw.is_empty() {
            // Check first operation
            if let Some(&(op, _)) = cigar_raw.first() {
                if raw_cigar_is_clipping(op) {
                    arrays.left_clippings[start] += 1;
                }
            }
            // Check last operation
            if let Some(&(op, _)) = cigar_raw.last() {
                if raw_cigar_is_clipping(op) {
                    // End position is where the clip happens
                    // (subtract 1 because end is exclusive)
                    let clip_pos = if end > 0 { 
                        let pos = end - 1;
                        if circular { pos % ref_length } else { pos }
                    } else { 
                        ref_length - 1 
                    };
                    arrays.right_clippings[clip_pos] += 1;
                }
            }
        }

        // --- Indels from CIGAR ---
        // Walk through CIGAR operations to find insertions and deletions
        let mut ref_pos = raw_start;
        for &(op, len) in cigar_raw {
            let c = op as u8 as char;
            match c {
                'I' => {
                    // Insertion: extra sequence in read, not in reference
                    // Record at current reference position
                    arrays.insertions[ref_pos % ref_length] += 1;
                    // Insertions don't consume reference (ref_pos stays same)
                }
                'D' => {
                    // Deletion: sequence in reference, not in read
                    // Spans multiple reference positions
                    let len = len as usize;
                    for j in 0..len {
                        arrays.deletions[(ref_pos + j) % ref_length] += 1;
                    }
                    ref_pos += len; // Deletions consume reference
                }
                _ => {
                    // Other operations (M, =, X, N) that consume reference
                    if raw_cigar_consumes_ref(op) {
                        ref_pos += len as usize;
                    }
                    // S, H, P don't consume reference (ignored here)
                }
            }
        }

        // --- Mismatches from MD tag ---
        // The MD tag describes mismatches in the alignment
        if let Some(md_bytes) = md_tag {
            let md = MdTag::new(md_bytes);
            // Iterator-based: no Vec allocation, just yields positions
            for pos in md.mismatch_positions_normalized(raw_start, ref_length) {
                arrays.mismatches[pos] += 1;
            }
        }
    }
}
