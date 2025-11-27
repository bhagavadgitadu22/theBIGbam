//! Feature calculation functions.
//!
//! This module contains all the feature calculation logic for coverage,
//! phage termini, and assembly check features.
//!
//! Optimized for single-pass processing to avoid multiple iterations over reads.

use crate::cigar::{raw_cigar_consumes_ref, raw_cigar_is_clipping, raw_has_match_at_position, MdTag};
use crate::circular::CircularArray;
use crate::types::{FeatureMap, SequencingType};

/// All feature arrays calculated in a single pass.
pub struct FeatureArrays {
    // Coverage
    pub coverage: Vec<u64>,

    // Phagetermini
    pub coverage_reduced: Vec<u64>,
    pub reads_starts: Vec<u64>,
    pub reads_ends: Vec<u64>,

    // Assemblycheck
    pub left_clippings: Vec<u64>,
    pub right_clippings: Vec<u64>,
    pub insertions: Vec<u64>,
    pub deletions: Vec<u64>,
    pub mismatches: Vec<u64>,
    pub sum_read_lengths: Vec<u64>,
    pub count_read_lengths: Vec<u64>,
    pub sum_insert_sizes: Vec<u64>,
    pub count_insert_sizes: Vec<u64>,
    pub bad_orientations: Vec<u64>,

    // Internal tracking for phagetermini strand separation
    start_plus: Vec<u64>,
    start_minus: Vec<u64>,
    end_plus: Vec<u64>,
    end_minus: Vec<u64>,

    // Coverage percentage tracking
    covered: Vec<bool>,
}

impl FeatureArrays {
    /// Create new feature arrays for a given reference length.
    #[inline]
    pub fn new(ref_length: usize) -> Self {
        Self {
            coverage: vec![0u64; ref_length],
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
            covered: vec![false; ref_length],
        }
    }

    /// Get reference length.
    #[inline]
    pub fn ref_length(&self) -> usize {
        self.coverage.len()
    }

    /// Calculate coverage percentage from the covered bitmap.
    #[inline]
    pub fn coverage_percentage(&self) -> f64 {
        let covered_bp = self.covered.iter().filter(|&&x| x).count();
        (covered_bp as f64 / self.ref_length() as f64) * 100.0
    }

    /// Finalize phagetermini strand arrays based on sequencing type.
    pub fn finalize_strands(&mut self, seq_type: SequencingType) {
        match seq_type {
            SequencingType::ShortPaired | SequencingType::ShortSingle => {
                // Short reads: starts from + strand, ends from - strand
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
    pub fn to_phagetermini_map(&self) -> FeatureMap {
        let mut results = FeatureMap::with_capacity(3);
        results.insert("coverage_reduced".to_string(), self.coverage_reduced.clone());
        results.insert("reads_starts".to_string(), self.reads_starts.clone());
        results.insert("reads_ends".to_string(), self.reads_ends.clone());
        results
    }

    /// Convert to FeatureMap for assemblycheck features.
    pub fn to_assemblycheck_map(&self, seq_type: SequencingType) -> FeatureMap {
        let mut results = FeatureMap::with_capacity(10);
        results.insert("left_clippings".to_string(), self.left_clippings.clone());
        results.insert("right_clippings".to_string(), self.right_clippings.clone());
        results.insert("insertions".to_string(), self.insertions.clone());
        results.insert("deletions".to_string(), self.deletions.clone());
        results.insert("mismatches".to_string(), self.mismatches.clone());

        if seq_type.is_long() {
            results.insert("sum_read_lengths".to_string(), self.sum_read_lengths.clone());
            results.insert("count_read_lengths".to_string(), self.count_read_lengths.clone());
        }

        if seq_type.is_short_paired() {
            results.insert("sum_insert_sizes".to_string(), self.sum_insert_sizes.clone());
            results.insert("count_insert_sizes".to_string(), self.count_insert_sizes.clone());
            results.insert("bad_orientations".to_string(), self.bad_orientations.clone());
        }

        results
    }
}

/// Configuration for which modules to process.
#[derive(Clone, Copy)]
pub struct ModuleFlags {
    pub coverage: bool,
    pub phagetermini: bool,
    pub assemblycheck: bool,
}

impl ModuleFlags {
    pub fn from_modules(modules: &[String]) -> Self {
        let assemblycheck = modules.iter().any(|m| m == "assemblycheck");
        Self {
            coverage: modules.iter().any(|m| m == "coverage") || assemblycheck,
            phagetermini: modules.iter().any(|m| m == "phagetermini"),
            assemblycheck,
        }
    }

    #[inline]
    pub fn needs_md(&self) -> bool {
        self.phagetermini || self.assemblycheck
    }
}

/// Process a single read and update all feature arrays in one pass.
///
/// This is the core optimized function - processes each read exactly once.
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
    cigar_raw: &[(u32, u32)],
    md_tag: Option<&[u8]>,
    seq_type: SequencingType,
    flags: ModuleFlags,
) {
    let ref_length = arrays.ref_length();
    let raw_start = ref_start as usize;
    let raw_end = ref_end as usize;
    let start = raw_start % ref_length;
    let end = raw_end % ref_length;

    // Update covered bitmap using efficient slice fill (always needed for coverage percentage)
    let cov_start = raw_start.min(ref_length);
    let cov_end = raw_end.min(ref_length);
    arrays.covered[cov_start..cov_end].fill(true);

    // Coverage (simple increment)
    if flags.coverage {
        arrays.coverage.increment_circular(start, end, 1);
    }

    // Phagetermini processing
    if flags.phagetermini {
        let check_end = seq_type.is_long();

        // Check if read starts/ends with match using raw CIGAR and MD tag (zero-allocation)
        let start_matches = raw_has_match_at_position(cigar_raw, md_tag, true);
        let end_matches = !check_end || raw_has_match_at_position(cigar_raw, md_tag, false);

        if start_matches && end_matches {
            arrays.coverage_reduced.increment_circular_inclusive(start, end, 1);

            if is_reverse {
                arrays.start_minus[start] += 1;
                arrays.end_minus[end] += 1;
            } else {
                arrays.start_plus[start] += 1;
                arrays.end_plus[end] += 1;
            }
        }
    }

    // Assemblycheck processing - combined loops for efficiency
    if flags.assemblycheck {
        // Read lengths (long reads) - single loop
        if seq_type.is_long() {
            let ql = query_length as u64;
            // Non-wrapping case (most common)
            if raw_end <= ref_length {
                for p in raw_start..raw_end {
                    arrays.sum_read_lengths[p] += ql;
                    arrays.count_read_lengths[p] += 1;
                }
            } else {
                // Wrapping case
                for pos in raw_start..raw_end {
                    let p = pos % ref_length;
                    arrays.sum_read_lengths[p] += ql;
                    arrays.count_read_lengths[p] += 1;
                }
            }
        }

        // Insert sizes and bad orientations (short-paired) - combined single loop
        if seq_type.is_short_paired() {
            let track_insert = is_read1 && template_length > 0;
            let tl = template_length as u64;
            let bad_orient = !is_proper_pair;

            // Non-wrapping case (most common for short reads)
            if raw_end <= ref_length {
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

        // Clippings from raw CIGAR (zero-allocation)
        if !cigar_raw.is_empty() {
            if let Some(&(op, _)) = cigar_raw.first() {
                if raw_cigar_is_clipping(op) {
                    arrays.left_clippings[start] += 1;
                }
            }
            if let Some(&(op, _)) = cigar_raw.last() {
                if raw_cigar_is_clipping(op) {
                    arrays.right_clippings[if end > 0 { end - 1 } else { ref_length - 1 }] += 1;
                }
            }
        }

        // Indels from raw CIGAR (zero-allocation)
        let mut ref_pos = raw_start;
        for &(op, len) in cigar_raw {
            let c = op as u8 as char;
            match c {
                'I' => {
                    arrays.insertions[ref_pos % ref_length] += 1;
                }
                'D' => {
                    let len = len as usize;
                    for j in 0..len {
                        arrays.deletions[(ref_pos + j) % ref_length] += 1;
                    }
                    ref_pos += len;
                }
                _ => {
                    if raw_cigar_consumes_ref(op) {
                        ref_pos += len as usize;
                    }
                }
            }
        }

        // Mismatches from MD tag (iterator-based, zero-allocation)
        if let Some(md_bytes) = md_tag {
            let md = MdTag::new(md_bytes);
            for pos in md.mismatch_positions_normalized(raw_start, ref_length) {
                arrays.mismatches[pos] += 1;
            }
        }
    }
}
