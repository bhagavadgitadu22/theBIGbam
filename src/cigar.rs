//! CIGAR and MD tag parsing for BAM alignment analysis.
//!
//! CIGAR strings describe how a read aligns to reference (e.g., "50M2I30M"):
//! - M: Match, I: Insertion, D: Deletion, S: Soft clip, H: Hard clip
//!
//! MD tags provide mismatch information (e.g., "10A5^GT3"):
//! - Numbers: matching bases
//! - Letters: mismatches
//! - ^: deletions

// ============================================================================
// CIGAR OPERATION TYPES
// ============================================================================

/// CIGAR operation types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Padding,
    SeqMatch,
    SeqMismatch,
}

impl CigarOp {
    #[inline]
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            'M' => Some(Self::Match),
            'I' => Some(Self::Insertion),
            'D' => Some(Self::Deletion),
            'N' => Some(Self::Skip),
            'S' => Some(Self::SoftClip),
            'H' => Some(Self::HardClip),
            'P' => Some(Self::Padding),
            '=' => Some(Self::SeqMatch),
            'X' => Some(Self::SeqMismatch),
            _ => None,
        }
    }

    #[inline]
    pub fn from_u32(op: u32) -> Option<Self> {
        Self::from_char(op as u8 as char)
    }

    #[inline]
    pub fn is_clipping(&self) -> bool {
        matches!(self, Self::SoftClip | Self::HardClip)
    }

    /// Check if this operation consumes the reference.
    ///
    /// "Consumes reference" means we advance along the reference genome.
    /// - M/=/X: move along both read and reference
    /// - D/N: move along reference only (gap in read)
    /// - I/S/H: don't move along reference
    ///
    /// This is important for calculating coverage positions.
    #[inline]
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Deletion | Self::Skip | Self::SeqMatch | Self::SeqMismatch
        )
    }

}

// ============================================================================
// MD TAG PARSING
// ============================================================================
// The MD tag is an optional BAM tag that encodes mismatch/deletion information.
// It allows us to determine base-level differences without the reference sequence.

/// MD tag parser for extracting mismatch positions.
///
/// # Bioinformatics Context
/// The MD (Mismatch Descriptor) tag shows where the read differs from reference:
/// - "10A5" = 10 matches, then mismatch (ref had A), then 5 matches
/// - "5^AC3" = 5 matches, then deletion of AC from reference, then 3 matches
///
/// The format is: `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`
/// - Numbers indicate consecutive matching bases
/// - Single uppercase letters indicate mismatches (the letter is the reference base)
/// - ^[A-Z]+ indicates deletions (letters are the deleted reference bases)
///
/// # Rust Concept: Lifetimes (`'a`)
/// `MdTag<'a>` has a lifetime parameter because it *borrows* the MD tag bytes
/// instead of copying them. The `'a` says "this MdTag lives only as long as
/// the byte slice it references". This is a zero-copy optimization.
///
/// In Python, you'd just store a reference and trust the GC. In Rust, lifetimes
/// make this explicit and guarantee memory safety at compile time.
pub struct MdTag<'a> {
    bytes: &'a [u8],  // Borrowed slice of bytes - no allocation, no copy!
}

impl<'a> MdTag<'a> {
    /// Create a new MD tag parser from raw bytes.
    ///
    /// # Rust Concept: Lifetime propagation
    /// Because we take `&'a [u8]` and return `Self` (which is `MdTag<'a>`),
    /// Rust knows the returned MdTag can only live as long as the input slice.
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { bytes }
    }

    /// Check if the tag indicates a match at the start.
    ///
    /// Returns true if the MD tag starts with a non-zero number (meaning there
    /// are matching bases before any mismatch or deletion).
    ///
    /// Used by PhageTerm to detect if read truly aligns at its start position.
    ///
    /// # Rust Concept: b'0' syntax
    /// `b'0'` is a byte literal (u8), while `'0'` is a char (4 bytes).
    /// We use byte literals because MD tags are stored as bytes.
    #[inline]
    pub fn has_match_at_start(&self) -> bool {
        !self.bytes.is_empty() && self.bytes[0].is_ascii_digit() && self.bytes[0] != b'0'
    }

    /// Check if the tag indicates a match at the end.
    ///
    /// Returns true if the MD tag ends with a non-zero number.
    #[inline]
    pub fn has_match_at_end(&self) -> bool {
        if self.bytes.is_empty() {
            return false;
        }
        // Find where the trailing number starts (parse backwards)
        let mut i = self.bytes.len();
        while i > 0 && self.bytes[i - 1].is_ascii_digit() {
            i -= 1;
        }
        // If we found digits at the end, check if the number > 0
        if i < self.bytes.len() {
            // Number > 0 if any digit is non-zero (handles "50", "100", etc.)
            self.bytes[i..].iter().any(|&b| b != b'0')
        } else {
            // No digits at end = ends with mismatch/deletion
            false
        }
    }

    /// Iterate over mismatch positions relative to reference start.
    ///
    /// Yields (reference_offset, mismatched_base) for each mismatch.
    /// Example: MD tag "5A3T2" yields [(5, 'A'), (9, 'T')]
    ///
    /// # Rust Concept: Lazy iterators
    /// This returns an iterator that computes positions on-demand.
    /// We don't allocate a vector of all mismatches upfront - we yield
    /// them one at a time as the caller iterates. Memory efficient!
    /// Iterate over mismatch positions normalized to circular genome.
    ///
    /// For circular genomes, positions wrap around using modulo.
    /// Used to count mismatches at each position on the circular reference.
    pub fn mismatch_positions_normalized(
        &self,
        ref_start: usize,
        ref_length: usize,
    ) -> MdMismatchNormalizedIter<'_> {
        MdMismatchNormalizedIter {
            bytes: self.bytes,
            pos: 0,
            ref_pos: ref_start,
            ref_length,
        }
    }
}

// ============================================================================
// MD TAG ITERATORS
// ============================================================================
// These custom iterators let us lazily yield mismatch positions from MD tags.
// They implement Rust's Iterator trait, enabling use with for loops and
// all iterator adapters (.map(), .filter(), .collect(), etc.)

/// Iterator over normalized mismatch positions in an MD tag.
///
/// Yields positions modulo `ref_length` for circular genome support.
/// This means position 105 on a genome of length 100 becomes position 5.
///
/// # Rust Concept: Iterator State Machine
/// An iterator in Rust is a struct that keeps track of its state between calls.
/// Each call to `next()` advances the state and yields the next item.
/// This is similar to Python generators but more explicit.
pub struct MdMismatchNormalizedIter<'a> {
    bytes: &'a [u8],   // The MD tag bytes we're parsing
    pos: usize,        // Current position in the byte slice
    ref_pos: usize,    // Current position along the reference genome
    ref_length: usize, // Length of reference (for modulo wrap-around)
}

impl<'a> Iterator for MdMismatchNormalizedIter<'a> {
    // What type does this iterator yield? Just the position (usize).
    type Item = usize;

    /// Get the next mismatch position.
    ///
    /// # Rust Concept: impl Iterator
    /// By implementing the `Iterator` trait with its `next()` method, our struct
    /// can be used in `for` loops and with all iterator methods:
    /// ```ignore
    /// for pos in md_tag.mismatch_positions_normalized(0, 100) {
    ///     // pos is yielded one at a time
    /// }
    /// ```
    fn next(&mut self) -> Option<Self::Item> {
        // Parse through the MD tag bytes until we find a mismatch or reach end
        while self.pos < self.bytes.len() {
            let c = self.bytes[self.pos];

            if c.is_ascii_digit() {
                // ----------------------------------------------------------
                // MATCH RUN: parse the number of matching bases
                // Example: "10" means 10 consecutive matches
                // ----------------------------------------------------------
                let mut num = 0usize;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                    // Convert ASCII digit to number: '5' -> 5
                    num = num * 10 + (self.bytes[self.pos] - b'0') as usize;
                    self.pos += 1;
                }
                self.ref_pos += num;  // Skip past the matching region
            } else if c == b'^' {
                // ----------------------------------------------------------
                // DELETION: skip the deleted bases (they're in reference, not read)
                // Example: "^AC" means reference has AC that read doesn't
                // ----------------------------------------------------------
                self.pos += 1;  // Skip the '^' character
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_uppercase() {
                    self.ref_pos += 1;  // Each deleted base advances reference position
                    self.pos += 1;
                }
            } else if c.is_ascii_uppercase() {
                // ----------------------------------------------------------
                // MISMATCH: found one! Yield its position
                // Example: "A" means reference had 'A' but read has something else
                // ----------------------------------------------------------
                let result = self.ref_pos % self.ref_length;  // Normalize for circular
                self.ref_pos += 1;
                self.pos += 1;
                return Some(result);  // Yield this position to the caller
            } else {
                // Skip any unexpected characters
                self.pos += 1;
            }
        }
        None  // End of MD tag - no more mismatches
    }
}

// ============================================================================
// RAW CIGAR HELPERS (ZERO-ALLOCATION)
// ============================================================================
// These functions work directly with rust-htslib's raw CIGAR format: &[(u32, u32)]
// where each tuple is (operation_char, length). They operate on slices directly
// with zero allocation, which matters in hot paths (called millions of times).

/// Check if raw CIGAR starts/ends with a match (no clipping/insertion).
///
/// Works directly on rust-htslib's raw format to avoid allocation.
#[inline]
pub fn raw_cigar_starts_with_match(cigar: &[(u32, u32)], at_start: bool) -> bool {
    let elem = if at_start { cigar.first() } else { cigar.last() };
    match elem {
        Some(&(op, _)) => {
            let c = op as u8 as char;
            // True match = not clipping (S/H) and not insertion (I)
            !matches!(c, 'S' | 'H' | 'I')
        }
        None => false,
    }
}

/// Check if a raw CIGAR operation is clipping (S or H).
#[inline]
pub fn raw_cigar_is_clipping(op: u32) -> bool {
    let c = op as u8 as char;
    matches!(c, 'S' | 'H')
}

/// Check if a raw CIGAR operation consumes the reference.
///
/// Operations that advance along the reference: M, D, N, =, X
/// Operations that don't: I, S, H, P
#[inline]
pub fn raw_cigar_consumes_ref(op: u32) -> bool {
    let c = op as u8 as char;
    matches!(c, 'M' | 'D' | 'N' | '=' | 'X')
}

/// Tolerant version of `raw_has_match_at_position` that treats small clips/insertions
/// (shorter than `min_clipping_length`) as near-matches for phagetermini.
///
/// If the read has an exact match, returns true immediately.
/// Otherwise, if the boundary operation is a small clip (S/H) or insertion (I),
/// checks the MD tag to see if the adjacent aligned bases match.
#[inline]
pub fn raw_has_near_match_at_position(
    cigar: &[(u32, u32)],
    md: Option<&[u8]>,
    at_start: bool,
    min_clipping_length: u32,
) -> bool {
    // Exact match passes immediately
    if raw_has_match_at_position(cigar, md, at_start) {
        return true;
    }

    if cigar.is_empty() {
        return false;
    }

    // Check if the boundary operation is a small clip/insertion
    if at_start {
        let (op, len) = cigar[0];
        let c = op as u8 as char;
        if !matches!(c, 'S' | 'H' | 'I') || len >= min_clipping_length {
            return false;
        }
        // Need at least one more operation after the clip
        if cigar.len() < 2 {
            return false;
        }
        // Next operation must not be another clip/insertion
        let (next_op, _) = cigar[1];
        let nc = next_op as u8 as char;
        if matches!(nc, 'S' | 'H' | 'I') {
            return false;
        }
        // Check MD tag for match at start
        match md {
            Some(bytes) if !bytes.is_empty() => MdTag::new(bytes).has_match_at_start(),
            _ => false,
        }
    } else {
        let (op, len) = cigar[cigar.len() - 1];
        let c = op as u8 as char;
        if !matches!(c, 'S' | 'H' | 'I') || len >= min_clipping_length {
            return false;
        }
        if cigar.len() < 2 {
            return false;
        }
        let (prev_op, _) = cigar[cigar.len() - 2];
        let pc = prev_op as u8 as char;
        if matches!(pc, 'S' | 'H' | 'I') {
            return false;
        }
        match md {
            Some(bytes) if !bytes.is_empty() => MdTag::new(bytes).has_match_at_end(),
            _ => false,
        }
    }
}

/// Return the boundary event length for a read at its start or end.
///
/// Returns `Some(0)` for exact match, `Some(len)` for a small clip/insertion
/// (shorter than `min_clipping_length`), or `None` if neither.
/// This mirrors `raw_has_near_match_at_position` but returns the actual length
/// instead of just a boolean, enabling per-position statistics.
#[inline]
pub fn raw_boundary_event_length(
    cigar: &[(u32, u32)],
    md: Option<&[u8]>,
    at_start: bool,
    min_clipping_length: u32,
) -> Option<u32> {
    // Exact match → event length 0
    if raw_has_match_at_position(cigar, md, at_start) {
        return Some(0);
    }

    if cigar.is_empty() {
        return None;
    }

    // Check if the boundary operation is a small clip/insertion
    if at_start {
        let (op, len) = cigar[0];
        let c = op as u8 as char;
        if !matches!(c, 'S' | 'H' | 'I') || len >= min_clipping_length {
            return None;
        }
        if cigar.len() < 2 {
            return None;
        }
        let (next_op, _) = cigar[1];
        let nc = next_op as u8 as char;
        if matches!(nc, 'S' | 'H' | 'I') {
            return None;
        }
        match md {
            Some(bytes) if !bytes.is_empty() => {
                if MdTag::new(bytes).has_match_at_start() { Some(len) } else { None }
            }
            _ => None,
        }
    } else {
        let (op, len) = cigar[cigar.len() - 1];
        let c = op as u8 as char;
        if !matches!(c, 'S' | 'H' | 'I') || len >= min_clipping_length {
            return None;
        }
        if cigar.len() < 2 {
            return None;
        }
        let (prev_op, _) = cigar[cigar.len() - 2];
        let pc = prev_op as u8 as char;
        if matches!(pc, 'S' | 'H' | 'I') {
            return None;
        }
        match md {
            Some(bytes) if !bytes.is_empty() => {
                if MdTag::new(bytes).has_match_at_end() { Some(len) } else { None }
            }
            _ => None,
        }
    }
}

/// Check if a read starts/ends with match using raw CIGAR and MD tag (zero-allocation).
///
/// This is the optimized version used in hot paths. Same logic as
/// `has_match_at_position` but works on raw slices to avoid allocation.
#[inline]
pub fn raw_has_match_at_position(cigar: &[(u32, u32)], md: Option<&[u8]>, at_start: bool) -> bool {
    // Check CIGAR first (cheap check)
    if !raw_cigar_starts_with_match(cigar, at_start) {
        return false;
    }

    // Then check MD tag
    match md {
        Some(bytes) if !bytes.is_empty() => {
            let md_tag = MdTag::new(bytes);
            if at_start {
                md_tag.has_match_at_start()
            } else {
                md_tag.has_match_at_end()
            }
        }
        _ => false,  // No MD tag or empty = can't confirm
    }
}
