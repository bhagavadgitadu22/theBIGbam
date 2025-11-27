//! CIGAR and MD tag parsing for BAM alignment analysis.
//!
//! # Overview
//!
//! ## CIGAR Strings
//! CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings describe how a
//! sequencing read aligns to a reference genome. A string like "50M2I30M" means:
//! - 50 bases Match/align to reference
//! - 2 bases Inserted in the read (not in reference)
//! - 30 bases Match/align to reference
//!
//! Common CIGAR operations:
//! - M: alignment Match (may include mismatches!)
//! - I: Insertion to the reference (bases in read, not in reference)
//! - D: Deletion from the reference (bases in reference, not in read)
//! - S: Soft clip (bases in read sequence but not aligned)
//! - H: Hard clip (bases not even in the read sequence)
//! - N: Skipped region (like introns in RNA-seq)
//! - =: Sequence match (exact match, unlike M)
//! - X: Sequence mismatch
//!
//! ## MD Tags
//! The MD tag provides mismatch and deletion information without needing the
//! reference sequence. Format: `10A5^GT3` means:
//! - 10 matching bases
//! - Mismatch to 'A'
//! - 5 matching bases
//! - Deletion of 'GT'
//! - 3 matching bases
//!
//! # Python Equivalent
//!
//! This module corresponds to the `starts_with_match()` function and MD tag
//! parsing logic in `calculating_data.py:301-315`:
//! ```python
//! def starts_with_match(cigar, md, start):
//!     op, length = cigar[0] if start else cigar[-1]
//!     if op in (4, 5):  # soft or hard clip
//!         return False
//!     if op == 1:  # insertion
//!         return False
//!     val = md[0] if start else md[-1]
//!     return val > 0
//! ```
//!
//! # Rust Concepts
//!
//! ## Enums with Derive Macros
//! `#[derive(Clone, Copy, Debug)]` automatically implements traits:
//! - `Clone`: enables `.clone()` to create a copy
//! - `Copy`: enables implicit copying (for small types like enums)
//! - `Debug`: enables `{:?}` formatting for debugging
//!
//! ## Lifetimes (`'a`)
//! `MdTag<'a>` has a lifetime parameter because it borrows data rather than
//! owning it. The `'a` indicates the struct lives only as long as the data it
//! references. This avoids copying the MD tag bytes.
//!
//! ## Custom Iterators
//! `MdMismatchIter` implements the `Iterator` trait to lazily yield mismatch
//! positions. This is memory-efficient: positions are computed on-demand rather
//! than building a complete list upfront.

// ============================================================================
// CIGAR OPERATION TYPES
// ============================================================================

/// CIGAR operation types - each represents a different alignment event.
///
/// # Rust Concept: Enums
/// In Rust, `enum` defines a type that can be one of several variants.
/// Unlike Python's string-based approach ("M", "I", etc.), this is type-safe:
/// the compiler ensures you handle all cases.
///
/// # Derive Macros Explained
/// - `Clone, Copy`: These are small (1 byte), so we can copy them freely
/// - `Debug`: Enables printing with `{:?}` for debugging
/// - `PartialEq, Eq`: Enables `==` comparison between operations
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    // Each variant maps to a single character in the CIGAR string
    Match,          // M - alignment match (can be sequence match or mismatch)
    Insertion,      // I - insertion to the reference
    Deletion,       // D - deletion from the reference
    Skip,           // N - skipped region from the reference (e.g., introns)
    SoftClip,       // S - soft clipping (sequence present in SEQ but not aligned)
    HardClip,       // H - hard clipping (sequence NOT present in SEQ)
    Padding,        // P - padding (silent deletion from padded reference)
    SeqMatch,       // = - sequence match (exact match, more specific than M)
    SeqMismatch,    // X - sequence mismatch (more specific than M)
}

// ============================================================================
// Rust Concept: impl blocks
// ============================================================================
// `impl CigarOp { ... }` adds methods to the CigarOp enum. In Python, these
// would be class methods. In Rust, the data definition (the enum above) is
// separated from its methods (the impl block below).

impl CigarOp {
    /// Parse a CIGAR operation from its character representation.
    ///
    /// # Rust Concept: #[inline]
    /// The `#[inline]` attribute suggests the compiler should inline this function
    /// at call sites (copy the code directly instead of making a function call).
    /// This is a performance optimization for small, frequently-called functions.
    ///
    /// # Rust Concept: Option<T>
    /// Returns `Some(CigarOp)` if valid, `None` if unknown character.
    /// This is safer than Python's approach of raising exceptions or returning None.
    #[inline]
    pub fn from_char(c: char) -> Option<Self> {
        // Rust Concept: match expression
        // Like Python's match/case but exhaustive - compiler ensures all cases handled
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
            _ => None,  // Unknown character - return None
        }
    }

    /// Parse from u32 (as stored internally from rust-htslib char()).
    ///
    /// rust-htslib stores CIGAR operations as u32 where the operation character
    /// is encoded in the lower bits. We extract it and convert.
    #[inline]
    pub fn from_u32(op: u32) -> Option<Self> {
        // Cast chain: u32 -> u8 (truncate to byte) -> char
        Self::from_char(op as u8 as char)
    }

    /// Check if this operation is a clipping operation (soft or hard).
    ///
    /// Clipping indicates bases at read ends that don't align:
    /// - SoftClip: bases present in read sequence but not used in alignment
    /// - HardClip: bases not even included in the read sequence
    ///
    /// # Rust Concept: matches! macro
    /// `matches!(value, pattern1 | pattern2)` returns true if value matches
    /// any of the patterns. More concise than a match expression that returns bool.
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

    /// Check if this operation consumes the query (read sequence).
    ///
    /// "Consumes query" means we advance along the read sequence.
    /// Used for calculating read length and sequence positions.
    #[inline]
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Insertion | Self::SoftClip | Self::SeqMatch | Self::SeqMismatch
        )
    }
}

// ============================================================================
// CIGAR ELEMENT (SINGLE OPERATION WITH LENGTH)
// ============================================================================

/// A single CIGAR operation with its length.
///
/// For example, in "50M2I30M":
/// - First element: CigarElement { op: Match, len: 50 }
/// - Second element: CigarElement { op: Insertion, len: 2 }
/// - Third element: CigarElement { op: Match, len: 30 }
#[derive(Clone, Copy, Debug)]
pub struct CigarElement {
    pub op: CigarOp,  // The operation type (M, I, D, etc.)
    pub len: u32,     // Number of bases this operation covers
}

impl CigarElement {
    /// Create from raw (op_char, length) tuple as returned by rust-htslib.
    ///
    /// # Rust Concept: Option::map()
    /// `option.map(|x| transform(x))` transforms the inner value if Some,
    /// leaves None unchanged. Like: `x.map(f)` == `if x.is_some() { Some(f(x.unwrap())) } else { None }`
    #[inline]
    pub fn from_raw(op: u32, len: u32) -> Option<Self> {
        // If from_u32 returns Some(op), transform it into Some(CigarElement)
        CigarOp::from_u32(op).map(|op| Self { op, len })
    }
}

// ============================================================================
// COMPLETE CIGAR STRING
// ============================================================================

/// A complete CIGAR string as a sequence of elements.
///
/// # Rust Concept: Tuple Struct
/// `struct Cigar(pub Vec<CigarElement>)` is a "tuple struct" - it's like a struct
/// but with unnamed fields accessed by index (self.0). Used when there's only
/// one field and we want a distinct type for type safety.
///
/// # Rust Concept: Default Trait
/// `#[derive(Default)]` gives us `Cigar::default()` which returns an empty Cigar.
/// Useful for initializing before filling with data.
#[derive(Clone, Debug, Default)]
pub struct Cigar(pub Vec<CigarElement>);

impl Cigar {
    /// Create from raw tuples of (op_char, length) as returned by rust-htslib.
    ///
    /// # Rust Concept: filter_map()
    /// `filter_map` combines filter and map: it takes a closure returning Option,
    /// keeps only the Some values, and unwraps them. Invalid operations are silently skipped.
    pub fn from_raw(raw: &[(u32, u32)]) -> Self {
        Self(
            raw.iter()
                .filter_map(|&(op, len)| CigarElement::from_raw(op, len))
                .collect(),
        )
    }

    /// Check if empty (no CIGAR operations).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()  // self.0 accesses the inner Vec
    }

    /// Get the first operation (or None if empty).
    ///
    /// # Rust Concept: Option<&T>
    /// Returns a reference to the element, not a copy. We don't need ownership,
    /// just to look at it. This avoids unnecessary copying.
    #[inline]
    pub fn first(&self) -> Option<&CigarElement> {
        self.0.first()
    }

    /// Get the last operation (or None if empty).
    #[inline]
    pub fn last(&self) -> Option<&CigarElement> {
        self.0.last()
    }

    /// Check if the first operation is a clipping (soft or hard clip).
    ///
    /// # Rust Concept: method chaining with Option
    /// `.map(|e| ...)` transforms Some(element) -> Some(bool)
    /// `.unwrap_or(false)` extracts the bool, defaulting to false if None
    #[inline]
    pub fn starts_with_clipping(&self) -> bool {
        self.first().map(|e| e.op.is_clipping()).unwrap_or(false)
    }

    /// Check if the last operation is a clipping.
    #[inline]
    pub fn ends_with_clipping(&self) -> bool {
        self.last().map(|e| e.op.is_clipping()).unwrap_or(false)
    }

    /// Check if an end starts with a match (not clipped, not insertion).
    ///
    /// This is used for PhageTerm analysis: we want to know if the read actually
    /// aligns at its start/end position, or if it's clipped/inserted there.
    ///
    /// `at_start`: if true, check the first op; if false, check the last op.
    pub fn starts_with_match(&self, at_start: bool) -> bool {
        let elem = if at_start { self.first() } else { self.last() };
        match elem {
            Some(e) => !matches!(e.op, CigarOp::SoftClip | CigarOp::HardClip | CigarOp::Insertion),
            None => false,
        }
    }

    /// Iterate over operations.
    ///
    /// # Rust Concept: impl Iterator
    /// `-> impl Iterator<Item = &CigarElement>` means "returns something that
    /// implements Iterator". The caller doesn't need to know the exact type.
    pub fn iter(&self) -> impl Iterator<Item = &CigarElement> {
        self.0.iter()
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
        !self.bytes.is_empty() && self.bytes.last().map(|&b| b.is_ascii_digit() && b != b'0').unwrap_or(false)
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
    pub fn mismatch_positions(&self) -> MdMismatchIter<'a> {
        MdMismatchIter {
            bytes: self.bytes,
            pos: 0,           // Current position in byte slice
            ref_offset: 0,    // Current position along reference
        }
    }

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

/// Iterator over mismatch positions in an MD tag.
///
/// Unlike the normalized version, this yields absolute offsets from read start,
/// along with the reference base that was mismatched.
pub struct MdMismatchIter<'a> {
    bytes: &'a [u8],
    pos: usize,
    ref_offset: usize,
}

impl<'a> Iterator for MdMismatchIter<'a> {
    // Yields (offset_from_read_start, reference_base) for each mismatch
    type Item = (usize, u8);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.bytes.len() {
            let c = self.bytes[self.pos];

            if c.is_ascii_digit() {
                // Parse match run - skip matching bases
                let mut num = 0usize;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                    num = num * 10 + (self.bytes[self.pos] - b'0') as usize;
                    self.pos += 1;
                }
                self.ref_offset += num;
            } else if c == b'^' {
                // Skip deletion bases
                self.pos += 1;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_uppercase() {
                    self.ref_offset += 1;
                    self.pos += 1;
                }
            } else if c.is_ascii_uppercase() {
                // Mismatch found - yield position and the reference base
                let offset = self.ref_offset;
                let base = c;  // The base is what the reference had
                self.ref_offset += 1;
                self.pos += 1;
                return Some((offset, base));
            } else {
                self.pos += 1;
            }
        }
        None
    }
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/// Helper to check if a read starts/ends with a match based on CIGAR and MD tag.
///
/// # Bioinformatics Context
/// For PhageTerm analysis, we need to know if a read truly aligns at its
/// start/end position. A read might report position 100, but if it starts
/// with soft clipping or a mismatch, position 100 isn't actually aligned.
///
/// This function checks BOTH:
/// 1. CIGAR says it starts with M/=/X (not S/H/I)
/// 2. MD tag confirms the first/last base is a match (not mismatch)
pub fn has_match_at_position(cigar: &Cigar, md: Option<&[u8]>, at_start: bool) -> bool {
    // Step 1: Check CIGAR - must not start with clip/insertion
    if !cigar.starts_with_match(at_start) {
        return false;
    }

    // Step 2: If MD tag exists, verify it also indicates a match
    // Rust Concept: Pattern matching with guards
    // `Some(bytes) if !bytes.is_empty()` only matches if Some AND condition is true
    match md {
        Some(bytes) if !bytes.is_empty() => {
            let md_tag = MdTag::new(bytes);
            if at_start {
                md_tag.has_match_at_start()
            } else {
                md_tag.has_match_at_end()
            }
        }
        Some(_) => false, // Empty MD tag - can't confirm match
        None => false,    // No MD tag - can't confirm match
    }
}

// ============================================================================
// RAW CIGAR HELPERS (ZERO-ALLOCATION)
// ============================================================================
// These functions work directly with rust-htslib's raw CIGAR format: &[(u32, u32)]
// where each tuple is (operation_char, length). This avoids allocating our
// CigarElement structs when we just need quick checks.
//
// # Performance Optimization
// The typed Cigar struct above is cleaner but requires allocation.
// These "raw" helpers operate on slices directly - zero allocation.
// In hot paths (called millions of times), this matters!

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
