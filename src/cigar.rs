//! CIGAR and MD tag parsing for BAM alignment analysis.
//!
//! Provides typed representations and efficient parsing of CIGAR strings
//! and MD tags commonly found in BAM files.

/// CIGAR operation types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,          // M - alignment match (can be sequence match or mismatch)
    Insertion,      // I - insertion to the reference
    Deletion,       // D - deletion from the reference
    Skip,           // N - skipped region from the reference
    SoftClip,       // S - soft clipping (sequence present in SEQ)
    HardClip,       // H - hard clipping (sequence NOT present in SEQ)
    Padding,        // P - padding (silent deletion from padded reference)
    SeqMatch,       // = - sequence match
    SeqMismatch,    // X - sequence mismatch
}

impl CigarOp {
    /// Parse a CIGAR operation from its character representation.
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

    /// Parse from u32 (as stored internally from rust-htslib char()).
    #[inline]
    pub fn from_u32(op: u32) -> Option<Self> {
        Self::from_char(op as u8 as char)
    }

    /// Check if this operation is a clipping operation (soft or hard).
    #[inline]
    pub fn is_clipping(&self) -> bool {
        matches!(self, Self::SoftClip | Self::HardClip)
    }

    /// Check if this operation consumes the reference.
    #[inline]
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Deletion | Self::Skip | Self::SeqMatch | Self::SeqMismatch
        )
    }

    /// Check if this operation consumes the query.
    #[inline]
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Insertion | Self::SoftClip | Self::SeqMatch | Self::SeqMismatch
        )
    }
}

/// A single CIGAR operation with its length.
#[derive(Clone, Copy, Debug)]
pub struct CigarElement {
    pub op: CigarOp,
    pub len: u32,
}

impl CigarElement {
    /// Create from raw (op_char, length) tuple.
    #[inline]
    pub fn from_raw(op: u32, len: u32) -> Option<Self> {
        CigarOp::from_u32(op).map(|op| Self { op, len })
    }
}

/// A complete CIGAR string as a sequence of elements.
#[derive(Clone, Debug, Default)]
pub struct Cigar(pub Vec<CigarElement>);

impl Cigar {
    /// Create from raw tuples of (op_char, length).
    pub fn from_raw(raw: &[(u32, u32)]) -> Self {
        Self(
            raw.iter()
                .filter_map(|&(op, len)| CigarElement::from_raw(op, len))
                .collect(),
        )
    }

    /// Check if empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Get the first operation.
    #[inline]
    pub fn first(&self) -> Option<&CigarElement> {
        self.0.first()
    }

    /// Get the last operation.
    #[inline]
    pub fn last(&self) -> Option<&CigarElement> {
        self.0.last()
    }

    /// Check if the first operation is a clipping.
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
    /// `at_start`: if true, check the first op; if false, check the last op.
    pub fn starts_with_match(&self, at_start: bool) -> bool {
        let elem = if at_start { self.first() } else { self.last() };
        match elem {
            Some(e) => !matches!(e.op, CigarOp::SoftClip | CigarOp::HardClip | CigarOp::Insertion),
            None => false,
        }
    }

    /// Iterate over operations.
    pub fn iter(&self) -> impl Iterator<Item = &CigarElement> {
        self.0.iter()
    }
}

/// MD tag parser for extracting mismatch positions.
///
/// The MD tag format is: `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`
/// - Numbers indicate matching bases
/// - Single uppercase letters indicate mismatches
/// - ^[A-Z]+ indicates deletions
pub struct MdTag<'a> {
    bytes: &'a [u8],
}

impl<'a> MdTag<'a> {
    /// Create a new MD tag parser.
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { bytes }
    }

    /// Check if the tag indicates a match at the start.
    ///
    /// Returns true if the MD tag exists, is non-empty, and starts with a non-zero number.
    #[inline]
    pub fn has_match_at_start(&self) -> bool {
        !self.bytes.is_empty() && self.bytes[0].is_ascii_digit() && self.bytes[0] != b'0'
    }

    /// Check if the tag indicates a match at the end.
    #[inline]
    pub fn has_match_at_end(&self) -> bool {
        !self.bytes.is_empty() && self.bytes.last().map(|&b| b.is_ascii_digit() && b != b'0').unwrap_or(false)
    }

    /// Iterate over mismatch positions relative to reference start.
    ///
    /// Yields (reference_offset, mismatched_base) for each mismatch.
    pub fn mismatch_positions(&self) -> MdMismatchIter<'a> {
        MdMismatchIter {
            bytes: self.bytes,
            pos: 0,
            ref_offset: 0,
        }
    }

    /// Iterate over mismatch positions (normalized to ref_length).
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

/// Iterator over normalized mismatch positions in an MD tag.
pub struct MdMismatchNormalizedIter<'a> {
    bytes: &'a [u8],
    pos: usize,
    ref_pos: usize,
    ref_length: usize,
}

impl<'a> Iterator for MdMismatchNormalizedIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.bytes.len() {
            let c = self.bytes[self.pos];

            if c.is_ascii_digit() {
                // Parse number and advance ref_pos
                let mut num = 0usize;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                    num = num * 10 + (self.bytes[self.pos] - b'0') as usize;
                    self.pos += 1;
                }
                self.ref_pos += num;
            } else if c == b'^' {
                // Deletion: skip the deletion bases
                self.pos += 1;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_uppercase() {
                    self.ref_pos += 1;
                    self.pos += 1;
                }
            } else if c.is_ascii_uppercase() {
                // Mismatch found
                let result = self.ref_pos % self.ref_length;
                self.ref_pos += 1;
                self.pos += 1;
                return Some(result);
            } else {
                self.pos += 1;
            }
        }
        None
    }
}

/// Iterator over mismatch positions in an MD tag.
pub struct MdMismatchIter<'a> {
    bytes: &'a [u8],
    pos: usize,
    ref_offset: usize,
}

impl<'a> Iterator for MdMismatchIter<'a> {
    type Item = (usize, u8); // (reference_offset, mismatched_base)

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.bytes.len() {
            let c = self.bytes[self.pos];

            if c.is_ascii_digit() {
                // Parse number and advance ref_offset
                let mut num = 0usize;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                    num = num * 10 + (self.bytes[self.pos] - b'0') as usize;
                    self.pos += 1;
                }
                self.ref_offset += num;
            } else if c == b'^' {
                // Deletion: skip the deletion bases
                self.pos += 1;
                while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_uppercase() {
                    self.ref_offset += 1;
                    self.pos += 1;
                }
            } else if c.is_ascii_uppercase() {
                // Mismatch found
                let offset = self.ref_offset;
                let base = c;
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

/// Helper to check if a read starts/ends with a match based on CIGAR and MD tag.
pub fn has_match_at_position(cigar: &Cigar, md: Option<&[u8]>, at_start: bool) -> bool {
    // Check CIGAR first
    if !cigar.starts_with_match(at_start) {
        return false;
    }

    // If MD tag exists, verify it also indicates a match
    match md {
        Some(bytes) if !bytes.is_empty() => {
            let md_tag = MdTag::new(bytes);
            if at_start {
                md_tag.has_match_at_start()
            } else {
                md_tag.has_match_at_end()
            }
        }
        Some(_) => false, // Empty MD tag
        None => false,    // No MD tag
    }
}

// ============================================================================
// Raw CIGAR helpers (zero-allocation)
// ============================================================================

/// Check if raw CIGAR starts/ends with a match (no clipping/insertion).
#[inline]
pub fn raw_cigar_starts_with_match(cigar: &[(u32, u32)], at_start: bool) -> bool {
    let elem = if at_start { cigar.first() } else { cigar.last() };
    match elem {
        Some(&(op, _)) => {
            let c = op as u8 as char;
            // Not clipping (S/H) and not insertion (I)
            !matches!(c, 'S' | 'H' | 'I')
        }
        None => false,
    }
}

/// Check if first/last element of raw CIGAR is clipping.
#[inline]
pub fn raw_cigar_is_clipping(op: u32) -> bool {
    let c = op as u8 as char;
    matches!(c, 'S' | 'H')
}

/// Check if CIGAR op consumes reference.
#[inline]
pub fn raw_cigar_consumes_ref(op: u32) -> bool {
    let c = op as u8 as char;
    matches!(c, 'M' | 'D' | 'N' | '=' | 'X')
}

/// Check if a read starts/ends with match using raw CIGAR and MD tag (zero-allocation).
#[inline]
pub fn raw_has_match_at_position(cigar: &[(u32, u32)], md: Option<&[u8]>, at_start: bool) -> bool {
    if !raw_cigar_starts_with_match(cigar, at_start) {
        return false;
    }

    match md {
        Some(bytes) if !bytes.is_empty() => {
            let md_tag = MdTag::new(bytes);
            if at_start {
                md_tag.has_match_at_start()
            } else {
                md_tag.has_match_at_end()
            }
        }
        _ => false,
    }
}
