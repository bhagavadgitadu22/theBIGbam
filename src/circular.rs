//! Circular array operations for genomic data.
//!
//! Handles circular genomes where reads can wrap around from the end
//! back to the beginning (e.g., position 9900-100 in a 10kb genome).

use std::ops::AddAssign;

// ============================================================================
// CIRCULAR ARRAY TRAIT
// ============================================================================

/// Trait for arrays that support circular (wrap-around) operations.
///
/// Circular genomes require special handling where a range can span
/// from near the end of the sequence back to the beginning.
pub trait CircularArray<T> {
    /// Increment values in a circular range [start, end).
    ///
    /// If start <= end, increments positions [start, end).
    /// If start > end, increments [start, len) and [0, end) (wrap-around).
    fn increment_circular(&mut self, start: usize, end: usize, delta: T);

    /// Increment values in a circular range [start, end] (inclusive).
    fn increment_circular_inclusive(&mut self, start: usize, end: usize, delta: T);

    /// Increment values in a non-circular range [start, end).
    fn increment_range(&mut self, start: usize, end: usize, delta: T);

    /// Increment values in a non-circular range [start, end] (inclusive).
    fn increment_range_inclusive(&mut self, start: usize, end: usize, delta: T);

    /// Iterate over positions in a circular range [start, end).
    fn circular_range(&self, start: usize, end: usize) -> CircularRangeIter;
}

// ============================================================================
// IMPLEMENTATION FOR VEC<T>
// ============================================================================

impl<T> CircularArray<T> for Vec<T>
where
    T: AddAssign + Copy,
{
    /// Increment positions in range [start, end), handling wrap-around.
    /// Start should already be modulo'd, but end can exceed len for wrapped reads.
    #[inline]
    fn increment_circular(&mut self, start: usize, end: usize, delta: T) {
        let len = self.len();
        let end_mod = end % len;
        
        if end <= len {
            // Normal case: no wrapping, end within bounds
            for pos in start..end {
                self[pos] += delta;
            }
        } else if start < end_mod {
            // Wrapped: read spans origin (e.g., start=19800, end=39732 → end_mod=0)
            // This shouldn't happen with our logic, but handle it safely
            for pos in start..len {
                self[pos] += delta;
            }
            for pos in 0..end_mod {
                self[pos] += delta;
            }
        } else {
            // Read wraps around (start > end_mod)
            for pos in start..len {
                self[pos] += delta;
            }
            for pos in 0..end_mod {
                self[pos] += delta;
            }
        }
    }

    /// Increment positions in range [start, end] (inclusive), handling wrap-around.
    #[inline]
    fn increment_circular_inclusive(&mut self, start: usize, end: usize, delta: T) {
        let len = self.len();
        if start <= end {
            for pos in start..=end.min(len.saturating_sub(1)) {
                self[pos] += delta;
            }
        } else {
            for pos in start..len {
                self[pos] += delta;
            }
            for pos in 0..=end.min(len.saturating_sub(1)) {
                self[pos] += delta;
            }
        }
    }

    /// Increment positions in range [start, end) without wrap-around.
    #[inline]
    fn increment_range(&mut self, start: usize, end: usize, delta: T) {
        for pos in start..end {
            self[pos] += delta;
        }
    }

    /// Increment positions in range [start, end] (inclusive) without wrap-around.
    #[inline]
    fn increment_range_inclusive(&mut self, start: usize, end: usize, delta: T) {
        for pos in start..=end {
            self[pos] += delta;
        }
    }

    /// Create an iterator over positions in a circular range.
    fn circular_range(&self, start: usize, end: usize) -> CircularRangeIter {
        CircularRangeIter::new(start, end, self.len())
    }
}

// ============================================================================
// CIRCULAR RANGE ITERATOR
// ============================================================================

/// Iterator over positions in a circular range.
pub struct CircularRangeIter {
    current: usize,
    end: usize,
    len: usize,
    wrapped: bool,
    done: bool,
}

impl CircularRangeIter {
    fn new(start: usize, end: usize, len: usize) -> Self {
        Self {
            current: start,
            end,
            len,
            wrapped: start > end,
            done: false,
        }
    }
}

impl Iterator for CircularRangeIter {
    type Item = usize;

    /// Yield the next position in the circular range.
    ///
    /// For a wrapped range (e.g., start=9900, end=100, len=10000):
    /// - First yields 9900, 9901, ..., 9999
    /// - Then yields 0, 1, ..., 99
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        if self.wrapped {
            // Wrap-around case: first [start, len), then [0, end)
            if self.current < self.len {
                let pos = self.current;
                self.current += 1;
                Some(pos)
            } else if self.current == self.len {
                self.current = 0;
                if self.current < self.end {
                    let pos = self.current;
                    self.current += 1;
                    Some(pos)
                } else {
                    self.done = true;
                    None
                }
            } else if self.current < self.end {
                let pos = self.current;
                self.current += 1;
                Some(pos)
            } else {
                self.done = true;
                None
            }
        } else {
            // Simple case: [start, end)
            if self.current < self.end {
                let pos = self.current;
                self.current += 1;
                Some(pos)
            } else {
                self.done = true;
                None
            }
        }
    }
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/// Normalize a position to within array bounds using modulo.
#[inline]
pub fn normalize_position(pos: i64, length: usize) -> usize {
    (pos as usize) % length
}

/// Create multiple zero-initialized arrays of the same length.
#[inline]
pub fn create_arrays<const N: usize>(length: usize) -> [Vec<u64>; N] {
    std::array::from_fn(|_| vec![0u64; length])
}
