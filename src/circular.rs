//! Circular array operations for genomic data.
//!
//! # Overview
//!
//! Many viral and bacterial genomes are circular, meaning reads can span the
//! origin and wrap around from the end of the sequence back to the beginning.
//! This module provides efficient operations for handling such wrap-around cases.
//!
//! For example, a genome of 10,000 bp with a read aligning from position 9,900
//! to position 100 requires special handling to correctly update coverage.
//!
//! # Python Equivalent
//!
//! In the original Python implementation (`calculating_data.py`), circular
//! handling was done inline:
//! ```python
//! # Python: calculate_coverage_numba() in calculating_data.py:277-287
//! if start_mod < end_mod:
//!     for j in range(start_mod, end_mod):
//!         coverage[j] += 1
//! else:
//!     # Handle wrap-around for circular contigs
//!     for j in range(start_mod, ref_length):
//!         coverage[j] += 1
//!     for j in range(0, end_mod):
//!         coverage[j] += 1
//! ```
//!
//! This module abstracts these operations into reusable traits and iterators.
//!
//! # Rust Concepts
//!
//! ## Traits
//! A trait defines shared behavior, similar to a Python abstract base class or
//! an interface. `CircularArray<T>` defines operations that any circular array
//! should support. The implementation `impl CircularArray<T> for Vec<T>` adds
//! these methods to Vec.
//!
//! ## Const Generics
//! `create_arrays<const N: usize>()` uses const generics - the array size N is
//! known at compile time, enabling efficient fixed-size array creation.

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
    ///
    /// # Python Equivalent
    /// ```python
    /// # From calculate_coverage_numba() in calculating_data.py:277-287
    /// if start_mod < end_mod:
    ///     for j in range(start_mod, end_mod):
    ///         coverage[j] += 1
    /// else:
    ///     for j in range(start_mod, ref_length):
    ///         coverage[j] += 1
    ///     for j in range(0, end_mod):
    ///         coverage[j] += 1
    /// ```
    #[inline]
    fn increment_circular(&mut self, start: usize, end: usize, delta: T) {
        if start <= end {
            // Normal case: [start, end)
            for pos in start..end {
                self[pos] += delta;
            }
        } else {
            // Wrap-around case: read spans the origin
            for pos in start..self.len() {
                self[pos] += delta;
            }
            for pos in 0..end {
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

    /// Create an iterator over positions in a circular range.
    fn circular_range(&self, start: usize, end: usize) -> CircularRangeIter {
        CircularRangeIter::new(start, end, self.len())
    }
}

// ============================================================================
// CIRCULAR RANGE ITERATOR
// ============================================================================

/// Iterator over positions in a circular range.
///
/// # Python Equivalent
/// ```python
/// positions = (np.arange(start, end) % ref_length)
/// ```
///
/// Unlike the numpy approach which allocates an array, this iterator yields
/// positions lazily without upfront memory allocation.
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
///
/// # Python Equivalent
/// ```python
/// pos = ref_starts[i] % ref_length
/// ```
#[inline]
pub fn normalize_position(pos: i64, length: usize) -> usize {
    (pos as usize) % length
}

/// Create multiple zero-initialized arrays of the same length.
///
/// # Python Equivalent
/// ```python
/// coverage = np.zeros(ref_length, dtype=np.uint64)
/// start_plus = np.zeros(ref_length, dtype=np.uint64)
/// # etc...
/// ```
#[inline]
pub fn create_arrays<const N: usize>(length: usize) -> [Vec<u64>; N] {
    std::array::from_fn(|_| vec![0u64; length])
}
