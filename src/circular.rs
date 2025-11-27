//! Circular array operations for genomic data.
//!
//! Provides efficient operations on arrays representing circular genomes,
//! where positions can wrap around from end to start.

use std::ops::AddAssign;

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

impl<T> CircularArray<T> for Vec<T>
where
    T: AddAssign + Copy,
{
    #[inline]
    fn increment_circular(&mut self, start: usize, end: usize, delta: T) {
        if start <= end {
            for pos in start..end {
                self[pos] += delta;
            }
        } else {
            // Wrap-around case
            for pos in start..self.len() {
                self[pos] += delta;
            }
            for pos in 0..end {
                self[pos] += delta;
            }
        }
    }

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

    fn circular_range(&self, start: usize, end: usize) -> CircularRangeIter {
        CircularRangeIter::new(start, end, self.len())
    }
}

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

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        if self.wrapped {
            // First iterate [start, len), then [0, end)
            if self.current < self.len {
                let pos = self.current;
                self.current += 1;
                Some(pos)
            } else if self.current == self.len {
                // Switch to second range
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

/// Normalize a position to within array bounds using modulo.
#[inline]
pub fn normalize_position(pos: i64, length: usize) -> usize {
    (pos as usize) % length
}

/// Create multiple zero-initialized arrays of the same length.
///
/// More efficient than multiple individual allocations.
#[inline]
pub fn create_arrays<const N: usize>(length: usize) -> [Vec<u64>; N] {
    std::array::from_fn(|_| vec![0u64; length])
}
