//! Circular array operations for genomic data.
//!
//! Handles circular genomes where reads can wrap around from the end
//! back to the beginning (e.g., position 9900-100 in a 10kb genome).

use std::ops::AddAssign;

// ============================================================================
// CIRCULAR ARRAY FUNCTIONS
// ============================================================================

/// Increment values in a circular range [start, end), handling wrap-around.
/// Start should already be modulo'd, but end can exceed len for wrapped reads.
#[inline]
pub fn increment_circular<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    let len = arr.len();
    let start_mod = start%len;
    let end_mod = end%len;
    
    if start_mod < end_mod {
        // Normal case: no wrapping, end within bounds
        for pos in start_mod..end_mod {
            arr[pos] += delta;
        }
    } else {
        // Wrapped case: read spans the origin
        // Increment from start to end of array
        for pos in start_mod..len {
            arr[pos] += delta;
        }
        // Then from beginning to wrapped end position
        for pos in 0..end_mod {
            arr[pos] += delta;
        }
    }
}

/// Increment values in a non-circular range [start, end).
#[inline]
pub fn increment_range<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    for pos in start..end {
        arr[pos] += delta;
    }
}
