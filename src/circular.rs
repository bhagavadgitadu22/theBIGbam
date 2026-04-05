//! Circular array operations for genomic data.
//!
//! Handles circular genomes where reads can wrap around from the end
//! back to the beginning (e.g., position 9900-100 in a 10kb genome).

use std::ops::AddAssign;

// ============================================================================
// CIRCULAR ARRAY FUNCTIONS
// ============================================================================

/// Increment values in a circular range [start, end), handling wrap-around.
/// /// Start should already be modulo'd, but end can exceed len for wrapped reads.
#[inline]
pub fn increment_circular<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    let len = arr.len();
    let start_mod = start % len;
    let end_mod = end % len;

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

/// Increment values in a circular range for long reads that may span multiple reference lengths.
///
/// For reads longer than the reference, this adds full_laps to all positions first,
/// then handles the remaining partial coverage with increment_circular.
///
/// IMPORTANT: Pass raw_start and raw_end (before any modulo) to correctly calculate read length.
#[inline]
pub fn increment_circular_long<T>(arr: &mut [T], raw_start: usize, raw_end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    let len = arr.len();
    let read_length = raw_end - raw_start;

    // Calculate how many complete laps around the reference
    let full_laps = read_length / len;

    // Add full laps to all positions
    if full_laps > 0 {
        for pos in 0..len {
            for _ in 0..full_laps {
                arr[pos] += delta;
            }
        }
    }

    // Handle remaining partial coverage
    let remaining = read_length % len;
    if remaining > 0 {
        let start_mod = raw_start % len;
        let new_end = start_mod + remaining;
        increment_circular(arr, start_mod, new_end, delta);
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

// ============================================================================
// DECREMENT COUNTERPARTS (u64, saturating)
// ============================================================================
// These mirror the increment_* functions but use saturating_sub to guard
// against u64 underflow. Used to remove CIGAR 'N' (splice) spans from
// coverage-like arrays that were populated by a flat range fill.

/// Decrement values in a circular range [start, end), handling wrap-around.
#[inline]
pub fn decrement_circular(arr: &mut [u64], start: usize, end: usize, delta: u64) {
    let len = arr.len();
    let start_mod = start % len;
    let end_mod = end % len;

    if start_mod < end_mod {
        for pos in start_mod..end_mod {
            arr[pos] = arr[pos].saturating_sub(delta);
        }
    } else {
        for pos in start_mod..len {
            arr[pos] = arr[pos].saturating_sub(delta);
        }
        for pos in 0..end_mod {
            arr[pos] = arr[pos].saturating_sub(delta);
        }
    }
}

/// Decrement values in a circular range for long reads that may span multiple reference lengths.
#[inline]
pub fn decrement_circular_long(arr: &mut [u64], raw_start: usize, raw_end: usize, delta: u64) {
    let len = arr.len();
    let read_length = raw_end - raw_start;

    let full_laps = read_length / len;
    if full_laps > 0 {
        for pos in 0..len {
            for _ in 0..full_laps {
                arr[pos] = arr[pos].saturating_sub(delta);
            }
        }
    }

    let remaining = read_length % len;
    if remaining > 0 {
        let start_mod = raw_start % len;
        let new_end = start_mod + remaining;
        decrement_circular(arr, start_mod, new_end, delta);
    }
}

/// Decrement values in a non-circular range [start, end).
#[inline]
pub fn decrement_range(arr: &mut [u64], start: usize, end: usize, delta: u64) {
    for pos in start..end {
        arr[pos] = arr[pos].saturating_sub(delta);
    }
}
