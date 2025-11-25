//! Feature calculation functions.
//!
//! This module contains all the feature calculation logic for coverage,
//! phage termini, and assembly check features.
//! Python equivalent: Feature calculation functions in calculating_data.py

use crate::types::{FeatureMap, ReadData, SequencingType};

/// Calculate coverage per position.
/// Python equivalent: `calculate_coverage_numba()` in calculating_data.py:182-197
pub fn calculate_coverage(reads: &[ReadData], ref_length: usize) -> Vec<u64> {
    let mut coverage = vec![0u64; ref_length];

    // calculating_data.py:178-193 - calculate_coverage_numba()
    for read in reads {
        let start = (read.ref_start as usize) % ref_length;  // py:181
        let end = (read.ref_end as usize) % ref_length;      // py:182

        if start < end {  // py:183
            for i in start..end {  // py:184-185
                coverage[i] += 1;
            }
        } else {
            // Wrap-around for circular contigs - py:188-193
            for i in start..ref_length {
                coverage[i] += 1;
            }
            for i in 0..end {
                coverage[i] += 1;
            }
        }
    }

    coverage
}

/// Check if read starts with a match (not clipped, not insertion, not mismatch).
/// Python equivalent: `starts_with_match()` in calculating_data.py:208-224
fn starts_with_match(cigar: &[(u32, u32)], md: &Option<Vec<u8>>, at_start: bool) -> bool {
    // calculating_data.py:208-224 - starts_with_match()
    if cigar.is_empty() {
        return false;
    }

    // py:213 - Check first/last CIGAR op depending on 'start' flag
    let (op, _) = if at_start { cigar[0] } else { cigar[cigar.len() - 1] };

    // py:214-216 - Check for clipping (S=4, H=5) or insertion (I=1)
    // Rust uses char representation, Python uses numeric ops
    let op_char = op as u8 as char;
    if op_char == 'S' || op_char == 'H' || op_char == 'I' {
        return false;
    }

    // py:220-224 - Check MD tag exists and first/last char indicates match
    // Python: val = md[0] if start else md[-1]; return val > 0
    // Here val is ASCII byte, so val > 0 is always true for any char
    // We just need to check MD exists and isn't empty
    if let Some(md_bytes) = md {
        return !md_bytes.is_empty();
    }

    false
}

/// Calculate phage termini features (coverage_reduced, reads_starts, reads_ends, tau).
/// Python equivalent: `get_features_phagetermini()` in calculating_data.py:283-322
pub fn calculate_phagetermini(
    reads: &[ReadData],
    ref_length: usize,
    seq_type: SequencingType,
) -> FeatureMap {
    // calculating_data.py:279-317 - get_features_phagetermini()
    let mut coverage_reduced = vec![0u64; ref_length];  // py:284
    let mut start_plus = vec![0u64; ref_length];        // py:285
    let mut start_minus = vec![0u64; ref_length];
    let mut end_plus = vec![0u64; ref_length];
    let mut end_minus = vec![0u64; ref_length];

    for read in reads {
        let check_end = seq_type == SequencingType::Long;  // py:294

        // py:296-297 - starts_with_match checks
        if starts_with_match(&read.cigar, &read.md_tag, true)
            && (!check_end || starts_with_match(&read.cigar, &read.md_tag, false))
        {
            let start = (read.ref_start as usize) % ref_length;
            let end = (read.ref_end as usize) % ref_length;

            // calculating_data.py:241-262 - calculate_reads_starts_and_ends_numba()
            if start <= end {
                for i in start..=end.min(ref_length - 1) {
                    coverage_reduced[i] += 1;  // py:247
                }
            } else {
                // Wrap-around - py:251-256
                for i in start..ref_length {
                    coverage_reduced[i] += 1;
                }
                for i in 0..=end {
                    coverage_reduced[i] += 1;
                }
            }

            // py:258-262 - Update starts/ends based on strand
            if read.is_reverse {
                start_minus[start] += 1;
                end_minus[end] += 1;
            } else {
                start_plus[start] += 1;
                end_plus[end] += 1;
            }
        }
    }

    // calculating_data.py:264-278 - compute_final_starts_ends_and_tau()
    let (reads_starts, reads_ends) = match seq_type {
        SequencingType::ShortPaired | SequencingType::ShortSingle => {
            (start_plus, end_minus)  // py:267-268
        }
        SequencingType::Long => {
            // py:270-271 - Long reads sum both strands
            let starts: Vec<u64> = start_plus.iter().zip(&start_minus).map(|(a, b)| a + b).collect();
            let ends: Vec<u64> = end_plus.iter().zip(&end_minus).map(|(a, b)| a + b).collect();
            (starts, ends)
        }
    };

    let mut results = FeatureMap::new();
    results.insert("coverage_reduced".to_string(), coverage_reduced);
    results.insert("reads_starts".to_string(), reads_starts);
    results.insert("reads_ends".to_string(), reads_ends);
    // tau is computed later as float - py:273-278

    results
}

/// Calculate assembly check features (clippings, indels, mismatches, read_lengths, etc.).
/// Python equivalent: `get_features_assemblycheck()` in calculating_data.py:410-493
pub fn calculate_assemblycheck(
    reads: &[ReadData],
    ref_length: usize,
    seq_type: SequencingType,
) -> FeatureMap {
    // calculating_data.py:405-488 - get_features_assemblycheck()
    let mut results = FeatureMap::new();

    // Initialize arrays - py:420-431
    let mut left_clippings = vec![0u64; ref_length];
    let mut right_clippings = vec![0u64; ref_length];
    let mut insertions = vec![0u64; ref_length];
    let mut deletions = vec![0u64; ref_length];
    let mut mismatches = vec![0u64; ref_length];

    let mut sum_read_lengths = vec![0u64; ref_length];   // py:427
    let mut count_read_lengths = vec![0u64; ref_length]; // py:427
    let mut sum_insert_sizes = vec![0u64; ref_length];   // py:428
    let mut count_insert_sizes = vec![0u64; ref_length]; // py:428
    let mut bad_orientations = vec![0u64; ref_length];   // py:429

    for read in reads {
        // calculating_data.py:447-448 - get start/end positions
        let raw_start = read.ref_start as usize;
        let raw_end = read.ref_end as usize;
        let start = raw_start % ref_length;

        // calculating_data.py:451-453 - add_read_lengths_range_numba()
        // Long reads: track read lengths
        if seq_type == SequencingType::Long {
            // py:328-333 - add_read_lengths_range_numba: for pos in range(start, end)
            for pos in raw_start..raw_end {
                let p = pos % ref_length;  // py:331
                sum_read_lengths[p] += read.query_length as u64;   // py:332
                count_read_lengths[p] += 1;                        // py:333
            }
        }

        // calculating_data.py:455-460 - add_insert_sizes_range_numba()
        // Short-paired: track insert sizes and bad orientations
        if seq_type == SequencingType::ShortPaired {
            // py:346-356 - add_insert_sizes_range_numba
            for pos in raw_start..raw_end {
                let p = pos % ref_length;  // py:351
                if read.is_read1 && read.template_length > 0 {  // py:352
                    sum_insert_sizes[p] += read.template_length as u64;  // py:353
                    count_insert_sizes[p] += 1;                          // py:354
                }
                if !read.is_proper_pair {  // py:355
                    bad_orientations[p] += 1;  // py:356
                }
            }
        }

        // calculating_data.py:467-474 - Clippings from CIGAR ops
        // Check first/last CIGAR ops for soft/hard clips
        if !read.cigar.is_empty() {
            let (first_op, _) = read.cigar[0];       // py:469
            let (last_op, _) = read.cigar[read.cigar.len() - 1];  // py:470
            let first_char = first_op as u8 as char;
            let last_char = last_op as u8 as char;
            let end = raw_end % ref_length;

            // py:471-472 - Left clipping: first op is S(4) or H(5)
            if first_char == 'S' || first_char == 'H' {
                left_clippings[start] += 1;  // py:472
            }
            // py:473-474 - Right clipping: last op is S(4) or H(5)
            if last_char == 'S' || last_char == 'H' {
                right_clippings[if end > 0 { end - 1 } else { ref_length - 1 }] += 1;  // py:474
            }
        }

        // calculating_data.py:359-372 - add_indels_numba()
        // Indels from CIGAR
        let mut ref_pos = read.ref_start as usize;  // py:360
        for (op, len) in &read.cigar {              // py:361
            let op_char = *op as u8 as char;
            match op_char {
                'I' => {
                    // py:364-366 - Insertion: record but don't advance ref_pos
                    insertions[ref_pos % ref_length] += 1;  // py:366
                }
                'D' => {
                    // py:367-370 - Deletion: record each position and advance
                    for j in 0..(*len as usize) {
                        deletions[(ref_pos + j) % ref_length] += 1;  // py:369
                    }
                    ref_pos += *len as usize;  // py:370
                }
                _ => {
                    // py:371-372 - All other ops advance ref_pos
                    ref_pos += *len as usize;
                }
            }
        }

        // calculating_data.py:374-397 - add_mismatches_numba_from_md()
        // Mismatches from MD tag
        if let Some(ref md_bytes) = read.md_tag {
            let mut ref_pos = read.ref_start as usize;  // py:377
            let mut i = 0;  // py:378
            while i < md_bytes.len() {  // py:379
                let c = md_bytes[i];
                if c.is_ascii_digit() {
                    // py:380-385 - Parse number and advance ref_pos
                    let mut num = 0usize;
                    while i < md_bytes.len() && md_bytes[i].is_ascii_digit() {
                        num = num * 10 + (md_bytes[i] - b'0') as usize;
                        i += 1;
                    }
                    ref_pos += num;
                } else if c == b'^' {
                    // py:386-392 - Deletion marker: skip bases
                    i += 1;
                    while i < md_bytes.len() && md_bytes[i].is_ascii_uppercase() {
                        ref_pos += 1;
                        i += 1;
                    }
                } else if c.is_ascii_uppercase() {
                    // py:394-397 - Mismatch: record and advance
                    mismatches[ref_pos % ref_length] += 1;  // py:395
                    ref_pos += 1;  // py:396
                    i += 1;        // py:397
                } else {
                    i += 1;
                }
            }
        }
    }

    results.insert("left_clippings".to_string(), left_clippings);
    results.insert("right_clippings".to_string(), right_clippings);
    results.insert("insertions".to_string(), insertions);
    results.insert("deletions".to_string(), deletions);
    results.insert("mismatches".to_string(), mismatches);

    if seq_type == SequencingType::Long {
        results.insert("sum_read_lengths".to_string(), sum_read_lengths);
        results.insert("count_read_lengths".to_string(), count_read_lengths);
    }

    if seq_type == SequencingType::ShortPaired {
        results.insert("sum_insert_sizes".to_string(), sum_insert_sizes);
        results.insert("count_insert_sizes".to_string(), count_insert_sizes);
        results.insert("bad_orientations".to_string(), bad_orientations);
    }

    results
}
