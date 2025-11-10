from pyexpat import features
import numpy as np
import pysam
from re import findall
from multiprocessing import Pool, cpu_count
import os
import csv
from numba import njit

import time
from functools import wraps
from collections import defaultdict

# global dictionary to accumulate times
function_times = defaultdict(float)

def track_time(func):
    """Decorator that records total time spent in each wrapped function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        function_times[func.__name__] += elapsed
        return result
    return wrapper

def print_timing_summary():
    print("\nTiming Summary", flush=True)
    total = sum(function_times.values())
    for name, t in sorted(function_times.items(), key=lambda x: -x[1]):
        print(f"{name:35s}: {t:.3f} s ({t/total*100:.1f}%)", flush=True)
    print(f"Total tracked time: {total:.3f} s", flush=True)

# Function dedicated to specific features
@njit
def calculate_coverage_numba(coverage, starts, ends, ref_length):
    for i in range(starts.size):
        start_mod = starts[i] % ref_length
        end_mod = ends[i] % ref_length

        if start_mod < end_mod:
            for j in range(start_mod, end_mod):
                coverage[j] += 1
        else:
            for j in range(start_mod, ref_length):
                coverage[j] += 1
            for j in range(0, end_mod):
                coverage[j] += 1

@track_time
def calculate_coverage(coverage, read, ref_length):
    blocks = np.array(read.get_blocks(), dtype=np.int32)
    if blocks.size == 0:
        return
    starts = blocks[:, 0]
    ends = blocks[:, 1]
    calculate_coverage_numba(coverage, starts, ends, ref_length)

def starts_with_match(read, start):
    """
    Return True if the first aligned base is a match (not clipped, not insertion, not mismatch).
    """
    if read.is_unmapped or read.cigartuples is None:
        return False

    # Check clipping at start
    first_op, _ = read.cigartuples[0]
    if not start:
        first_op, _ = read.cigartuples[-1]  # check end if start is False
    if first_op in (4, 5):  # soft or hard clip
        return False
    if first_op == 1:  # insertion
        return False

    # Now check MD tag
    try:
        md = read.get_tag("MD")
    except KeyError:
        raise ValueError("No MD tag present in BAM")

    md_status = md[0].isdigit()
    if not start:
        md_status = md[-1].isdigit()  # check end if start is False
    return md_status

@track_time
def calculate_reads_starts_and_ends(read, temporary_dict, ref_length):
    if starts_with_match(read, True) and starts_with_match(read, False):
        start = read.reference_start % ref_length
        end = (read.reference_end - 1) % ref_length

        s = read.reference_start
        e = read.reference_end
        if s < ref_length and e > ref_length:
            temporary_dict["coverage_reduced"][s:ref_length] += 1
            temporary_dict["coverage_reduced"][0:e % ref_length] += 1
        else:
            temporary_dict["coverage_reduced"][s % ref_length:e % ref_length] += 1

        if read.is_reverse:
            temporary_dict["start_minus"][start] += 1
            temporary_dict["end_minus"][end] += 1
        else:
            temporary_dict["start_plus"][start] += 1
            temporary_dict["end_plus"][end] += 1

@track_time
def compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length):
    feature_dict["coverage_reduced"] = temporary_dict["coverage_reduced"]
    if sequencing_type == "short":
        feature_dict["reads_starts"] = temporary_dict["start_plus"]
        feature_dict["reads_ends"] = temporary_dict["end_minus"]
    else:
        feature_dict["reads_starts"] = temporary_dict["start_plus"] + temporary_dict["start_minus"]
        feature_dict["reads_ends"] = temporary_dict["end_plus"] + temporary_dict["end_minus"]

    tau = np.zeros(ref_length, dtype=np.float32)
    mask = feature_dict["coverage_reduced"] > 0
    tau[mask] = (feature_dict["reads_starts"][mask] + feature_dict["reads_ends"][mask]) / feature_dict["coverage_reduced"][mask]
    feature_dict["tau"] = tau

@track_time
def compute_read_lengths(sum_read_lengths, counts_read_lengths, read, ref_length):
    positions = np.fromiter(read.get_reference_positions(), dtype=np.int32)
    if positions.size == 0:
        return
    positions %= ref_length
    np.add.at(sum_read_lengths, positions, read.query_length)
    np.add.at(counts_read_lengths, positions, 1)

@njit
def add_insert_sizes_numba(sum_insert_sizes, counts_insert_sizes, bad_orientations, positions, insert_len, is_read1, proper_pair, insert_flag, bad_flag, ref_length):
    for i in range(positions.size):
        pos = positions[i] % ref_length
        if insert_flag and is_read1 and insert_len > 0:
            sum_insert_sizes[pos] += insert_len
            counts_insert_sizes[pos] += 1
        if bad_flag and not proper_pair:
            bad_orientations[pos] += 1

@track_time
def compute_inserts_characteristics(sum_insert_sizes, counts_insert_sizes, bad_orientations, features, read, ref_length):
    if not read.is_paired or read.mate_is_unmapped:
        return
    positions = np.fromiter(read.get_reference_positions(), dtype=np.int32)
    if positions.size == 0:
        return
    insert_len = abs(read.template_length)
    insert_flag = 1 if "insert_sizes" in features else 0
    bad_flag = 1 if "bad_orientations" in features else 0
    add_insert_sizes_numba(sum_insert_sizes, counts_insert_sizes, bad_orientations, positions, insert_len, read.is_read1, read.is_proper_pair, insert_flag, bad_flag, ref_length)

@track_time
def compute_clippings(feature_dict, features, read, ref_length):
    if read.cigartuples:
        first_op, _ = read.cigartuples[0]
        last_op, _ = read.cigartuples[-1]
        if "left_clippings" in features and first_op in (4, 5):
            feature_dict["left_clippings"][read.reference_start % ref_length] += 1
        if "right_clippings" in features and last_op in (4, 5):
            feature_dict["right_clippings"][(read.reference_end - 1) % ref_length] += 1

@njit
def add_indels_numba(deletions, insertions, cigartuples, ref_start, ref_length):
    ref_pos = ref_start
    for i in range(cigartuples.shape[0]):
        op = cigartuples[i, 0]
        length = cigartuples[i, 1]
        if op == 1 and insertions is not None:  # insertion
            pos = ref_pos % ref_length
            insertions[pos] += 1
        elif op == 2 and deletions is not None:  # deletion
            for j in range(length):
                deletions[(ref_pos + j) % ref_length] += 1
            ref_pos += length
        else:
            ref_pos += length

@track_time
def compute_indels(feature_dict, features, read, ref_length):
    if read.cigartuples is None:
        return feature_dict
    cigartuples = np.array(read.cigartuples, dtype=np.int32)
    deletions = feature_dict["deletions"] if "deletions" in features else None
    insertions = feature_dict["insertions"] if "insertions" in features else None
    add_indels_numba(deletions, insertions, cigartuples, read.reference_start, ref_length)
    return feature_dict

@njit
def add_mismatches_numba_from_md(mismatches, md_chars, md_len, ref_start, ref_length):
    ref_pos = ref_start
    i = 0
    while i < md_len:
        c = md_chars[i]
        # --- Number: consecutive matches ---
        if 48 <= c <= 57:  # '0'-'9'
            num = 0
            while i < md_len and 48 <= md_chars[i] <= 57:
                num = num * 10 + (md_chars[i] - 48)
                i += 1
            ref_pos += num
        # --- Deletion from reference '^' --- 
        elif c == 94:  # '^'
            i += 1
            while i < md_len and (65 <= md_chars[i] <= 90):  # 'A'-'Z'
                ref_pos += 1
                i += 1
        # --- Mismatch --- 
        else:  # 'A'-'Z'
            mismatches[ref_pos % ref_length] += 1
            ref_pos += 1
            i += 1

@track_time
def compute_mismatches(feature_dict, features, read, ref_length):
    if "mismatches" not in features or not read.has_tag("MD"):
        return
    md_tag = read.get_tag("MD")
    # Convert string to numpy array of bytes for Numba
    md_chars = np.frombuffer(md_tag.encode('ascii'), dtype=np.uint8)
    add_mismatches_numba_from_md(feature_dict["mismatches"], md_chars, md_chars.size, read.reference_start, ref_length)
    
@track_time
def compute_final_lengths(sum_lengths, count_lengths, ref_length):
    count_lengths = np.maximum(count_lengths, 1)  # avoid division by zero
    arr = sum_lengths / count_lengths
    arr[count_lengths == 0] = 0
    return arr.astype(np.float32)

### Main logic to get features per position
@track_time
def get_features(reads_mapped, features, ref_length, sequencing_type):    
    # Initialize arrays
    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {}
    
    # Temporary arrays
    if "tau" in features:
        for f in ["start_plus", "start_minus", "end_plus", "end_minus", "coverage_reduced"]:
            temporary_dict[f] = np.zeros(ref_length, dtype=np.uint64)
    if "read_lengths" in features:
        temporary_dict["sum_read_lengths"] = np.zeros(ref_length, dtype=np.uint64)
        temporary_dict["count_read_lengths"] = np.zeros(ref_length, dtype=np.uint64)
    if "insert_sizes" in features:
        temporary_dict["sum_insert_sizes"] = np.zeros(ref_length, dtype=np.uint64)
        temporary_dict["count_insert_sizes"] = np.zeros(ref_length, dtype=np.uint64)
    if "bad_orientations" in features:
        feature_dict["bad_orientations"] = np.zeros(ref_length, dtype=np.float32)

    for read in reads_mapped:
        if read.is_unmapped:
            continue

        if "coverage" in features:
            calculate_coverage(feature_dict["coverage"], read, ref_length)
        if "tau" in features:
            calculate_reads_starts_and_ends(read, temporary_dict, ref_length)
        if "read_lengths" in features:
            compute_read_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"], read, ref_length)
        if {"insert_sizes", "bad_orientations"}.intersection(features):
            compute_inserts_characteristics(
                temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"],
                feature_dict.get("bad_orientations", None),
                features, read, ref_length
            )
        if {"left_clippings", "right_clippings"}.intersection(features):
            compute_clippings(feature_dict, features, read, ref_length)
        if {"insertions", "deletions"}.intersection(features):
            compute_indels(feature_dict, features, read, ref_length)
        if "mismatches" in features:
            compute_mismatches(feature_dict, features, read, ref_length)

    if "tau" in features:
        compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length)
    if "read_lengths" in features:
        feature_dict["read_lengths"] = compute_final_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"], ref_length)
    if "insert_sizes" in features:
        feature_dict["insert_sizes"] = compute_final_lengths(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"], ref_length)

    return feature_dict

### Summarise data for each feature
def merge_sorted_unique(*arrays):
    """Merge multiple sorted 1D arrays into a unique sorted array (very fast)."""
    arrays = [arr for arr in arrays if len(arr) > 0]
    if not arrays:
        return np.array([], dtype=int)
    merged = np.concatenate(arrays)
    merged.sort(kind='mergesort')  # stable and fast for pre-sorted subarrays
    # Drop duplicates efficiently
    return merged[np.concatenate(([True], np.diff(merged) != 0))]

@track_time
def compress_signal(feature_values, ref_length, step=50, z_thresh=3, deriv_thresh=3, max_points=10000):
    """
    Fast signal compression while preserving spikes and sharp transitions.
    
    Parameters
    ----------
    x, y : np.ndarray, Input signal data (same length).
    step : int, Keep every Nth point.
    z_thresh : float, Z-score threshold for keeping high/low outlier values.
    deriv_thresh : float, Z-score threshold for keeping large derivative changes.
    max_points : int, Hard limit on total number of kept points.
    
    Returns
    -------
    x_new, y_new : np.ndarray, Compressed signal retaining key features.
    """
    # Value outliers
    y_mean = np.mean(feature_values)
    y_std = np.std(feature_values) or 1e-9
    val_outliers = np.abs(feature_values - y_mean) > z_thresh * y_std

    # Derivative outliers
    dy = np.diff(feature_values, prepend=feature_values[0])
    dy_std = np.std(dy) or 1e-9
    der_outliers = np.abs(dy) > deriv_thresh * dy_std

    # Regular subsampling
    n = len(feature_values)
    keep_idx = merge_sorted_unique(
        np.arange(0, n, step, dtype=int),
        np.nonzero(val_outliers | der_outliers)[0],
        np.array([n - 1], dtype=int)
    )

    # Apply hard limit if needed
    if len(keep_idx) > max_points:
        step_lim = len(keep_idx) // max_points
        keep_idx = keep_idx[::step_lim]
        if keep_idx[-1] != n - 1:
            keep_idx = np.append(keep_idx, n - 1)

    # Initialize x array of ref_length size
    x = np.arange(1, ref_length + 1)
    return {"x": x[keep_idx], "y": feature_values[keep_idx]}

### Calculating features per contig per sample
@track_time
def calculating_features_per_contig_per_sample(reads_mapped, locus_size, feature_list, sequencing_type, window_size):
    sequencing_type = "short" if sequencing_type.startswith("short") else "long"

    # Read bam once to calculate all features
    feature_values = get_features(reads_mapped, feature_list, locus_size, sequencing_type)

    # Averaging and masking data per feature
    data = {}
    for feature in feature_list:
        feature_compressed = compress_signal(feature_values[feature], locus_size)
        if feature_compressed is not None:
            data[feature] = feature_compressed
                
    return data
    
### Distributing the calculation for all the contigs in the bam file
def calculating_features_per_sample(mapping_file, feature_list, sequencing_type, window_size):
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    references = bam_file.references
    lengths = [l // 2 for l in bam_file.lengths]  # need to divide by 2 because each contig was doubled for mapping to deal with circularity

    all_data = {}
    for ref, length in zip(references, lengths):
        # Fetch reads only for this contig
        reads_iter = bam_file.fetch(ref)
        contig_data = calculating_features_per_contig_per_sample(reads_iter, length, feature_list, sequencing_type, window_size)
        all_data[ref] = contig_data

    bam_file.close()
    print_timing_summary()
    return all_data

### Save data computed per sample
# Will be replaced by calls to database in future version
def save_data_dictionary(bam_name, contigs_list, output_file):
    """
    Save a nested data dictionary to a CSV file.
    """
    # Define CSV headers
    fieldnames = ["sample", "contig", "variable", "position", "value"]

    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for contig_name, variables_list in contigs_list.items():
            for variable_name, content in variables_list.items():
                xs = content.get("x", [])
                ys = content.get("y", [])
                # Ensure lengths match
                for pos, val in zip(xs, ys):
                    writer.writerow({
                        "sample": bam_name,
                        "contig": contig_name,
                        "variable": variable_name,
                        "position": pos,
                        "value": val
                    })

### Distributing the calculation for all the bam files (samples)
def _process_single_sample(bam_file, feature_list, sequencing_type, window_size, output_prefix, n_contig_cores=None):
    sample_name = os.path.basename(bam_file).replace(".bam", "")

    # Calculate features per contig (can still be parallelized with n_contig_cores)
    sample_data = calculating_features_per_sample(bam_file, feature_list, sequencing_type, window_size)

    # Save to CSV
    output_name = f"{output_prefix}_{os.path.splitext(sample_name)[0]}_features.csv"
    save_data_dictionary(sample_name, sample_data, output_name)
    return sample_name, sample_data

def calculating_all_features_parallel(bam_files, feature_list, sequencing_type, window_size, output_prefix, n_sample_cores=None, n_contig_cores=None):
    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)
    print(f"Using {n_sample_cores} cores to process {len(bam_files)} samples in parallel...", flush=True)

    args_list = [(bam, feature_list, sequencing_type, window_size, output_prefix, n_contig_cores) for bam in bam_files]

    with Pool(processes=n_sample_cores) as pool:
        results = pool.starmap(_process_single_sample, args_list)

    all_samples_data = {sample: data for sample, data in results}
    print("Finished all samples.", flush=True)

    return all_samples_data