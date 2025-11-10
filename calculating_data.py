import argparse, sys, os, csv
from multiprocessing import Pool, cpu_count
from pyexpat import features
from numba import njit
import numpy as np
import pysam

import time
from functools import wraps
from collections import defaultdict

### Measuring time spent in each function
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

### Find sequencing type from BAM file
def find_sequencing_type_from_bam(bam_file, n_reads_check=100):
    """
    Infer sequencing type from a BAM file.
    Rules:
      - If any read > 1000 bp, reads are "long"
      - Elif any read.is_paired, reads are "short-paired"
      - Else reads are "short-single"
    """
    bam = pysam.AlignmentFile(bam_file, "rb")

    for i, read in enumerate(bam.fetch(until_eof=True)):
        if read.is_unmapped:
            continue
        if read.query_length and read.query_length > 1000:
            bam.close()
            return "long"
        if read.is_paired:
            return "short-paired"

        if i + 1 >= n_reads_check:
            break

    bam.close()
    return "short-single"

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

def compress_signal(feature_values, ref_length, step=50, z_thresh=3, deriv_thresh=3, max_points=10000):
    """
    step : int, Keep every Nth point
    z_thresh : float, Z-score threshold for keeping high/low outlier values
    deriv_thresh : float, Z-score threshold for keeping large derivative changes
    max_points : int, Hard limit on total number of kept points
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

### Save data computed per sample
# Will be replaced by calls to database in future version
def save_feature_values_per_contig_per_sample(feature, feature_values, output_dir, sample_name, ref_name):
    # Define CSV headers
    fieldnames = ["sample", "contig", "position", "value"]

    output_file = os.path.join(output_dir, feature+"_values_for_"+ref_name+"_in_"+sample_name+".csv")
    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        xs = feature_values.get("x", [])
        ys = feature_values.get("y", [])
        # Ensure lengths match
        for pos, val in zip(xs, ys):
            writer.writerow({
                "sample": sample_name,
                "contig": ref_name,
                "position": pos,
                "value": val
            })

### Functions of coverage module
@njit
def calculate_coverage_numba(coverage, starts, ends, ref_length):
    n = starts.size
    for i in range(n):
        start_mod = starts[i] % ref_length
        end_mod = ends[i] % ref_length

        if start_mod < end_mod:
            for j in range(start_mod, end_mod):
                coverage[j] += 1
        else:
            # Handle wrap-around for circular contigs
            for j in range(start_mod, ref_length):
                coverage[j] += 1
            for j in range(0, end_mod):
                coverage[j] += 1

def get_feature_coverage(ref_starts, ref_ends, ref_length, output_dir, sample_name, ref_name):
    coverage = np.zeros(ref_length, dtype=np.uint64)
    calculate_coverage_numba(coverage, ref_starts, ref_ends, ref_length)

    coverage_compact = compress_signal(coverage, ref_length)
    save_feature_values_per_contig_per_sample("coverage", coverage_compact, output_dir, sample_name, ref_name)

### Functions of phagetermini module
def starts_with_match(cigar, md, start):
    """
    Return True if the first aligned base is a match (not clipped, not insertion, not mismatch).
    """
    # Check clipping or insertion at start/end
    op, length = cigar[0] if start else cigar[-1]
    if op in (4, 5):  # soft or hard clip
        return False
    if op == 1:  # insertion
        return False

    # Check MD tag for match at start/end
    # MD string: digits represent matches, letters/deletions mismatches
    val = md[0] if start else md[-1]
    return val > 0

def calculate_reads_starts_and_ends(start, end, is_reverse, temporary_dict, ref_length):
    start = start % ref_length
    end = end % ref_length

    # Update coverage
    if start <= end:
        temporary_dict["coverage_reduced"][start:end+1] += 1
    else:
        temporary_dict["coverage_reduced"][start:ref_length] += 1
        temporary_dict["coverage_reduced"][0:end+1] += 1

    # Update strand-specific starts/ends
    if is_reverse:
        temporary_dict["start_minus"][start] += 1
        temporary_dict["end_minus"][end] += 1
    else:
        temporary_dict["start_plus"][start] += 1
        temporary_dict["end_plus"][end] += 1

def compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length):
    feature_dict["coverage_reduced"] = temporary_dict["coverage_reduced"]
    if sequencing_type == "short-paired" or sequencing_type == "short-single":
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
def get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, 
                              sequencing_type, ref_length, output_dir, sample_name, ref_name):
    features = ["coverage_reduced", "reads_starts", "reads_ends", "tau"]
    temporary_features = ["coverage_reduced", "start_plus", "start_minus", "end_plus", "end_minus"]

    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    # Iterate through preprocessed reads
    for i in range(len(ref_starts)):
        cigar = cigars[i]
        md = md_list[i]
        if starts_with_match(cigar, md, start=True) and starts_with_match(cigar, md, start=False):
            start = ref_starts[i]
            end = ref_ends[i]
            is_reverse_flag = is_reverse[i]
            calculate_reads_starts_and_ends(start, end, is_reverse_flag, temporary_dict, ref_length)

    # Compute tau
    compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length)

    # Save features
    for feature in features:
        feature_values = compress_signal(feature_dict[feature], ref_length)
        save_feature_values_per_contig_per_sample(feature, feature_values, output_dir, sample_name, ref_name)

### Functions of assemblycheck module
@njit
def add_read_lengths_numba(sum_read_lengths, counts_read_lengths, positions, read_len, ref_length):
    for i in range(positions.size):
        pos = positions[i] % ref_length
        sum_read_lengths[pos] += read_len
        counts_read_lengths[pos] += 1

@track_time
def compute_read_lengths(sum_read_lengths, counts_read_lengths, read, ref_length):
    positions = np.fromiter(read.get_reference_positions(), dtype=np.int32)
    if positions.size == 0:
        return
    add_read_lengths_numba(sum_read_lengths, counts_read_lengths, positions, read.query_length, ref_length)

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
def compute_final_lengths(sum_lengths, count_lengths):
    count_lengths = np.maximum(count_lengths, 1)  # avoid division by zero
    arr = sum_lengths / count_lengths
    arr[count_lengths == 0] = 0
    return arr.astype(np.float32)

@track_time
def get_features_assemblycheck(reads_mapped, ref_length, sequencing_type, output_dir, sample_name, ref_name):    
    features = []
    temporary_features = []

    # Select features depending on sequencing type
    if sequencing_type == "long":
        features.extend(["read_lengths"])
        temporary_features.extend(["sum_read_lengths", "count_read_lengths"])
    if sequencing_type == "short-paired":
        features.extend(["insert_sizes", "bad_orientations"])
        temporary_features.extend(["sum_insert_sizes", "count_insert_sizes", "bad_orientation"])
    features.extend(["left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])

    # Initialize arrays
    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    for read in reads_mapped:
        if sequencing_type == "long":
            compute_read_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"], read, ref_length)
        if sequencing_type == "short-paired":
            compute_inserts_characteristics(
                temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"],
                feature_dict.get("bad_orientations", None),
                features, read, ref_length
            )
        compute_clippings(feature_dict, features, read, ref_length)
        compute_indels(feature_dict, features, read, ref_length)
        compute_mismatches(feature_dict, features, read, ref_length)

    if sequencing_type == "long":
        feature_dict["read_lengths"] = compute_final_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"])
    if sequencing_type == "short-paired":
        feature_dict["insert_sizes"] = compute_final_lengths(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"])

    for feature in features:
        feature_values = compress_signal(feature_dict[feature], ref_length)
        save_feature_values_per_contig_per_sample(feature, feature_values, output_dir, sample_name, ref_name)

### Calculating features per contig per sample
@track_time
def preprocess_reads(reads_mapped, modules):
    ref_starts, ref_ends = [], []
    query_lengths, template_lengths = [], []
    is_read1, is_proper_pair, is_paired, is_reverse = [], [], [], []
    cigars, has_md, md_list, md_lengths = [], [], [], []

    for read in reads_mapped:
        if read.is_unmapped:
            continue

        ref_starts.append(read.reference_start)
        ref_ends.append(read.reference_end)
        query_lengths.append(read.query_length)
        template_lengths.append(abs(read.template_length))
        is_read1.append(read.is_read1)
        is_proper_pair.append(read.is_proper_pair)
        is_paired.append(read.is_paired)
        is_reverse.append(read.is_reverse)

        # CIGAR as array of ops/lengths
        cigars.append(np.array(read.cigartuples, dtype=np.int32) if read.cigartuples else np.zeros((0,2), dtype=np.int32))

        # Only store MD tags if mismatches will be computed
        if {"phagetermini", "assemblycheck"}.intersection(modules):
            if read.has_tag("MD"):
                md_bytes = np.frombuffer(read.get_tag("MD").encode("ascii"), dtype=np.uint8)
                md_list.append(md_bytes)
                md_lengths.append(len(md_bytes))
                has_md.append(True)
            else:
                md_list.append(np.zeros(0, dtype=np.uint8))
                md_lengths.append(0)
                has_md.append(False)

    # Convert all lists to numpy arrays
    ref_starts = np.array(ref_starts, dtype=np.int32)
    ref_ends = np.array(ref_ends, dtype=np.int32)
    query_lengths = np.array(query_lengths, dtype=np.int32)
    template_lengths = np.array(template_lengths, dtype=np.int32)
    is_read1 = np.array(is_read1, dtype=np.bool_)
    is_proper_pair = np.array(is_proper_pair, dtype=np.bool_)
    is_paired = np.array(is_paired, dtype=np.bool_)
    is_reverse = np.array(is_reverse, dtype=np.bool_)
    cigars = np.array(cigars, dtype=object)
    has_md = np.array(has_md, dtype=np.bool_)
    md_lengths = np.array(md_lengths, dtype=np.int32)
    md_list = np.array(md_list, dtype=object)

    return (ref_starts, ref_ends, query_lengths, template_lengths,
            is_read1, is_proper_pair, is_paired, is_reverse,
            cigars, has_md, md_list, md_lengths)

@track_time
def calculating_features_per_contig_per_sample(module_list, bam_file, ref, locus_size, sequencing_type, output_dir, sample_name, ref_name):
    reads_mapped = bam_file.fetch(ref)
    # Save relevant info from bam file
    (ref_starts, ref_ends, query_lengths, template_lengths,
     is_read1, is_proper_pair, is_paired, is_reverse,
     cigars, has_md, md_list, md_lengths) = preprocess_reads(reads_mapped, module_list)

    # Calculate all features
    if "coverage" in module_list:
        get_feature_coverage(ref_starts, ref_ends, locus_size, output_dir, sample_name, ref_name)
    if "phagetermini" in module_list:
        get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, 
                                  sequencing_type, locus_size, output_dir, sample_name, ref_name)
    if "assemblycheck" in module_list:
        reads_mapped = bam_file.fetch(ref)
        get_features_assemblycheck(reads_mapped, locus_size, sequencing_type, output_dir, sample_name, ref_name)
    
### Distributing the calculation for all the contigs in the bam file
def calculating_features_per_sample(module_list, mapping_file, sequencing_type, output_dir, sample_name):
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    references = bam_file.references
    lengths = [l // 2 for l in bam_file.lengths]  # need to divide by 2 because each contig was doubled for mapping to deal with circularity

    for ref, length in zip(references, lengths):
        calculating_features_per_contig_per_sample(module_list, bam_file, ref, length, sequencing_type, output_dir, sample_name, ref)

    bam_file.close()

### Distributing the calculation for all the bam files (samples)
def _process_single_sample(list_modules, bam_file, output_dir):
    sample_name = os.path.basename(bam_file).replace(".bam", "")
    sequencing_type = find_sequencing_type_from_bam(bam_file)

    # Calculate features per contig (can still be parallelized with n_contig_cores)
    calculating_features_per_sample(list_modules, bam_file, sequencing_type, output_dir, sample_name)

    # Temporary debugging
    print_timing_summary()

def calculating_all_features_parallel(list_modules, bam_files, output_dir, n_sample_cores=None):
    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)
    print(f"Using {n_sample_cores} cores to process {len(bam_files)} samples in parallel...", flush=True)

    args_list = [(list_modules, bam, output_dir) for bam in bam_files]

    with Pool(processes=n_sample_cores) as pool:
        pool.starmap(_process_single_sample, args_list)
    print("Finished all samples.", flush=True)

### Main function
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to store output csv files")
    args = parser.parse_args()

    # Get list of BAM files
    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    # Getting list of features requested
    # If starts in feature_list replace it by starts_plus, starts_minus, ends_plus, ends_minus
    requested_modules = args.modules.split(",")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    n_cores = int(args.threads)

    # Calculating values for all requested features from mapping files
    print("### Calculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(requested_modules, bam_files, output_dir, n_cores)

if __name__ == "__main__":
    main()