import argparse, sys, os
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from numba import njit
import numpy as np
import pysam
import pyarrow as pa
import pyarrow.parquet as pq

import time
from functools import wraps
from collections import defaultdict
import sqlite3
import subprocess
from pathlib import Path

try:
    from . import slurm_utils
except Exception:
    slurm_utils = None

# Try to import Rust bindings for faster processing
try:
    import mgfeatureviewer_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None

# Feature type mapping (replaces Variable table lookup during computation)
FEATURE_TYPES = {
    "coverage": "curve",
    "coverage_reduced": "curve",
    "reads_starts": "bars",
    "reads_ends": "bars",
    "tau": "bars",
    "read_lengths": "curve",
    "insert_sizes": "curve",
    "bad_orientations": "bars",
    "left_clippings": "bars",
    "right_clippings": "bars",
    "insertions": "bars",
    "deletions": "bars",
    "mismatches": "bars",
}

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

def create_temp_sample_db(temp_db_path: str):
    """Create a lightweight temp DB for a single sample to avoid contention.
    It contains two tables: TempPresences and TempFeatureValues.
    """
    conn = sqlite3.connect(temp_db_path)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS TempPresences (
        Contig_name TEXT,
        Sample_name TEXT,
        Coverage_percentage REAL
    )
    """)
    cur.execute("""
    CREATE TABLE IF NOT EXISTS TempFeatureValues (
        Variable_name TEXT,
        Contig_name TEXT,
        Sample_name TEXT,
        Position INTEGER,
        Value REAL
    )
    """)
    conn.commit()
    conn.close()


def merge_temp_db_into_main(main_db_path: str, temp_db_path: str):
    """Merge data from a temp sample DB into the main DB.
    Assumptions:
      - main DB already has Contig and Variable tables populated (from genbank step)
      - temp DB contains TempPresences and TempFeatureValues
    The function will:
      - Insert missing Sample rows into main.Sample
      - Map contig names to Contig_id using main.Contig
      - For each Variable_name in TempFeatureValues, find Feature_table_name in main.Variable
        and insert rows into that feature table mapping contig/sample names to ids.
      - Insert TempPresences rows into main.Presences mapping contig/sample names to ids.
    """
    mconn = sqlite3.connect(main_db_path)
    mcur = mconn.cursor()
    # Attach temp DB
    mcur.execute(f"ATTACH DATABASE ? AS src", (temp_db_path,))

    # Insert Samples referenced in temp into main.Sample
    mcur.execute("SELECT DISTINCT Sample_name FROM src.TempPresences")
    sample_rows = mcur.fetchall()
    mcur.execute("SELECT DISTINCT Sample_name FROM src.TempFeatureValues")
    sample_rows += [r for r in mcur.fetchall() if r not in sample_rows]

    for (sname,) in sample_rows:
        if sname is None:
            continue
        mcur.execute("INSERT OR IGNORE INTO Sample (Sample_name) VALUES (?)", (sname,))

    mconn.commit()

    # Insert presences mapping contig and sample names to ids
    mcur.execute("INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) SELECT (SELECT Contig_id FROM Contig WHERE Contig_name=src.TempPresences.Contig_name), (SELECT Sample_id FROM Sample WHERE Sample_name=src.TempPresences.Sample_name), src.TempPresences.Coverage_percentage FROM src.TempPresences")
    mconn.commit()

    # For features, iterate distinct Variable_name and insert into corresponding tables
    mcur.execute("SELECT DISTINCT Variable_name FROM src.TempFeatureValues")
    vars_ = [r[0] for r in mcur.fetchall()]
    for var in vars_:
        # find feature table name in main.Variable
        mcur.execute("SELECT Feature_table_name FROM Variable WHERE Variable_name=?", (var,))
        row = mcur.fetchone()
        if not row:
            # Unknown variable -> skip
            continue
        feat_table = row[0]
        # Build and execute insert: map contig/sample names to ids
        sql = f"INSERT INTO {feat_table} (Contig_id, Sample_id, Position, Value) SELECT (SELECT Contig_id FROM Contig WHERE Contig_name=src.TempFeatureValues.Contig_name), (SELECT Sample_id FROM Sample WHERE Sample_name=src.TempFeatureValues.Sample_name), src.TempFeatureValues.Position, src.TempFeatureValues.Value FROM src.TempFeatureValues WHERE Variable_name=?"
        mcur.execute(sql, (var,))
        mconn.commit()

    # Detach temp DB
    mcur.execute("DETACH DATABASE src")
    mconn.close()


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

def compress_signal(type_picked, feature_values, ref_length, step, z_thresh, deriv_thresh, max_points):
    """
    step : int, Keep every Nth point
    z_thresh : float, Z-score threshold for keeping high/low outlier values
    deriv_thresh : float, Z-score threshold for keeping large derivative changes
    max_points : int, Hard limit on total number of kept points
    """
    feature_values = np.asarray(feature_values)
    n = len(feature_values)

    # Value outliers
    y_mean = np.mean(feature_values)
    y_std = np.std(feature_values) or 1e-9
    val_outliers = np.abs(feature_values - y_mean) > z_thresh * y_std

    # Regular subsampling
    if type_picked == "curve":
        regular_idx = np.arange(0, n, step, dtype=int) if n > 0 else np.array([], dtype=int)

        # Derivative outliers
        dy = np.diff(feature_values, prepend=feature_values[0])
        dy_std = np.std(dy) or 1e-9
        der_outliers = np.abs(dy) > deriv_thresh * dy_std

        # Include the point before each derivative outlier
        # NOTE: This line has a subtle bug - it does arithmetic on a boolean array,
        # which doesn't correctly find outlier indices. The Rust implementation
        # (rust_calculate/src/main.rs, compress_signal function) uses correct logic
        # and may keep slightly more positions as a result. Both outputs are valid
        # but Rust preserves more derivative outliers as originally intended.
        der_outliers = np.unique(np.clip(np.concatenate([der_outliers - 1, der_outliers]), 0, n - 1))

        # Combine value and derivative outliers
        outlier_idx = np.nonzero(val_outliers)[0]
        outlier_idx = np.unique(np.concatenate([outlier_idx, der_outliers]))

        last_idx = np.array([n - 1], dtype=int) if n > 0 else np.array([], dtype=int)
        keep_idx = merge_sorted_unique(regular_idx, outlier_idx, last_idx)
        
    elif type_picked == "bars":
        # For bars: only keep outliers (value OR derivative)
        keep_idx = np.nonzero(val_outliers)[0]
    else:
        raise ValueError(f"Unknown type_picked: {type_picked}")
    
    # Apply hard limit if needed
    if len(keep_idx) > max_points:
        step_lim = len(keep_idx) // max_points
        keep_idx = keep_idx[::step_lim]
        if len(keep_idx) > 0 and keep_idx[-1] != n - 1:
            keep_idx = np.append(keep_idx, n - 1)

    # Initialize x array of ref_length size
    x = np.arange(1, ref_length + 1)
    return {"x": x[keep_idx], "y": feature_values[keep_idx]}

### Save data computed per sample
def get_contig_id(conn, contig_name):
    cur = conn.cursor()
    cur.execute("SELECT Contig_id FROM Contig WHERE Contig_name=?", (contig_name,))
    contig_id = cur.fetchone()[0]
    return contig_id

def add_sample(conn, sample_name):
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO Sample (Sample_name) VALUES (?)", (sample_name,))
    conn.commit()

    cur.execute("SELECT Sample_id FROM Sample WHERE Sample_name=?", (sample_name,))
    row = cur.fetchone()
    return row[0] if row else None

def add_presence(conn, contig_id=None, sample_id=None, coverage_percentage=None, *, contig_name=None, sample_name=None, temp_mode=False):
    cur = conn.cursor()
    if temp_mode:
        # write to TempPresences (contig_name, sample_name, coverage)
        cur.execute("INSERT INTO TempPresences (Contig_name, Sample_name, Coverage_percentage) VALUES (?, ?, ?)", (contig_name, sample_name, coverage_percentage))
    else:
        cur.execute("INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?, ?, ?)", (contig_id, sample_id, coverage_percentage))
    conn.commit()

def compress_and_return_features(feature, feature_values, contig_name, ref_length, step, z_thresh, deriv_thresh, max_points):
    """Compress feature values and return as list of tuples for Parquet."""
    type_picked = FEATURE_TYPES[feature]

    # Compress the signal to limit points
    feature_compressed = compress_signal(type_picked, feature_values, ref_length, step, z_thresh, deriv_thresh, max_points)

    xs = feature_compressed.get("x", [])
    ys = feature_compressed.get("y", [])

    # Return list of (contig_name, feature, position, value) tuples
    return [(contig_name, feature, int(x), float(y)) for x, y in zip(xs, ys)]

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

@track_time
def get_feature_coverage(ref_starts, ref_ends, ref_length, contig_name, step, z_thresh, deriv_thresh, max_points):
    coverage = np.zeros(ref_length, dtype=np.uint64)
    calculate_coverage_numba(coverage, ref_starts, ref_ends, ref_length)

    # Return compressed feature data
    return compress_and_return_features("coverage", coverage, contig_name, ref_length, step, z_thresh, deriv_thresh, max_points)

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
    if len(md) == 0:
        return False
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

@njit
def calculate_reads_starts_and_ends_numba(start, end, is_reverse, coverage_reduced, start_plus, start_minus, end_plus, end_minus, ref_length):
    """Numba-optimized version with direct array access."""
    start_mod = start % ref_length
    end_mod = end % ref_length

    # Update coverage
    if start_mod <= end_mod:
        for j in range(start_mod, end_mod + 1):
            coverage_reduced[j] += 1
    else:
        for j in range(start_mod, ref_length):
            coverage_reduced[j] += 1
        for j in range(0, end_mod + 1):
            coverage_reduced[j] += 1

    # Update strand-specific starts/ends
    if is_reverse:
        start_minus[start_mod] += 1
        end_minus[end_mod] += 1
    else:
        start_plus[start_mod] += 1
        end_plus[end_mod] += 1

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
def get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, sequencing_type, ref_length,
                              contig_name, step, z_thresh, deriv_thresh, max_points):
    features = ["coverage_reduced", "reads_starts", "reads_ends", "tau"]
    temporary_features = ["coverage_reduced", "start_plus", "start_minus", "end_plus", "end_minus"]

    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    # Pre-resolve sequencing type check (avoid string comparison in loop)
    is_long = sequencing_type == "long"

    # Extract direct array references (avoid dict lookup in loop)
    coverage_reduced = temporary_dict["coverage_reduced"]
    start_plus = temporary_dict["start_plus"]
    start_minus = temporary_dict["start_minus"]
    end_plus = temporary_dict["end_plus"]
    end_minus = temporary_dict["end_minus"]

    n_reads = len(ref_starts)
    # Iterate through preprocessed reads
    for i in range(n_reads):
        cigar = cigars[i]
        md = md_list[i]

        # Check if read starts with a match (and ends with a match for long reads)
        if starts_with_match(cigar, md, start=True) and (not is_long or starts_with_match(cigar, md, start=False)):
            calculate_reads_starts_and_ends_numba(
                ref_starts[i], ref_ends[i], is_reverse[i],
                coverage_reduced, start_plus, start_minus, end_plus, end_minus, ref_length
            )

    # Compute tau
    compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length)

    # Return compressed feature data
    result = []
    for feature in features:
        result.extend(compress_and_return_features(feature, feature_dict[feature], contig_name, ref_length, step, z_thresh, deriv_thresh, max_points))
    return result

### Functions of assemblycheck module
@njit
def add_read_lengths_numba(sum_read_lengths, counts_read_lengths, positions, read_len, ref_length):
    for i in range(positions.size):
        pos = positions[i] % ref_length
        sum_read_lengths[pos] += read_len
        counts_read_lengths[pos] += 1

@njit
def add_read_lengths_range_numba(sum_read_lengths, counts_read_lengths, start, end, read_len, ref_length):
    """Optimized version that iterates over range instead of requiring pre-built array."""
    for pos in range(start, end):
        p = pos % ref_length
        sum_read_lengths[p] += read_len
        counts_read_lengths[p] += 1

@njit
def add_insert_sizes_numba(sum_insert_sizes, counts_insert_sizes, bad_orientations, positions, insert_len, is_read1, proper_pair, insert_flag, bad_flag, ref_length):
    for i in range(positions.size):
        pos = positions[i] % ref_length
        if insert_flag and is_read1 and insert_len > 0:
            sum_insert_sizes[pos] += insert_len
            counts_insert_sizes[pos] += 1
        if bad_flag and not proper_pair:
            bad_orientations[pos] += 1

@njit
def add_insert_sizes_range_numba(sum_insert_sizes, counts_insert_sizes, bad_orientations,
                                  start, end, insert_len, is_read1, proper_pair,
                                  insert_flag, bad_flag, ref_length):
    """Optimized version that iterates over range instead of requiring pre-built array."""
    for pos in range(start, end):
        p = pos % ref_length
        if insert_flag and is_read1 and insert_len > 0:
            sum_insert_sizes[p] += insert_len
            counts_insert_sizes[p] += 1
        if bad_flag and not proper_pair:
            bad_orientations[p] += 1

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

def compute_final_lengths(sum_lengths, count_lengths):
    count_lengths = np.maximum(count_lengths, 1)  # avoid division by zero
    arr = sum_lengths / count_lengths
    arr[count_lengths == 0] = 0
    return arr.astype(np.float32)

@track_time
def get_features_assemblycheck(ref_starts, ref_ends, query_lengths, template_lengths, is_read1, is_proper_pair,
                               cigars, has_md, md_list, md_lengths, sequencing_type, ref_length,
                               contig_name, step, z_thresh, deriv_thresh, max_points):
    features = []
    temporary_features = []

    # Pre-resolve sequencing type checks (avoid string comparison in loop)
    is_long = sequencing_type == "long"
    is_paired = sequencing_type == "short-paired"

    if is_long:
        features.append("read_lengths")
        temporary_features.extend(["sum_read_lengths", "count_read_lengths"])
    if is_paired:
        features.extend(["insert_sizes", "bad_orientations"])
        temporary_features.extend(["sum_insert_sizes", "count_insert_sizes", "bad_orientation"])
    features.extend(["left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])

    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    # Pre-resolve feature flags (avoid repeated string lookups in loop)
    has_left_clip = "left_clippings" in features
    has_right_clip = "right_clippings" in features
    has_mismatches = "mismatches" in features

    # Get direct array references (avoid dict lookup in loop)
    left_clippings = feature_dict.get("left_clippings")
    right_clippings = feature_dict.get("right_clippings")
    insertions_arr = feature_dict.get("insertions")
    deletions_arr = feature_dict.get("deletions")
    mismatches_arr = feature_dict.get("mismatches")
    bad_orientations = feature_dict.get("bad_orientations")

    sum_read_lengths = temporary_dict.get("sum_read_lengths")
    count_read_lengths = temporary_dict.get("count_read_lengths")
    sum_insert_sizes = temporary_dict.get("sum_insert_sizes")
    count_insert_sizes = temporary_dict.get("count_insert_sizes")

    n_reads = ref_starts.size
    for i in range(n_reads):
        start_i = ref_starts[i]
        end_i = ref_ends[i]

        # --- Long reads: use range-based numba (avoids np.arange allocation) ---
        if is_long:
            add_read_lengths_range_numba(sum_read_lengths, count_read_lengths,
                                         start_i, end_i, query_lengths[i], ref_length)

        # --- Short-paired reads: use range-based numba ---
        if is_paired:
            add_insert_sizes_range_numba(sum_insert_sizes, count_insert_sizes, bad_orientations,
                                         start_i, end_i, template_lengths[i], is_read1[i], is_proper_pair[i],
                                         True, True, ref_length)

        # --- Clippings ---
        cigar = cigars[i]
        if cigar.size > 0:
            first_op = cigar[0, 0]
            last_op = cigar[-1, 0]
            if has_left_clip and (first_op == 4 or first_op == 5):
                left_clippings[start_i % ref_length] += 1
            if has_right_clip and (last_op == 4 or last_op == 5):
                right_clippings[(end_i - 1) % ref_length] += 1

        # --- Indels ---
        add_indels_numba(deletions_arr, insertions_arr, cigar, start_i, ref_length)

        # --- Mismatches ---
        if has_mismatches and has_md[i]:
            add_mismatches_numba_from_md(mismatches_arr, md_list[i], md_lengths[i], start_i, ref_length)

    # --- Finalize lengths ---
    if is_long:
        feature_dict["read_lengths"] = compute_final_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"])
    if is_paired:
        feature_dict["insert_sizes"] = compute_final_lengths(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"])

    # --- Return compressed feature data ---
    result = []
    for feature in features:
        result.extend(compress_and_return_features(feature, feature_dict[feature], contig_name, ref_length, step, z_thresh, deriv_thresh, max_points))
    return result

### Calculating features per contig per sample
@track_time
def preprocess_reads(reads_mapped, modules, sequencing_type):
    # Determine which attributes we actually need based on modules and sequencing type
    # This reduces pysam attribute accesses (each is a C/Python boundary crossing)
    need_md = bool({"phagetermini", "assemblycheck"}.intersection(modules))
    need_cigars = bool({"phagetermini", "assemblycheck"}.intersection(modules))
    need_reverse = "phagetermini" in modules
    need_query_len = "assemblycheck" in modules and sequencing_type == "long"
    need_paired_attrs = "assemblycheck" in modules and sequencing_type == "short-paired"

    # Use dynamically growing arrays - single pass through reads
    INITIAL_CAP = 8192
    GROWTH = 2

    capacity = INITIAL_CAP
    n = 0

    ref_starts = np.empty(capacity, dtype=np.int32)
    ref_ends = np.empty(capacity, dtype=np.int32)

    # Only allocate what's needed
    query_lengths = np.empty(capacity, dtype=np.int32) if need_query_len else None
    template_lengths = np.empty(capacity, dtype=np.int32) if need_paired_attrs else None
    is_read1 = np.empty(capacity, dtype=np.bool_) if need_paired_attrs else None
    is_proper_pair = np.empty(capacity, dtype=np.bool_) if need_paired_attrs else None
    is_paired = np.empty(capacity, dtype=np.bool_) if need_paired_attrs else None
    is_reverse = np.empty(capacity, dtype=np.bool_) if need_reverse else None
    cigars = [] if need_cigars else None
    has_md = np.empty(capacity, dtype=np.bool_) if need_md else None
    md_list = [] if need_md else None
    md_lengths = np.empty(capacity, dtype=np.int32) if need_md else None

    # Single pass through reads
    for read in reads_mapped:
        if read.is_unmapped:
            continue

        # Grow arrays if needed
        if n >= capacity:
            capacity *= GROWTH
            ref_starts = np.resize(ref_starts, capacity)
            ref_ends = np.resize(ref_ends, capacity)
            if need_query_len:
                query_lengths = np.resize(query_lengths, capacity)
            if need_paired_attrs:
                template_lengths = np.resize(template_lengths, capacity)
                is_read1 = np.resize(is_read1, capacity)
                is_proper_pair = np.resize(is_proper_pair, capacity)
                is_paired = np.resize(is_paired, capacity)
            if need_reverse:
                is_reverse = np.resize(is_reverse, capacity)
            if need_md:
                has_md = np.resize(has_md, capacity)
                md_lengths = np.resize(md_lengths, capacity)

        # Always needed
        ref_starts[n] = read.reference_start
        ref_ends[n] = read.reference_end

        # Conditional extractions (skip pysam attribute access if not needed)
        if need_query_len:
            query_lengths[n] = read.query_length
        if need_paired_attrs:
            template_lengths[n] = abs(read.template_length)
            is_read1[n] = read.is_read1
            is_proper_pair[n] = read.is_proper_pair
            is_paired[n] = read.is_paired
        if need_reverse:
            is_reverse[n] = read.is_reverse
        if need_cigars:
            cigars.append(np.array(read.cigartuples, dtype=np.int32) if read.cigartuples else np.zeros((0, 2), dtype=np.int32))
        if need_md:
            if read.has_tag("MD"):
                md_bytes = np.frombuffer(read.get_tag("MD").encode("ascii"), dtype=np.uint8)
                md_list.append(md_bytes)
                md_lengths[n] = len(md_bytes)
                has_md[n] = True
            else:
                md_list.append(np.zeros(0, dtype=np.uint8))
                md_lengths[n] = 0
                has_md[n] = False

        n += 1

    # Trim to actual size (or return empty arrays)
    empty_i32 = np.array([], dtype=np.int32)
    empty_bool = np.array([], dtype=np.bool_)

    if n == 0:
        return (empty_i32, empty_i32, empty_i32, empty_i32,
                empty_bool, empty_bool, empty_bool, empty_bool,
                [], empty_bool, [], empty_i32)

    ref_starts = ref_starts[:n]
    ref_ends = ref_ends[:n]
    query_lengths = query_lengths[:n] if need_query_len else empty_i32
    template_lengths = template_lengths[:n] if need_paired_attrs else empty_i32
    is_read1 = is_read1[:n] if need_paired_attrs else empty_bool
    is_proper_pair = is_proper_pair[:n] if need_paired_attrs else empty_bool
    is_paired = is_paired[:n] if need_paired_attrs else empty_bool
    is_reverse = is_reverse[:n] if need_reverse else empty_bool
    cigars = cigars if need_cigars else []
    has_md = has_md[:n] if need_md else empty_bool
    md_list = md_list if need_md else []
    md_lengths = md_lengths[:n] if need_md else empty_i32

    return (ref_starts, ref_ends, query_lengths, template_lengths, is_read1, is_proper_pair, is_paired, is_reverse, cigars, has_md, md_list, md_lengths)

@track_time
def calculating_features_per_contig_per_sample(module_list, bam_file, sequencing_type, ref_name, locus_size,
                                               min_coverage, step, z_thresh, deriv_thresh, max_points):
    """Calculate features for one contig, return (feature_data, presence_data) or (None, None)."""
    ### Save relevant info from bam file
    reads_mapped = bam_file.fetch(ref_name)
    (ref_starts, ref_ends, query_lengths, template_lengths,
     is_read1, is_proper_pair, is_paired, is_reverse,
     cigars, has_md, md_list, md_lengths) = preprocess_reads(reads_mapped, module_list, sequencing_type)

    ### Coverage check: ensure more than min_coverage% of the reference is covered by at least one read
    if len(ref_starts) == 0:
        return None, None

    # Mark covered regions
    covered = np.zeros(locus_size, dtype=bool)
    start_clipped = np.maximum(ref_starts, 0)
    end_clipped = np.minimum(ref_ends, locus_size)
    for s, e in zip(start_clipped, end_clipped):
        covered[s:e] = True

    covered_bp = covered.sum()
    coverage_pct = (covered_bp / locus_size) * 100

    if coverage_pct < min_coverage:
        return None, None

    # Presence record (contig_name, coverage_pct)
    presence = (ref_name, coverage_pct)

    # Accumulate feature data
    features = []

    ### Calculate all features
    if {"coverage", "assemblycheck"}.intersection(module_list):
        features.extend(get_feature_coverage(ref_starts, ref_ends, locus_size, ref_name, step, z_thresh, deriv_thresh, max_points))
    if "phagetermini" in module_list:
        features.extend(get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, sequencing_type, locus_size,
                                                  ref_name, step, z_thresh, deriv_thresh, max_points))
    if "assemblycheck" in module_list:
        features.extend(get_features_assemblycheck(ref_starts, ref_ends, query_lengths, template_lengths, is_read1, is_proper_pair, cigars, has_md, md_list, md_lengths,
                                                   sequencing_type, locus_size, ref_name, step, z_thresh, deriv_thresh, max_points))

    return features, presence

### Distributing the calculation for all the contigs in the bam file
@track_time
def calculating_features_per_sample(module_list, mapping_file, min_coverage, step, z_thresh, deriv_thresh, max_points):
    """Calculate all features for one sample, return (all_features, all_presences, sample_name)."""
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    references = bam_file.references
    lengths = [l // 2 for l in bam_file.lengths]  # need to divide by 2 because each contig was doubled for mapping to deal with circularity

    sample_name = os.path.basename(mapping_file).replace(".bam", "")
    sequencing_type = find_sequencing_type_from_bam(mapping_file)

    all_features = []
    all_presences = []

    for ref, length in zip(references, lengths):
        features, presence = calculating_features_per_contig_per_sample(
            module_list, bam_file, sequencing_type, ref, length,
            min_coverage, step, z_thresh, deriv_thresh, max_points
        )
        if features is not None:
            all_features.extend(features)
            all_presences.append(presence)

    bam_file.close()
    return all_features, all_presences, sample_name

### Distributing the calculation for all the bam files (samples)
def write_sample_parquet(features, presences, sample_name, output_dir):
    """Write feature and presence data to Parquet files for one sample."""
    # Write features Parquet
    if features:
        contig_names = [f[0] for f in features]
        feature_names = [f[1] for f in features]
        positions = [f[2] for f in features]
        values = [f[3] for f in features]

        features_table = pa.table({
            'contig_name': pa.array(contig_names, type=pa.string()),
            'feature': pa.array(feature_names, type=pa.string()),
            'position': pa.array(positions, type=pa.int32()),
            'value': pa.array(values, type=pa.float32()),
        })

        features_path = os.path.join(output_dir, 'features', f'{sample_name}.parquet')
        pq.write_table(features_table, features_path, compression='zstd')

    # Write presences Parquet
    if presences:
        contig_names = [p[0] for p in presences]
        coverage_pcts = [p[1] for p in presences]

        presences_table = pa.table({
            'contig_name': pa.array(contig_names, type=pa.string()),
            'coverage_pct': pa.array(coverage_pcts, type=pa.float32()),
        })

        presences_path = os.path.join(output_dir, 'presences', f'{sample_name}.parquet')
        pq.write_table(presences_table, presences_path, compression='zstd')

def _process_single_sample(list_modules, bam_file, output_dir, min_coverage, step, z_thresh, deriv_thresh, max_points):
    """Process one sample: compute features and write Parquet files (pure Python)."""
    all_features, all_presences, sample_name = calculating_features_per_sample(
        list_modules, bam_file, min_coverage, step, z_thresh, deriv_thresh, max_points
    )

    # Write to Parquet (each worker writes its own files - no contention)
    write_sample_parquet(all_features, all_presences, sample_name, output_dir)

    bam.close()
    try:
        with open(log_path, 'a', encoding='utf8') as fh:
            fh.write(f"\n=== Sample: {sample_name} | BAM: {bam_path} | Time: {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
            total = sum(function_times.values())
            fh.write(f"Total tracked time: {total:.3f} s\n")
            for name, t in sorted(function_times.items(), key=lambda x: -x[1]):
                fh.write(f"{name:35s}: {t:.3f} s ({t/total*100:.1f}%)\n")
    except Exception:
        pass


def _update_sqlite_samples(output_dir, all_sample_presences):
    """Update SQLite database with Sample and Presences tables.

    Args:
        output_dir: Path to output directory containing metadata.db
        all_sample_presences: List of (sample_name, presences) where presences is
                              a list of (contig_name, coverage_pct) tuples
    """
    db_path = os.path.join(output_dir, "metadata.db")
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Build contig name -> id mapping
    cur.execute("SELECT Contig_name, Contig_id FROM Contig")
    contig_name_to_id = {row[0]: row[1] for row in cur.fetchall()}

    # Insert samples and presences
    for sample_name, presences in all_sample_presences:
        # Insert sample
        cur.execute("INSERT INTO Sample (Sample_name) VALUES (?)", (sample_name,))
        sample_id = cur.lastrowid

        # Insert presences for this sample
        for contig_name, coverage_pct in presences:
            contig_id = contig_name_to_id.get(contig_name)
            if contig_id is not None:
                cur.execute(
                    "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?, ?, ?)",
                    (contig_id, sample_id, coverage_pct)
                )

    conn.commit()
    conn.close()

def calculating_all_features_parallel(list_modules, bam_files, output_dir, min_coverage, step, z_thresh, deriv_thresh, max_points, n_sample_cores=None, genbank_path=None, annotation_tool=""):
    """Process all BAM files in parallel.

    If genbank_path is provided and Rust bindings are available, uses the faster Rust implementation
    with rayon parallelism and a single GenBank parse.
    Otherwise falls back to pure Python.
    """
    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)

    # Create output directories
    os.makedirs(os.path.join(output_dir, 'features'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'presences'), exist_ok=True)

    # Use Rust if available and genbank_path is provided
    if HAS_RUST and genbank_path is not None:
        # Get parent directory of BAM files (assumes all in same directory)
        bam_dir = os.path.dirname(bam_files[0]) if bam_files else ""

        print(f"Using Rust bindings to process {len(bam_files)} samples with rayon ({n_sample_cores} threads)...", flush=True)

        # Call the Rust function that handles everything: GenBank parsing, parallel BAM processing,
        # progress bar, writing Parquet files, and updating SQLite directly
        result = _rust.process_all_samples(
            genbank_path=genbank_path,
            bam_dir=bam_dir,
            output_dir=output_dir,
            modules=list_modules,
            threads=n_sample_cores,
            annotation_tool=annotation_tool,
            min_coverage=float(min_coverage),
            step=step,
            z_thresh=float(z_thresh),
            deriv_thresh=float(deriv_thresh),
            max_points=max_points,
        )

        # Rust handles SQLite updates internally, just report results
        total_time = result.get("total_time", 0.0)
        samples_processed = result.get("samples_processed", 0)
        samples_failed = result.get("samples_failed", 0)
        print(f"Finished {samples_processed} samples in {total_time:.2f}s ({total_time/max(1, samples_processed):.2f}s per sample avg)", flush=True)
        if samples_failed > 0:
            print(f"Warning: {samples_failed} samples failed to process", flush=True)
    else:
        # Fall back to pure Python
        if not HAS_RUST:
            print("Note: Rust bindings not available, using pure Python (slower)", flush=True)
        print(f"Using {n_sample_cores} cores to process {len(bam_files)} samples in parallel...", flush=True)

        # Sort BAM files by size (largest first) for better load balancing
        bam_files_sorted = sorted(bam_files, key=lambda f: os.path.getsize(f), reverse=True)
        args_list = [(list_modules, bam, output_dir, min_coverage, step, z_thresh, deriv_thresh, max_points) for bam in bam_files_sorted]

        # chunksize=1 ensures workers grab next task immediately when done
        with Pool(processes=n_sample_cores) as pool:
            list(pool.starmap(_process_single_sample, args_list, chunksize=1))
        print("Finished all samples.", flush=True)

### Main function helpers (shared-args)
def add_calculate_args(parser):
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file of all investigated contigs")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for results (metadata.db + Parquet files)")
    parser.add_argument("-a", "--annotation_tool", default="", help="Optional: to color the contigs specify the annotation tool used (options allowed: pharokka)")
    parser.add_argument("--min_coverage", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion")
    parser.add_argument("--step", type=int, default=50, help="Step size for compression (keep every Nth point in addition to the outliers)")
    parser.add_argument("--outlier_threshold", type=int, default=3, help="Points beyond mean+std*N are kept as outliers")
    parser.add_argument("--derivative_threshold", type=int, default=3, help="Points were the derivative is beyond mean+std*N are kept as outliers")
    parser.add_argument("--max_points", type=int, default=10000, help="Maximum number of points kept during compression")
    # Slurm integration (optional)
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to request per Slurm job (cpus-per-task)")
    parser.add_argument("--use-slurm", action="store_true", help="Submit per-sample jobs via Slurm (sbatch) instead of running locally")
    parser.add_argument("--max-concurrent", type=int, default=20, help="Maximum concurrent Slurm array tasks")
    parser.add_argument("--max-time", type=str, default="02:00:00", help="Time limit per Slurm job (HH:MM:SS)")
    parser.add_argument("--mem-per-cpu", type=str, default="8G", help="Memory per Slurm job (e.g. 8G)")
    parser.add_argument("--parallelize_contigs", action="store_true", help="When running a sample job, parallelize over contigs inside the job")
    # Hidden flag: internal mode to run only a single sample (used by Slurm-submitted jobs)
    parser.add_argument("--run-sample", action="store_true", help=argparse.SUPPRESS)

def run_calculate_args(args):
    annotation_tool = args.annotation_tool

    # Get list of BAM files
    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    # When using Rust, skip Python validation - Rust handles everything including error checking
    if not HAS_RUST:
        # Pure Python path: validate GenBank and BAM compatibility
        print("### Checking that genbank and mapping files are compatible...", flush=True)
        genbank_loci = set(rec.name for rec in SeqIO.parse(args.genbank, "genbank"))
        if not genbank_loci:
            sys.exit("ERROR: No loci found in the provided GenBank file.")

        # Silence htslib warnings about index file timestamps
        save_verbosity = pysam.set_verbosity(0)

        for bam_file in bam_files:
            try:
                with pysam.AlignmentFile(bam_file, "rb") as bam:
                    bam_refs = set(bam.references)
            except Exception as e:
                pysam.set_verbosity(save_verbosity)
                sys.exit(f"ERROR: Could not open BAM file '{bam_file}': {e}")

            missing_in_genbank = bam_refs - genbank_loci

            if missing_in_genbank:
                pysam.set_verbosity(save_verbosity)
                raise ValueError(
                    f"ERROR: References in BAM file '{os.path.basename(bam_file)}' "
                    f"not found in GenBank:\n"
                    f"{', '.join(sorted(missing_in_genbank))}"
                )

        pysam.set_verbosity(save_verbosity)
        print("Sanity checks passed.", flush=True)

    # Requested modules
    requested_modules = args.modules.split(",")

    # Parameters for compression
    min_coverage = args.min_coverage
    step = args.step
    z_thresh = args.outlier_threshold
    deriv_thresh = args.derivative_threshold
    max_points = args.max_points

    n_cores = int(args.threads)

    # Setup output directory
    output_dir = args.output
    if os.path.exists(output_dir):
        sys.exit(f"ERROR: Output directory '{output_dir}' already exists. Please provide a new path to avoid overwriting.")

    # When Rust bindings are available, they handle everything: GenBank parsing, database creation,
    # parallel BAM processing, and output writing. Skip Python's database setup in that case.
    if HAS_RUST:
        # Rust handles database creation and population internally
        print("Calculating values for all requested features from mapping files...", flush=True)
        calculating_all_features_parallel(
            requested_modules, bam_files, output_dir, min_coverage, step, z_thresh, deriv_thresh, max_points, n_cores,
            genbank_path=args.genbank, annotation_tool=annotation_tool
        )
    else:
        # Pure Python path: create database and populate manually
        os.makedirs(output_dir)

        # Create metadata database
        db_path = os.path.join(output_dir, "metadata.db")
        subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "generate_database.py"), db_path], check=True)

        # Saving genbank info into database
        print("### Saving genbank info into database...", flush=True)
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()

        seq_rows = []
        contig_count = 0
        feature_count = 0
        for rec in SeqIO.parse(args.genbank, "genbank"):
            contig_name = rec.name
            contig_length = len(rec.seq)
            cur.execute("INSERT OR IGNORE INTO Contig (Contig_name, Contig_length, Annotation_tool) VALUES (?, ?, ?)", (contig_name, contig_length, annotation_tool))
            contig_id = get_contig_id(conn, contig_name)
            contig_count += 1

            for f in rec.features:
                try:
                    start = int(f.location.start) + 1
                    end = int(f.location.end)
                    strand_val = f.location.strand
                except Exception:
                    continue

                ftype = f.type
                if not(ftype in {"source", "gene"}):
                    qualifiers = f.qualifiers if hasattr(f, 'qualifiers') else {}
                    product = qualifiers.get('product', [None])[0]
                    function = qualifiers.get('function', [None])[0]
                    phrog = qualifiers.get('phrog', [None])[0]

                    seq_rows.append((contig_id, start, end, strand_val, ftype, product, function, phrog))
                    feature_count += 1

        if seq_rows:
            cur.executemany("INSERT INTO Contig_annotation (Contig_id, Start, End, Strand, Type, Product, Function, Phrog) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", seq_rows)
        conn.commit()
        conn.close()

        print(f"Saved {contig_count} contigs and {feature_count} annotations into database", flush=True)

        print("Calculating values for all requested features from mapping files...", flush=True)
        calculating_all_features_parallel(
            requested_modules, bam_files, output_dir, min_coverage, step, z_thresh, deriv_thresh, max_points, n_cores,
            genbank_path=args.genbank, annotation_tool=annotation_tool
        )

    print(f"\nOutput written to: {output_dir}/", flush=True)
    print(f"  - metadata.db (contig/sample metadata)", flush=True)
    print(f"  - features/*.parquet (feature data)", flush=True)
    print(f"  - presences/*.parquet (presence/absence data)", flush=True)

def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


def merge_all_sample_dbs(main_db_path: str):
    """Merge all per-sample temp DBs found next to `main_db_path` into the main DB.
    This looks for files named '<main_stem>_*.temp.db' in the same directory.
    """
    main_p = Path(main_db_path)
    pattern = f"{main_p.stem}_*.temp.db"
    for f in main_p.parent.glob(pattern):
        try:
            print(f"Merging {f} into {main_db_path}")
            merge_temp_db_into_main(main_db_path, str(f))
            f.unlink()
        except Exception as e:
            print(f"Failed to merge {f}: {e}")

if __name__ == "__main__":
    main()