import argparse, sys, os
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from numba import njit
import numpy as np
import pysam

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

def compress_and_write_features(feature, feature_values, db_conn, sample_id=None, contig_id=None, ref_length=None, step=None, z_thresh=None, deriv_thresh=None, max_points=None, *, sample_name=None, contig_name=None, temp_mode=False):
    # Read variable metadata from DB instead of constants
    cur = db_conn.cursor()

    # When writing into a temp DB, we do not need feature_table resolution
    if not temp_mode:
        cur.execute("SELECT Type, Feature_table_name FROM Variable WHERE Variable_name=?", (feature,))
        row = cur.fetchone()
        type_picked, feature_table = row
    else:
        # For temp mode, we still need the variable Type to compress appropriately
        cur.execute("SELECT Type FROM Variable WHERE Variable_name=?", (feature,))
        row = cur.fetchone()
        type_picked = row[0] if row else 'curve'

    # Compress the signal to limit points
    feature_compressed = compress_signal(type_picked, feature_values, ref_length, step, z_thresh, deriv_thresh, max_points)

    xs = feature_compressed.get("x", [])
    ys = feature_compressed.get("y", [])

    if temp_mode:
        # Insert into generic TempFeatureValues table for later merging
        to_insert = [(feature, contig_name, sample_name, int(x), float(y)) for x, y in zip(xs, ys)]
        if to_insert:
            cur.executemany("INSERT INTO TempFeatureValues (Variable_name, Contig_name, Sample_name, Position, Value) VALUES (?, ?, ?, ?, ?)", to_insert)
            db_conn.commit()
    else:
        to_insert = [(contig_id, sample_id, int(x), float(y)) for x, y in zip(xs, ys)]
        if to_insert:
            cur.executemany(f"INSERT INTO {feature_table} (Contig_id, Sample_id, Position, Value) VALUES (?, ?, ?, ?)", to_insert)
            db_conn.commit()

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
def get_feature_coverage(ref_starts, ref_ends, ref_length, db_conn, sample_id=None, contig_id=None, step=None, z_thresh=None, deriv_thresh=None, max_points=None, *, sample_name=None, contig_name=None, temp_mode=False):
    coverage = np.zeros(ref_length, dtype=np.uint64)
    calculate_coverage_numba(coverage, ref_starts, ref_ends, ref_length)

    # Save into DB
    compress_and_write_features("coverage", coverage, db_conn, sample_id, contig_id, ref_length, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)

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
def get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, sequencing_type, ref_length, 
                              db_conn, sample_id=None, contig_id=None, step=None, z_thresh=None, deriv_thresh=None, max_points=None, *, sample_name=None, contig_name=None, temp_mode=False):
    features = ["coverage_reduced", "reads_starts", "reads_ends", "tau"]
    temporary_features = ["coverage_reduced", "start_plus", "start_minus", "end_plus", "end_minus"]

    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    # Iterate through preprocessed reads
    for i in range(len(ref_starts)):
        cigar = cigars[i]
        md = md_list[i]

        # Check if read starts with a match (and ends with a match for long reads)
        if starts_with_match(cigar, md, start=True) and (sequencing_type != "long" or starts_with_match(cigar, md, start=False)):
            start = ref_starts[i]
            end = ref_ends[i]
            is_reverse_flag = is_reverse[i]
            calculate_reads_starts_and_ends(start, end, is_reverse_flag, temporary_dict, ref_length)

    # Compute tau
    compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length)

    # Save features into DB
    for feature in features:
        compress_and_write_features(feature, feature_dict[feature], db_conn, sample_id, contig_id, ref_length, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)

### Functions of assemblycheck module
@njit
def add_read_lengths_numba(sum_read_lengths, counts_read_lengths, positions, read_len, ref_length):
    for i in range(positions.size):
        pos = positions[i] % ref_length
        sum_read_lengths[pos] += read_len
        counts_read_lengths[pos] += 1

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
                               db_conn, sample_id=None, contig_id=None, step=None, z_thresh=None, deriv_thresh=None, max_points=None, *, sample_name=None, contig_name=None, temp_mode=False):
    features = []
    temporary_features = []

    if sequencing_type == "long":
        features.append("read_lengths")
        temporary_features.extend(["sum_read_lengths", "count_read_lengths"])
    if sequencing_type == "short-paired":
        features.extend(["insert_sizes", "bad_orientations"])
        temporary_features.extend(["sum_insert_sizes", "count_insert_sizes", "bad_orientation"])
    features.extend(["left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])

    feature_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in features}
    temporary_dict = {f: np.zeros(ref_length, dtype=np.uint64) for f in temporary_features}

    n_reads = ref_starts.size
    for i in range(n_reads):
        # --- Long reads ---
        if sequencing_type == "long":
            positions = np.arange(ref_starts[i], ref_ends[i], dtype=np.int32) % ref_length
            add_read_lengths_numba(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"],
                                   positions, query_lengths[i], ref_length)

        # --- Short-paired reads ---
        if sequencing_type == "short-paired":
            positions = np.arange(ref_starts[i], ref_ends[i], dtype=np.int32) % ref_length
            insert_len = template_lengths[i]
            insert_flag = 1 if "insert_sizes" in features else 0
            bad_flag = 1 if "bad_orientations" in features else 0
            add_insert_sizes_numba(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"],
                                   feature_dict.get("bad_orientations", None),
                                   positions, insert_len, is_read1[i], is_proper_pair[i], insert_flag, bad_flag, ref_length)

        # --- Clippings ---
        if cigars[i].size > 0:
            first_op, _ = cigars[i][0]
            last_op, _ = cigars[i][-1]
            if "left_clippings" in features and first_op in (4,5):
                feature_dict["left_clippings"][ref_starts[i] % ref_length] += 1
            if "right_clippings" in features and last_op in (4,5):
                feature_dict["right_clippings"][(ref_ends[i]-1) % ref_length] += 1

        # --- Indels ---
        deletions = feature_dict["deletions"] if "deletions" in features else None
        insertions = feature_dict["insertions"] if "insertions" in features else None
        add_indels_numba(deletions, insertions, cigars[i], ref_starts[i], ref_length)

        # --- Mismatches ---
        if "mismatches" in features and has_md[i]:
            add_mismatches_numba_from_md(feature_dict["mismatches"], md_list[i], md_lengths[i], ref_starts[i], ref_length)

    # --- Finalize lengths ---
    if sequencing_type == "long":
        feature_dict["read_lengths"] = compute_final_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"])
    if sequencing_type == "short-paired":
        feature_dict["insert_sizes"] = compute_final_lengths(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"])

    # --- Save features into DB ---
    for feature in features:
        compress_and_write_features(feature, feature_dict[feature], db_conn, sample_id, contig_id, ref_length, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)

### Calculating features per contig per sample
@track_time
def preprocess_reads(reads_mapped, modules):
    ref_starts, ref_ends = [], []
    query_lengths, template_lengths = [], []
    is_read1, is_proper_pair, is_paired, is_reverse = [], [], [], []
    cigars, has_md, md_list, md_lengths = [], [], [], []

    need_md = {"phagetermini", "assemblycheck"}.intersection(modules)
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
        if need_md:
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
    # Keep `cigars` as a Python list of numpy arrays instead of an
    # object-dtype numpy array. Passing elements extracted from an
    # object array into numba-decorated functions causes numba to see
    # non-precise object types (array(pyobject,...)) and fail typing.
    # A list preserves the concrete numpy array elements when indexed.
    cigars = cigars
    has_md = np.array(has_md, dtype=np.bool_)
    md_lengths = np.array(md_lengths, dtype=np.int32)
    # Similarly, keep `md_list` as a list of numpy uint8 arrays so numba
    # functions that receive `md_list[i]` get a concrete array type.
    md_list = md_list

    return (ref_starts, ref_ends, query_lengths, template_lengths, is_read1, is_proper_pair, is_paired, is_reverse, cigars, has_md, md_list, md_lengths)

@track_time
def calculating_features_per_contig_per_sample(module_list, bam_file, sequencing_type, ref_name, locus_size, db_conn, sample_id=None, contig_id=None, 
                                               min_coverage=None, step=None, z_thresh=None, deriv_thresh=None, max_points=None, *, sample_name=None, contig_name=None, temp_mode=False):
    ### Save relevant info from bam file
    reads_mapped = bam_file.fetch(ref_name)
    (ref_starts, ref_ends, query_lengths, template_lengths,
     is_read1, is_proper_pair, is_paired, is_reverse,
     cigars, has_md, md_list, md_lengths) = preprocess_reads(reads_mapped, module_list)
    
    ### Coverage check: ensure more than min_coverage% of the reference is covered by at least one read
    if len(ref_starts) == 0:
        # If no reads fail coverage threshold
        return None
    
    # Mark covered regions
    # np.minimum / np.maximum keep indexes inside bounds
    covered = np.zeros(locus_size, dtype=bool)
    start_clipped = np.maximum(ref_starts, 0)
    end_clipped = np.minimum(ref_ends, locus_size)
    for s, e in zip(start_clipped, end_clipped):
        covered[s:e] = True

    covered_bp = covered.sum()
    coverage_pct = (covered_bp / locus_size) * 100

    if coverage_pct < min_coverage:
        # If not enough coverage skip this contig
        return None
    
    # Add presence record (either into main DB by id, or temp DB by names)
    if temp_mode:
        add_presence(db_conn, coverage_percentage=coverage_pct, contig_name=contig_name, sample_name=sample_name, temp_mode=True)
    else:
        add_presence(db_conn, contig_id=contig_id, sample_id=sample_id, coverage_percentage=coverage_pct)

    ### Calculate all features
    if {"coverage", "assemblycheck"}.intersection(module_list):
        get_feature_coverage(ref_starts, ref_ends, locus_size, db_conn, sample_id, contig_id, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)
    if "phagetermini" in module_list:
        get_features_phagetermini(ref_starts, ref_ends, is_reverse, cigars, md_list, sequencing_type, locus_size, 
                                  db_conn, sample_id, contig_id, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)
    if "assemblycheck" in module_list:
        get_features_assemblycheck(ref_starts, ref_ends, query_lengths, template_lengths, is_read1, is_proper_pair, cigars, has_md, md_list, md_lengths, 
                                   sequencing_type, locus_size, db_conn, sample_id, contig_id, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=contig_name, temp_mode=temp_mode)

### Distributing the calculation for all the contigs in the bam file
@track_time
def calculating_features_per_sample(module_list, mapping_file, db_conn, min_coverage, step, z_thresh, deriv_thresh, max_points):
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    references = bam_file.references
    lengths = [l // 2 for l in bam_file.lengths]  # need to divide by 2 because each contig was doubled for mapping to deal with circularity

    # Save sample and contig names into DB
    sample_name = os.path.basename(mapping_file).replace(".bam", "")
    sequencing_type = find_sequencing_type_from_bam(mapping_file)

    # Detect if db_conn is a temp DB by checking for TempFeatureValues table
    cur = db_conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='TempFeatureValues'")
    is_temp = cur.fetchone() is not None

    if is_temp:
        # temp DB: do not try to add sample/contig ids here. Write using names.
        for ref, length in zip(references, lengths):
            calculating_features_per_contig_per_sample(module_list, bam_file, sequencing_type, ref, length, db_conn, None, None,
                                                       min_coverage, step, z_thresh, deriv_thresh, max_points, sample_name=sample_name, contig_name=ref, temp_mode=True)
    else:
        sample_id = add_sample(db_conn, sample_name)
        for ref, length in zip(references, lengths):
            contig_id = get_contig_id(db_conn, ref)
            calculating_features_per_contig_per_sample(module_list, bam_file, sequencing_type, ref, length, db_conn, sample_id, contig_id, 
                                                       min_coverage, step, z_thresh, deriv_thresh, max_points)
    bam_file.close()

### Distributing the calculation for all the bam files (samples)
def _process_single_sample(list_modules, bam_file, db_path, min_coverage, step, z_thresh, deriv_thresh, max_points, n_contig_cores=None):
    # Backwards-compatible: `bam_file` may be a path string or a path-like object.
    bam_path = str(bam_file)

    # Reset timing counters for this sample so logs contain per-sample data only
    function_times.clear()

    # Prepare logfile path in same directory as the DB
    db_p = Path(db_path)
    log_path = db_p.parent / (db_p.stem + "_times.log")

    # Each worker may optionally parallelize over contigs. This function can be
    # called with an additional final argument `n_contig_cores` (int).
    # If caller passed a tuple with an extra element (when invoked via starmap
    # with extra arg), handle it gracefully. But in normal use callers should
    # pass n_contig_cores as a keyword in Python invocations.
    # Open the BAM to inspect contigs
    bam = pysam.AlignmentFile(bam_path, "rb")
    references = bam.references
    lengths = [l // 2 for l in bam.lengths]

    # Create a per-sample temp DB to avoid concurrent writers to the main DB
    db_main = db_path
    dbp = Path(db_main)
    sample_name = os.path.basename(bam_path).replace('.bam', '')
    temp_db = dbp.parent / f"{dbp.stem}_{sample_name}.temp.db"
    create_temp_sample_db(str(temp_db))

    # Open temp connection and run calculations (workers write into temp DB)
    conn_temp = sqlite3.connect(str(temp_db))
    # Copy Variable metadata from main DB into temp DB so compression can lookup Type
    try:
        with sqlite3.connect(db_path) as mconn:
            mcur = mconn.cursor()
            mcur.execute("SELECT Variable_name, Type FROM Variable")
            vars_meta = mcur.fetchall()
        tcur = conn_temp.cursor()
        tcur.execute("CREATE TABLE IF NOT EXISTS Variable (Variable_name TEXT PRIMARY KEY, Type TEXT, Feature_table_name TEXT)")
        if vars_meta:
            rows_to_insert = [(v[0], v[1]) for v in vars_meta]
            tcur.executemany("INSERT OR REPLACE INTO Variable (Variable_name, Type) VALUES (?, ?)", rows_to_insert)
            conn_temp.commit()
    except Exception:
        # If copying fails, continue; compression will fall back to defaults
        pass
    try:
        if not n_contig_cores or int(n_contig_cores) <= 1:
            # sequential: delegate to calculating_features_per_sample which will detect temp DB and write by names
            calculating_features_per_sample(list_modules, bam_path, conn_temp, min_coverage, step, z_thresh, deriv_thresh, max_points)
        else:
            # Parallelize contigs using a process Pool. Each worker will write into the shared temp DB (safe because each worker uses sqlite connection separately)
            from multiprocessing import Pool as MPPool

            def _process_single_contig_worker(args_tuple):
                (modules, bam_p, ref_name, length_val, temp_db_p, samp_name, min_cov, stp, zt, dt, mp, log_p) = args_tuple
                # Each worker keeps its own timing counters; after finishing the contig
                # append a brief timing snapshot to the shared logfile.
                conn_w = sqlite3.connect(temp_db_p)
                try:
                    bam_w = pysam.AlignmentFile(bam_p, "rb")
                    try:
                        seqtype = find_sequencing_type_from_bam(bam_p)
                        # Reset local function timing counters for this worker
                        function_times.clear()
                        calculating_features_per_contig_per_sample(modules, bam_w, seqtype, ref_name, length_val, conn_w, None, None, min_cov, stp, zt, dt, mp, sample_name=samp_name, contig_name=ref_name, temp_mode=True)
                        # Snapshot timings
                        try:
                            total_local = sum(function_times.values())
                            with open(log_p, 'a', encoding='utf8') as fh:
                                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tPID:{os.getpid()}\tSample:{samp_name}\tContig:{ref_name}\tTotal:{total_local:.3f}s\n")
                                for nm, tv in sorted(function_times.items(), key=lambda x: -x[1]):
                                    fh.write(f"\t{nm}: {tv:.3f}s\n")
                        except Exception:
                            pass
                    finally:
                        bam_w.close()
                finally:
                    conn_w.close()

            # Build tasks for contig workers
            tasks = []
            for ref, length in zip(references, lengths):
                tasks.append((list_modules, bam_path, ref, length, str(temp_db), sample_name, min_coverage, step, z_thresh, deriv_thresh, max_points, str(log_path)))

            with MPPool(processes=min(int(n_contig_cores), max(1, len(tasks)))) as pool:
                pool.map(_process_single_contig_worker, tasks)
    finally:
        conn_temp.close()

    # Merge temp DB into main DB and remove temp DB
    try:
        merge_temp_db_into_main(db_main, str(temp_db))
    except Exception as e:
        print(f"Warning: failed to merge temp DB {temp_db} into main DB {db_main}: {e}")
    try:
        Path(temp_db).unlink()
    except Exception:
        pass

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

def calculating_all_features_parallel(list_modules, bam_files, db_path, min_coverage, step, z_thresh, deriv_thresh, max_points, n_sample_cores=None):
    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)
    print(f"Using {n_sample_cores} cores to process {len(bam_files)} samples in parallel...", flush=True)

    args_list = [(list_modules, bam, db_path, min_coverage, step, z_thresh, deriv_thresh, max_points) for bam in bam_files]

    with Pool(processes=n_sample_cores) as pool:
        pool.starmap(_process_single_sample, args_list)
    print("Finished all samples.", flush=True)

### Main function helpers (shared-args)
def add_calculate_args(parser):
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file of all investigated contigs")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-d", "--db", required=True, help="Path to sqlite database file to store results")
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
    # If invoked in single-sample mode (used by Slurm jobs), simply process one bam
    if getattr(args, "run_sample", False):
        if not os.path.exists(args.db):
            sys.exit(f"ERROR: Database '{args.db}' does not exist for single-sample run.")
        requested_modules = args.modules.split(",")
        n_contig_cores = int(args.threads) if getattr(args, 'threads', None) else None
        _process_single_sample(requested_modules, args.bam_files, args.db, args.min_coverage, args.step, args.outlier_threshold, args.derivative_threshold, args.max_points, n_contig_cores)
        return

    # Sanity checks and preparation
    print("### Checking that genbank and mapping files are compatible...", flush=True)
    genbank_loci = set(rec.name for rec in SeqIO.parse(args.genbank, "genbank"))
    annotation_tool = args.annotation_tool
    if not genbank_loci:
        sys.exit("ERROR: No loci found in the provided GenBank file.")

    # Get list of BAM files
    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    # Check that references in BAM headers match loci in GenBank
    for bam_file in bam_files:
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                bam_refs = set(bam.references)
        except Exception as e:
            sys.exit(f"ERROR: Could not open BAM file '{bam_file}': {e}")

        missing_in_genbank = bam_refs - genbank_loci
        missing_in_bam = genbank_loci - bam_refs

        if missing_in_genbank:
            raise ValueError(
                f"ERROR: References in BAM file '{os.path.basename(bam_file)}' "
                f"not found in GenBank: '{os.path.basename(args.genbank)}'\n"
                f"{', '.join(sorted(missing_in_genbank))}"
            )

        if missing_in_bam:
            print(
                f"Warning: Some GenBank loci not present in BAM '{os.path.basename(bam_file)}': "
                f"{', '.join(sorted(missing_in_bam))}"
            )
    print("Sanity checks passed: No BAM contained unexpected reference names.", flush=True)

    # Requested modules
    requested_modules = args.modules.split(",")

    # Parameters for compression
    min_coverage = args.min_coverage
    step = args.step
    z_thresh = args.outlier_threshold
    deriv_thresh = args.derivative_threshold
    max_points = args.max_points

    n_cores = int(args.threads)

    # Running generate_database.py to create database
    db_path = args.db

    directory, filename = os.path.split(db_path)
    if not filename or not filename.endswith(".db"):
        sys.exit(f"ERROR: Invalid database path '{db_path}'. Must provide a file ending in '.db'.")
    if directory and not os.path.isdir(directory):
        sys.exit(f"ERROR: Directory '{directory}' does not exist. Please create it or choose a valid path.")
    if os.path.exists(db_path):
        sys.exit(f"ERROR: Database file '{db_path}' already exists. Please provide a new path to avoid overwriting.")

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
    # If requested, submit per-sample Slurm jobs instead of running locally
    if getattr(args, "use_slurm", False):
        if slurm_utils is None:
            sys.exit("ERROR: Slurm utilities not available in this installation.")

        outdir = Path(db_path).parent or Path(".")
        # Create a CSV listing bam files (one per row) so submit_array can iterate
        csv_path = outdir / (Path(db_path).stem + "_bam_list.csv")
        with open(csv_path, "w", newline="") as fh:
            for bam in bam_files:
                fh.write(f"{bam}\n")

        # Build the module CLI string to run a single sample. $bam will be substituted by the array script.
        # For Slurm runs we write per-sample temp DBs (to avoid contention) named by the bam basename.
        db_dir = outdir
        db_stem = Path(db_path).stem
        # use $readbase (set by the array script) to create per-sample db
        module_cli = (
            f"{sys.executable} -m mgfeatureviewer.calculating_data --run-sample --threads {args.threads_per_job} --bam_files $bam --db {db_dir}/{db_stem}_$readbase.temp.db "
            f"--modules {args.modules} --min_coverage {min_coverage} --step {step} --outlier_threshold {z_thresh} "
            f"--derivative_threshold {deriv_thresh} --max_points {max_points} {'--parallelize_contigs' if args.parallelize_contigs else ''}"
        )

        jobid = slurm_utils.submit_array(csv_path, outdir, module_cli, array_size=len(bam_files), concurrency=args.max_concurrent,
                                         cpus_per_task=args.threads, mem=args.mem_per_cpu, time=args.max_time, columns=["bam"])
        print(f"Submitted Slurm array job {jobid} to process {len(bam_files)} samples (csv: {csv_path})", flush=True)
        print("Note: Slurm tasks will write per-sample temp DBs named '<dbstem>_<readbase>.temp.db' in the DB directory. After jobs finish, run the provided merge helper to merge them into the main DB.")
        return

    # Otherwise run locally across samples
    calculating_all_features_parallel(requested_modules, bam_files, db_path, min_coverage, step, z_thresh, deriv_thresh, max_points, n_cores)

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