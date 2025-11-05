import numpy as np
import pysam
from re import findall
from multiprocessing import Pool, cpu_count
import constants
import os
import csv

# Function dedicated to specific features
def reduce_position(pos, ref_length):
    return pos % ref_length

def calculate_coverage(coverage, read, ref_length):
    for block_start, block_end in read.get_blocks():
        if block_start < ref_length and block_end > ref_length:  # circular boundary
            coverage[block_start:ref_length] += 1
            coverage[0:reduce_position(block_end, ref_length)] += 1
        else:
            coverage[reduce_position(block_start, ref_length):reduce_position(block_end, ref_length)] += 1

    return coverage

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

def calculate_reads_starts_and_ends(read, temporary_dict, ref_length):
    if starts_with_match(read, start=True) and starts_with_match(read, start=False):
        start = read.reference_start
        end = read.reference_end

        if start < ref_length and end > ref_length:
            temporary_dict["coverage_reduced"][start:ref_length] += 1
            temporary_dict["coverage_reduced"][0:reduce_position(end, ref_length)] += 1
        else:
            temporary_dict["coverage_reduced"][reduce_position(start, ref_length):reduce_position(end, ref_length)] += 1

        start_pos = reduce_position(read.reference_start, ref_length)
        end_pos = reduce_position(read.reference_end - 1, ref_length)  # end is exclusive

        if read.is_reverse:
            temporary_dict["start_minus"][start_pos] += 1
            temporary_dict["end_minus"][end_pos] += 1
        else:
            temporary_dict["start_plus"][start_pos] += 1
            temporary_dict["end_plus"][end_pos] += 1

    return temporary_dict

def compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length):
    # Combine per sequencing_type
    feature_dict["coverage_reduced"] = temporary_dict["coverage_reduced"]
    if sequencing_type == "short":
        feature_dict["reads_starts"] = temporary_dict["start_plus"]
        feature_dict["reads_ends"] = temporary_dict["end_minus"]
    elif sequencing_type == "long":
        feature_dict["reads_starts"] = temporary_dict["start_plus"] + temporary_dict["start_minus"]
        feature_dict["reads_ends"] = temporary_dict["end_plus"] + temporary_dict["end_minus"]
    else:
        raise ValueError("Sequencing_type must be 'short' or 'long'")

    # --- COMPUTE TAU ---
    tau = np.zeros(ref_length, dtype=np.float32)
    mask = feature_dict["coverage_reduced"] > 0
    tau[mask] = (feature_dict["reads_starts"][mask] + feature_dict["reads_ends"][mask]) / feature_dict["coverage_reduced"][mask]
    feature_dict["tau"] = tau

    return feature_dict

def compute_read_lengths(sum_read_lengths, counts_read_lengths, read, ref_length):
    read_length = read.query_length
    for pos in read.get_reference_positions():
        pos_reduced = reduce_position(pos, ref_length)
        sum_read_lengths[pos_reduced] += read_length
        counts_read_lengths[pos_reduced] += 1

    return sum_read_lengths, counts_read_lengths

def compute_inserts_characteristics(sum_insert_sizes, counts_insert_sizes, features, read, ref_length):
    arr_bad_orientations = np.zeros(ref_length, dtype=np.float32)
    for pos in read.get_reference_positions():
        pos_reduced = reduce_position(pos, ref_length)  # wrap-around safety

        if read.is_paired and not read.mate_is_unmapped:
            if "insert_sizes" in features and read.is_read1 and abs(read.template_length) > 0:
                sum_insert_sizes[pos_reduced] += abs(read.template_length)
                counts_insert_sizes[pos_reduced] += 1

            if "bad_orientations" in features and not read.is_proper_pair:
                arr_bad_orientations[pos_reduced] += 1
    
    return sum_insert_sizes, counts_insert_sizes, arr_bad_orientations

def compute_clippings(feature_dict, features, read, ref_length):
    if read.cigartuples:
        first_op, first_len = read.cigartuples[0]
        last_op, last_len = read.cigartuples[-1]

        if "left_clippings" in features:
            if first_op in (4, 5):
                pos_reduced = read.reference_start % ref_length
                feature_dict["left_clippings"][pos_reduced] += 1

        if "right_clippings" in features:
            if last_op in (4, 5):
                pos_reduced = (read.reference_end - 1) % ref_length
                feature_dict["right_clippings"][pos_reduced] += 1

    return feature_dict

def compute_indels(feature_dict, features, read, ref_length):
    ref_pos = read.reference_start
    for cigar_op, length in read.cigartuples:
        if "insertions" in features:  # Insertion
            if cigar_op == 1:
                pos_reduced = ref_pos % ref_length
                feature_dict["insertions"][pos_reduced] += 1
        elif "deletions" in features: # Deletion
            if cigar_op == 2:
                for i in range(length):
                    pos_reduced = (ref_pos + i) % ref_length
                    feature_dict["deletions"][pos_reduced] += 1
                ref_pos += length
        else:
            ref_pos += length

    return feature_dict

def compute_mismatches(feature_dict, features, read, ref_length):
    if "mismatches" in features:
        if read.has_tag("MD"):
            ref_pos = read.reference_start
            md_tag = read.get_tag("MD")
            tokens = findall(r'\d+|\^[A-Z]+|[A-Z]', md_tag)
            for token in tokens:
                if token.isdigit():
                    ref_pos += int(token)
                elif token.startswith("^"):
                    ref_pos += len(token) - 1
                else:
                    pos_reduced = ref_pos % ref_length
                    feature_dict["mismatches"][pos_reduced] += 1
                    ref_pos += 1
    
    return feature_dict

def compute_final_lengths(sum_lengths, count_lengths, ref_length):
    arr = np.zeros(ref_length, dtype=np.float32)
    for i in range(ref_length):
        total = sum_lengths[i]
        count = count_lengths[i]
        arr[i] = total / count if count > 0 else 0
    return arr

### Main logic to get features per position
def get_features(reads_mapped, features, ref_length, sequencing_type):    
    # Initialize arrays for all requested features
    feature_dict = {feature: np.zeros(ref_length, dtype=np.uint64) for feature in features}
    
    # Temporary array
    temporary_feature = []
    temporary_dict = {}
    if "tau" in features:
        temporary_feature.extend(["start_plus", "start_minus", "end_plus", "end_minus", "coverage_reduced"])
    if "read_lengths" in features:
        temporary_feature.extend(["sum_read_lengths", "count_read_lengths"])
    if "insert_sizes" in features:
        temporary_feature.extend(["sum_insert_sizes", "count_insert_sizes"])
    new_dict = {feature: np.zeros(ref_length, dtype=np.uint64) for feature in temporary_feature}
    temporary_dict.update(new_dict)

    for read in reads_mapped:
        if read.is_unmapped:
            continue

        # --- COVERAGE ---
        if "coverage" in features:
            feature_dict["coverage"] = calculate_coverage(feature_dict["coverage"], read, ref_length)

        # --- START/END BY STRAND ---
        if "tau" in features:
            temporary_dict = calculate_reads_starts_and_ends(read, temporary_dict, ref_length)

        # --- Read lengths / insert sizes / bad orientations ---
        if "read_lengths" in features:
            temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"] = compute_read_lengths(
                temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"], read, ref_length
            )
        if {"insert_sizes", "bad_orientations"}.intersection(features):
            temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"], feature_dict["bad_orientations"] = compute_inserts_characteristics(
                temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"], features, read, ref_length
            )

        # --- CIGAR parsing of the read endings ---
        if {"left_clippings", "right_clippings"}.intersection(features):
            feature_dict = compute_clippings(feature_dict, features, read, ref_length)

        # --- CIGAR parsing within the read for indels ---
        if {"insertions", "deletions"}.intersection(features):
            feature_dict = compute_indels(feature_dict, features, read, ref_length)

        # --- Mismatches can only obtained from MD tag ---
        if "mismatches" in features:
            feature_dict = compute_mismatches(feature_dict, features, read, ref_length)

    # --- FINALIZE START/END AND TAU ---
    if "tau" in features:
        feature_dict = compute_final_starts_ends_and_tau(feature_dict, temporary_dict, sequencing_type, ref_length)

    # --- FINALIZE READ LENGTHS AND INSERT SIZES ---
    if "read_lengths" in features:
        feature_dict["read_lengths"] = compute_final_lengths(temporary_dict["sum_read_lengths"], temporary_dict["count_read_lengths"], ref_length)
    if "insert_sizes" in features:
        feature_dict["insert_sizes"] = compute_final_lengths(temporary_dict["sum_insert_sizes"], temporary_dict["count_insert_sizes"], ref_length)

    return feature_dict

### Summarise data for each feature
def smooth_values(feature_values, ref_length, window_size):
    # Aggregate by window
    n_windows = (ref_length + window_size - 1) // window_size
    window_cov = np.add.reduceat(feature_values, np.arange(0, ref_length, window_size))
    window_cov = window_cov[:n_windows]

    # Compute actual counts per window (last one may be smaller)
    counts = np.full(n_windows, window_size)
    remainder = ref_length % window_size
    if remainder:
        counts[-1] = remainder

    # Convert to mean coverage per window safely
    coverage_list = (window_cov / counts).tolist()

    return coverage_list

def summarise_data(feature, feature_averaged, coverage_averaged, locus_size, window_size, min_ratio = 0.1):
    # Define window midpoints and y-values
    if window_size == 1:
        # For window size 1, use integer positions starting from 1
        xx = np.arange(1, locus_size + 1)
    else:
        # For larger windows, use midpoints including the last window if it fits
        xx = np.arange(window_size / 2, locus_size + window_size / 2, window_size)
        
    feature_averaged = np.asarray(feature_averaged, dtype=float)
    coverage_averaged = np.asarray(coverage_averaged, dtype=float)

    type_picked = constants.FEATURE_SUBPLOTS[feature]["type_picked"]
    if type_picked == "bars":
        with np.errstate(divide='ignore', invalid='ignore'):
            ratios = np.where(coverage_averaged > 0, feature_averaged / coverage_averaged, np.nan)

        # Mask values where ratio < min_ratio
        mask = ratios >= min_ratio
        xx = xx[mask]
        feature_averaged = feature_averaged[mask]

    if len(xx) == 0:
        return None  # skip this feature subplot
    return {"x": xx, "y": feature_averaged}

### Calculating features per contig per sample
def calculating_features_per_contig_per_sample(reads_mapped, locus_size, feature_list, sequencing_type, window_size):
    sequencing_type = "short" if sequencing_type.startswith("short") else "long"

    # Read bam once to calculate all features
    feature_values = get_features(reads_mapped, feature_list, locus_size, sequencing_type)

    # Averaging and masking data per feature
    coverage_averaged = smooth_values(feature_values["coverage"], locus_size, window_size)
    coverage_masked = summarise_data("coverage", coverage_averaged, coverage_averaged, locus_size, window_size)
    data = {"coverage": coverage_masked}  # always include coverage data

    for feature in feature_list:
        if feature != "coverage":
            feature_averaged = smooth_values(feature_values[feature], locus_size, window_size)
            feature_masked = summarise_data(feature, feature_averaged, coverage_averaged, locus_size, window_size)
            if feature_masked is not None:
                data[feature] = feature_masked
                
    return data
    
### Distributing the calculation for all the contigs in the bam file
def calculating_features_per_sample(mapping_file, feature_list, sequencing_type, window_size):
    print("Calculating features per sample...", flush=True)
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
    print(f"Processing sample: {sample_name}", flush=True)
    
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