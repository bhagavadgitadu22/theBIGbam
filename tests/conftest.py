"""
Pytest configuration and fixtures for theBIGbam tests.

This module provides fixtures for generating test BAM files from FASTQ data.
BAMs are generated once and cached in the tests/ directory.
"""

import os
import subprocess
import pytest

# Test data directory
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(TESTS_DIR, "test_10000bp.fasta")

# FASTQ files and their mapper presets
# Format: (r1_file, r2_file or None, mapper)
FASTQ_FILES = [
    # Paired-end short reads
    ("50_read_pairs_for_test_10kbp_R1.fastq", "50_read_pairs_for_test_10kbp_R2.fastq", "minimap2-sr-secondary"),
    ("50_read_pairs_for_test_10kbp_inverted_R1.fastq", "50_read_pairs_for_test_10kbp_inverted_R2.fastq", "minimap2-sr-secondary"),
    ("5000_read_pairs_for_test_10kbp_concatenated_100_times_R1.fastq", "5000_read_pairs_for_test_10kbp_concatenated_100_times_R2.fastq", "minimap2-sr-secondary"),
    # Long reads (single file)
    ("1000_long_reads_for_test_10kbp.fastq", None, "minimap2-ont"),
    ("100_long_reads_for_test_10kbp_concatenated_100_times.fastq", None, "minimap2-ont"),
]

def get_bam_name(fastq_name: str, circular: bool) -> str:
    """Generate BAM filename from FASTQ name."""
    base = fastq_name.replace(".fastq", "").replace(".fq", "")
    # Remove R1/R2 suffix for paired reads
    base = base.replace("_R1", "").replace("_R2", "")
    # Shorten long names
    base = base.replace("_for_test_10kbp", "").replace("_concatenated_100_times", "_concat")
    suffix = "_circular" if circular else ""
    return f"{base}{suffix}.bam"

def generate_bam(r1_path: str, r2_path: str | None, mapper: str, circular: bool, output_bam: str) -> None:
    """Generate a BAM file using thebigbam mapping-per-sample command."""
    cmd = [
        "thebigbam", "mapping-per-sample",
        "-r1", r1_path,
        "-a", REFERENCE,
        "--mapper", mapper,
        "-o", output_bam,
    ]
    if r2_path:
        cmd.extend(["-r2", r2_path])
    if circular:
        cmd.append("--circular")

    print(f"Generating BAM: {os.path.basename(output_bam)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"BAM generation failed: {result.stderr}")

@pytest.fixture(scope="session")
def test_bams():
    """
    Fixture that ensures all test BAM files exist.

    Generates BAMs from FASTQ files if they don't exist.
    BAMs are stored in linear_bams/ and circular_bams/ subdirectories.
    Returns a dict with 'circular' and 'linear' BAM file lists.
    """
    # Create output directories
    linear_dir = os.path.join(TESTS_DIR, "linear_bams")
    circular_dir = os.path.join(TESTS_DIR, "circular_bams")
    os.makedirs(linear_dir, exist_ok=True)
    os.makedirs(circular_dir, exist_ok=True)

    bam_files = {"circular": [], "linear": []}

    for r1_name, r2_name, mapper in FASTQ_FILES:
        r1_path = os.path.join(TESTS_DIR, r1_name)
        r2_path = os.path.join(TESTS_DIR, r2_name) if r2_name else None

        if not os.path.exists(r1_path):
            pytest.skip(f"FASTQ file not found: {r1_path}")
        if r2_path and not os.path.exists(r2_path):
            pytest.skip(f"FASTQ file not found: {r2_path}")

        # Generate both linear and circular BAMs
        for circular in [False, True]:
            bam_name = get_bam_name(r1_name, circular)
            out_dir = circular_dir if circular else linear_dir
            bam_path = os.path.join(out_dir, bam_name)
            bai_path = bam_path + ".bai"

            # Generate if BAM or index doesn't exist
            if not os.path.exists(bam_path) or not os.path.exists(bai_path):
                generate_bam(r1_path, r2_path, mapper, circular, bam_path)

            if os.path.exists(bam_path):
                key = "circular" if circular else "linear"
                bam_files[key].append(bam_path)

    return bam_files

@pytest.fixture(scope="session")
def tests_dir():
    """Return the tests directory path."""
    return TESTS_DIR

@pytest.fixture(scope="session")
def reference_fasta():
    """Return the reference FASTA path."""
    return REFERENCE
