"""
Pytest configuration and fixtures for theBIGbam tests.

This module provides fixtures for generating test BAM files from FASTQ data.
BAMs are generated once and cached in the tests/ directory.
"""

import os
import shutil
import subprocess
import sys

import pytest

# Test data directory
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(TESTS_DIR, "fasta", "test_10000bp.fasta")


def pytest_collection_modifyitems(items):
    """Classify every test into one explicit execution tier."""
    for item in items:
        if item.get_closest_marker("integration") or item.get_closest_marker("browser"):
            continue
        item.add_marker(pytest.mark.fast)


def thebigbam_command(*args: str) -> list[str]:
    """Build a CLI command from the environment running pytest.

    Restricting the lookup to the active interpreter's scripts directory
    prevents a different Conda environment's entry point from leaking in
    through PATH.
    """
    scripts_dir = os.path.dirname(sys.executable)
    executable = shutil.which("thebigbam", path=scripts_dir)
    if executable is None:
        pytest.fail(
            f"thebigbam is not installed in the active test environment ({scripts_dir}); "
            "install the project there with `python -m pip install -e .`"
        )
    return [executable, *args]


# FASTQ files and their mapper presets
# Format: (r1_file, r2_file or None, mapper)
FASTQ_FILES = [
    # Paired-end short reads
    ("50_read_pairs_for_test_10kbp_R1.fastq", "50_read_pairs_for_test_10kbp_R2.fastq", "minimap2-sr-secondary"),
    (
        "50_read_pairs_for_test_10kbp_inverted_R1.fastq",
        "50_read_pairs_for_test_10kbp_inverted_R2.fastq",
        "minimap2-sr-secondary",
    ),
    (
        "5000_read_pairs_for_test_10kbp_concatenated_100_times_R1.fastq",
        "5000_read_pairs_for_test_10kbp_concatenated_100_times_R2.fastq",
        "minimap2-sr-secondary",
    ),
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
    cmd = thebigbam_command(
        "mapping-per-sample",
        "-r1",
        r1_path,
        "-a",
        REFERENCE,
        "--mapper",
        mapper,
        "-o",
        output_bam,
    )
    if r2_path:
        cmd.extend(["-r2", r2_path])
    if circular:
        cmd.append("--circular")

    print(f"Generating BAM: {os.path.basename(output_bam)}", flush=True)
    print(f"CMD: {cmd}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"BAM generation failed: {result.stderr}")


@pytest.fixture(scope="session")
def test_bams(tmp_path_factory):
    """
    Fixture that ensures all test BAM files exist.

    Generates BAMs from FASTQ files if they don't exist.
    BAMs are stored in a session-temporary build directory.
    Returns a dict with 'circular' and 'linear' BAM file lists.
    """
    # Create output directories
    build_dir = str(tmp_path_factory.mktemp("mapping-artifacts"))
    linear_dir = os.path.join(build_dir, "linear_bams")
    circular_dir = os.path.join(build_dir, "circular_bams")
    os.makedirs(linear_dir, exist_ok=True)
    os.makedirs(circular_dir, exist_ok=True)

    bam_files = {"circular": [], "linear": []}

    for r1_name, r2_name, mapper in FASTQ_FILES:
        r1_path = os.path.join(TESTS_DIR, "fastq", r1_name)
        r2_path = os.path.join(TESTS_DIR, "fastq", r2_name) if r2_name else None

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

    return {**bam_files, "linear_dir": linear_dir, "circular_dir": circular_dir}


@pytest.fixture(scope="session")
def tests_dir():
    """Return the tests directory path."""
    return TESTS_DIR


@pytest.fixture(scope="session")
def reference_fasta():
    """Return the reference FASTA path."""
    return REFERENCE


@pytest.fixture(scope="session")
def linear_db(tmp_path_factory, test_bams):
    """Create linear database from linear BAMs. Returns DB path."""
    db_path = str(tmp_path_factory.mktemp("databases") / "test_10kbp_linear.db")
    linear_dir = test_bams["linear_dir"]
    cmd = thebigbam_command(
        "calculate",
        "-b",
        linear_dir,
        "-a",
        REFERENCE,
        "-m",
        "coverage,phagetermini",
        "-o",
        db_path,
        "-t",
        "4",
        "--coverage_percentage",
        "0",
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"Calculate failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
    return db_path


@pytest.fixture(scope="session")
def circular_db(tmp_path_factory, test_bams):
    """Create circular database from circular BAMs. Returns DB path."""
    db_path = str(tmp_path_factory.mktemp("databases") / "test_10kbp_circular.db")
    circular_dir = test_bams["circular_dir"]
    cmd = thebigbam_command(
        "calculate",
        "-b",
        circular_dir,
        "-a",
        REFERENCE,
        "-m",
        "coverage,phagetermini",
        "-o",
        db_path,
        "-t",
        "4",
        "--coverage_percentage",
        "0",
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"Calculate failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
    return db_path


@pytest.fixture(scope="session")
def linear_db_simple(tmp_path_factory, test_bams):
    """Create linear database without percentage options. Returns DB path."""
    db_path = str(tmp_path_factory.mktemp("databases") / "test_10kbp_linear_simple.db")
    linear_dir = test_bams["linear_dir"]
    cmd = thebigbam_command(
        "calculate",
        "-b",
        linear_dir,
        "-a",
        REFERENCE,
        "-m",
        "coverage,phagetermini",
        "-o",
        db_path,
        "-t",
        "4",
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"Calculate failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
    return db_path


@pytest.fixture(scope="session")
def circular_db_simple(tmp_path_factory, test_bams):
    """Create circular database without percentage options. Returns DB path."""
    db_path = str(tmp_path_factory.mktemp("databases") / "test_10kbp_circular_simple.db")
    circular_dir = test_bams["circular_dir"]
    cmd = thebigbam_command(
        "calculate",
        "-b",
        circular_dir,
        "-a",
        REFERENCE,
        "-m",
        "coverage,phagetermini",
        "-o",
        db_path,
        "-t",
        "4",
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"Calculate failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
    return db_path
