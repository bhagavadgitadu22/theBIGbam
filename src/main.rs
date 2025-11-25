//! # MGFeatureViewer Rust Calculator
//!
//! This is a 1:1 Rust port of the Python calculation code in `mgfeatureviewer/calculating_data.py`.
//! It produces identical output values, with one minor difference noted below.
//!
//! ## Module Structure
//!
//! - `types.rs`      - All structs and enums
//! - `bam_reader.rs` - BAM file reading and sequencing type detection
//! - `features.rs`   - Feature calculations (coverage, phagetermini, assemblycheck)
//! - `compress.rs`   - Signal compression for storage
//! - `db.rs`         - SQLite database operations
//! - `parquet.rs`    - Parquet file writing
//! - `genbank.rs`    - GenBank file parsing
//! - `processing.rs` - Shared processing logic for CLI and Python bindings
//!
//! ## Function Mapping (Rust → Python)
//!
//! | Rust Function                    | Python Function                              | Line  |
//! |----------------------------------|----------------------------------------------|-------|
//! | `detect_sequencing_type`         | `find_sequencing_type_from_bam`              | ~55   |
//! | `process_reads_for_contig`       | `preprocess_reads`                           | ~497  |
//! | `calculate_coverage`             | `calculate_coverage_numba`                   | ~182  |
//! | `calculate_phagetermini`         | `get_features_phagetermini`                  | ~283  |
//! | `calculate_assemblycheck`        | `get_features_assemblycheck`                 | ~410  |
//! | `compress_signal`                | `compress_signal`                            | ~91   |
//! | `process_sample`                 | `calculating_features_per_sample`            | ~651  |
//!
//! ## Known Difference: Derivative Outlier Detection
//!
//! The Python `compress_signal` function (line ~116) has a subtle bug in derivative outlier
//! detection where it performs arithmetic on a boolean array. This causes Python to keep
//! fewer positions after compression. The Rust implementation correctly identifies derivative
//! outliers, so it may keep slightly more positions. **All values at common positions are
//! identical** - Rust just preserves more data points as originally intended.
//!
//! ## Verification
//!
//! Run `python compare_values_only.py <sqlite_db> <parquet_dir>` to verify outputs match:
//! - All 11 features tested across all samples
//! - 100% exact match at common positions (0 difference)
//! - Metadata tables match exactly
//!
//! See also: `PYTHON_COMPARISON.md` for side-by-side code examples.

use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

use mgfeatureviewer_rs::{run_all_samples, ProcessConfig};

#[derive(Parser, Debug)]
#[command(name = "mgfeatureviewer_calc")]
#[command(about = "Calculate features from BAM files (Rust implementation)")]
struct Args {
    /// Number of threads
    #[arg(short = 't', long)]
    threads: usize,

    /// Path to GenBank file
    #[arg(short = 'g', long)]
    genbank: PathBuf,

    /// Path to BAM file or directory
    #[arg(short = 'b', long)]
    bam_files: PathBuf,

    /// Modules to compute (comma-separated: coverage,phagetermini,assemblycheck)
    #[arg(short = 'm', long)]
    modules: String,

    /// Output directory
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Annotation tool (optional)
    #[arg(short = 'a', long, default_value = "")]
    annotation_tool: String,

    /// Minimum coverage percentage for contig inclusion
    #[arg(long, default_value = "50")]
    min_coverage: f64,

    /// Step size for compression
    #[arg(long, default_value = "50")]
    step: usize,

    /// Z-score threshold for outliers
    #[arg(long, default_value = "3.0")]
    outlier_threshold: f64,

    /// Derivative threshold for outliers
    #[arg(long, default_value = "3.0")]
    derivative_threshold: f64,

    /// Maximum points after compression
    #[arg(long, default_value = "10000")]
    max_points: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Check output directory doesn't exist
    if args.output.exists() {
        anyhow::bail!("Output directory already exists: {:?}", args.output);
    }

    // Parse modules
    let modules: Vec<String> = args.modules.split(',').map(|s| s.trim().to_string()).collect();

    let config = ProcessConfig {
        threads: args.threads,
        min_coverage: args.min_coverage,
        step: args.step,
        z_thresh: args.outlier_threshold,
        deriv_thresh: args.derivative_threshold,
        max_points: args.max_points,
    };

    // Call the shared processing function
    let result = run_all_samples(
        &args.genbank,
        &args.bam_files,
        &args.output,
        &modules,
        &args.annotation_tool,
        &config,
    )?;

    if result.samples_failed > 0 {
        eprintln!("\nWarning: {} samples failed to process", result.samples_failed);
    }

    Ok(())
}
