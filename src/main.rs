//! # MGFeatureViewer Rust Calculator - DEPRECATED
//!
//! **⚠️ This standalone Rust CLI is deprecated and no longer built.**
//!
//! Use the Python CLI instead:
//!     mgfeatureviewer calculate -g assembly.gbk -b bams/ -m coverage -o output.db
//!
//! This file is kept for reference/debugging but is not compiled into a binary.
//! The Python CLI (`mgfeatureviewer calculate`) calls the same Rust code via PyO3 bindings.
//!
//! ## Why deprecated?
//! - Duplicate argument parsing (this file vs. calculating_data.py)
//! - Maintenance burden - every parameter change needed updating in 2 places
//! - Python CLI provides unified interface for all commands (calculate, plot, serve, etc.)
//!
//! ## For developers
//! If you need to test Rust code directly without Python, you can temporarily re-enable
//! the binary in Cargo.toml by uncommenting the [[bin]] section.

#![allow(dead_code)]

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

    /// Output database file (.db)
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Annotation tool (optional)
    #[arg(short = 'a', long, default_value = "")]
    annotation_tool: String,

    /// Minimum alignment-length coverage proportion for contig inclusion (default 50%% change threshold)
    #[arg(long, default_value = "50")]
    min_coverage: f64,

    /// RLE compression ratio (default 10%% change threshold)
    #[arg(long, default_value = "10")]
    compress_ratio: f64,

    /// Circular genome flag: set if assembly was doubled during mapping (enables modulo logic)
    #[arg(long, default_value = "false")]
    circular: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Check output file doesn't exist
    if args.output.exists() {
        anyhow::bail!("Output file already exists: {:?}", args.output);
    }

    // Parse modules
    let modules: Vec<String> = args.modules.split(',').map(|s| s.trim().to_string()).collect();

    let config = ProcessConfig {
        threads: args.threads,
        min_coverage: args.min_coverage,
        compress_ratio: args.compress_ratio,
        circular: args.circular,
    };

    // Call the shared processing function
    let result = run_all_samples(
        &args.genbank,
        &args.bam_files,
        &args.output,
        &modules,
        &args.annotation_tool,
        &config,
        true,  // create_indexes
    )?;

    if result.samples_failed > 0 {
        eprintln!("\nWarning: {} samples failed to process", result.samples_failed);
    }

    Ok(())
}
