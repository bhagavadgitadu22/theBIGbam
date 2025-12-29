//! MGFeatureViewer Rust Library
//!
//! This module provides both a Rust library and Python bindings (via PyO3).
//!
//! ## Python Usage
//!
//! ```python
//! import mgfeatureviewer_rs as mgfv
//!
//! # Process all BAM files in parallel (main API)
//! result = mgfv.process_all_samples(
//!     genbank_path="annotation.gbk",
//!     bam_dir="/path/to/bams",
//!     output_dir="/path/to/output",
//!     modules=["coverage", "phagetermini", "assemblycheck"],
//!     threads=8,
//!     annotation_tool="pharokka",
//! )
//! # result = {"samples_processed": 10, "samples_failed": 0, "total_time": 123.4}
//! ```

pub mod bam_reader;
pub mod cigar;
pub mod circular;
pub mod compress;
pub mod db;
pub mod features;
pub mod genbank;
pub mod processing;
pub mod processing_completeness;
pub mod processing_phage_packaging;
pub mod types;

// Re-export main types for library users
pub use bam_reader::{detect_sequencing_type, process_contig_streaming};
pub use cigar::{Cigar, CigarElement, CigarOp, MdTag};
pub use circular::{increment_circular, increment_range};
pub use compress::compress_signal_with_reference;
pub use features::{process_read, FeatureArrays, ModuleFlags};
pub use genbank::parse_genbank;
pub use processing::{process_sample, run_all_samples, ProcessConfig, ProcessResult};
pub use processing_phage_packaging::PhageTerminiConfig;
pub use types::{
    mean_std, ContigInfo, FeaturePoint, PlotType, PresenceData, SequencingType,
    ASSEMBLYCHECK_FEATURES, PHAGETERMINI_FEATURES, VARIABLES,
};

// ============================================================================
// PyO3 Python Bindings
// ============================================================================

#[cfg(feature = "python")]
mod python {
    use pyo3::prelude::*;
    use pyo3::types::PyDict;
    use std::path::Path;

    /// Process all BAM files in parallel with rayon and progress bar.
    ///
    /// This is the main API for Python. It parses GenBank once, then processes
    /// all BAM files in parallel using rayon.
    ///
    /// Args:
    ///     genbank_path: Path to the GenBank annotation file (empty string to skip)
    ///     bam_files: List of BAM file paths to process
    ///     output_db: Output database file path (.db)
    ///     modules: List of modules to compute: "coverage", "phagetermini", "assemblycheck"
    ///     threads: Number of threads to use
    ///     annotation_tool: Annotation tool name (e.g., "pharokka")
    ///     min_coverage: Minimum coverage percentage for contig inclusion (default 50.0)
    ///     compress_ratio: Compression ratio threshold (default 10.0)
    ///     circular: Whether the genome is circular (default False)
    ///     create_indexes: Whether to create database indexes (default True)
    ///
    /// Returns:
    ///     dict with keys:
    ///         - "samples_processed": int
    ///         - "samples_failed": int
    ///         - "total_time": float (seconds)
    #[pyfunction]
    #[pyo3(signature = (genbank_path, bam_files, output_db, modules, threads, sequencing_type=None, annotation_tool="", min_coverage=50.0, curve_ratio=10.0, bar_ratio=10.0, circular=false, create_indexes=true, max_samples_in_memory=10, autoblast_file=""))]
    fn process_all_samples<'py>(
        py: Python<'py>,
        genbank_path: &str,
        bam_files: Vec<String>,
        output_db: &str,
        modules: Vec<String>,
        threads: usize,
        sequencing_type: Option<&str>,
        annotation_tool: &str,
        min_coverage: f64,
        curve_ratio: f64,
        bar_ratio: f64,
        circular: bool,
        create_indexes: bool,
        max_samples_in_memory: usize,
        autoblast_file: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        use crate::processing::{run_all_samples, ProcessConfig};
        use crate::processing_phage_packaging::PhageTerminiConfig;
        use std::path::PathBuf;

        // Parse sequencing type: Some(type) if user provided valid value, None for per-sample auto-detection
        let seq_type = sequencing_type.and_then(|s| ProcessConfig::parse_sequencing_type(s));

        let config = ProcessConfig {
            threads,
            min_coverage,
            curve_ratio,
            bar_ratio,
            circular,
            sequencing_type: seq_type,
            phagetermini_config: PhageTerminiConfig::default(),
        };

        // Convert string paths to PathBuf
        let bam_paths: Vec<PathBuf> = bam_files.iter().map(|s| PathBuf::from(s)).collect();

        // Release the GIL and call the shared processing function
        let process_result = py.allow_threads(|| {
            run_all_samples(
                Path::new(genbank_path),
                &bam_paths,
                Path::new(output_db),
                &modules,
                annotation_tool,
                &config,
                create_indexes,
                max_samples_in_memory,
                Path::new(autoblast_file),
            )
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{}", e)))
        })?;

        let result = PyDict::new(py);
        result.set_item("samples_processed", process_result.samples_processed)?;
        result.set_item("samples_failed", process_result.samples_failed)?;
        result.set_item("total_time", process_result.total_time_secs)?;
        result.set_item("processing_time", process_result.processing_time_secs)?;
        result.set_item("writing_time", process_result.writing_time_secs)?;
        Ok(result)
    }

    /// MGFeatureViewer Rust bindings for Python.
    #[pymodule]
    fn mgfeatureviewer_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(process_all_samples, m)?)?;
        m.add("__doc__", "MGFeatureViewer Rust bindings - fast BAM processing")?;
        Ok(())
    }
}
