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
pub mod compress;
pub mod db;
pub mod features;
pub mod genbank;
pub mod processing;
pub mod types;

// Re-export main types for library users
pub use bam_reader::{detect_sequencing_type, process_reads_for_contig};
pub use compress::compress_signal;
pub use features::{calculate_assemblycheck, calculate_coverage, calculate_phagetermini};
pub use genbank::parse_genbank;
pub use processing::{process_sample, run_all_samples, ProcessConfig, ProcessResult};
pub use types::{ContigInfo, FeaturePoint, PlotType, PresenceData, ReadData, SequencingType};

// ============================================================================
// PyO3 Python Bindings
// ============================================================================

#[cfg(feature = "python")]
mod python {
    use super::*;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;
    use std::path::Path;

    /// Process all BAM files in parallel with rayon and progress bar.
    ///
    /// This is the main API for Python. It parses GenBank once, then processes
    /// all BAM files in parallel using rayon.
    ///
    /// Args:
    ///     genbank_path: Path to the GenBank annotation file
    ///     bam_dir: Path to directory containing BAM files (or single BAM file)
    ///     output_db: Output database file path (.db)
    ///     modules: List of modules to compute: "coverage", "phagetermini", "assemblycheck"
    ///     threads: Number of threads to use
    ///     annotation_tool: Annotation tool name (e.g., "pharokka")
    ///     min_coverage: Minimum coverage percentage for contig inclusion (default 50.0)
    ///     step: Compression step size (default 50)
    ///     z_thresh: Z-score threshold for outliers (default 3.0)
    ///     deriv_thresh: Derivative threshold for outliers (default 3.0)
    ///     max_points: Maximum points per feature (default 10000)
    ///
    /// Returns:
    ///     dict with keys:
    ///         - "samples_processed": int
    ///         - "samples_failed": int
    ///         - "total_time": float (seconds)
    #[pyfunction]
    #[pyo3(signature = (genbank_path, bam_dir, output_db, modules, threads, annotation_tool="", min_coverage=50.0, step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000, python_compat=false))]
    fn process_all_samples<'py>(
        py: Python<'py>,
        genbank_path: &str,
        bam_dir: &str,
        output_db: &str,
        modules: Vec<String>,
        threads: usize,
        annotation_tool: &str,
        min_coverage: f64,
        step: usize,
        z_thresh: f64,
        deriv_thresh: f64,
        max_points: usize,
        python_compat: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        use crate::processing::{run_all_samples, ProcessConfig};

        let config = ProcessConfig {
            threads,
            min_coverage,
            step,
            z_thresh,
            deriv_thresh,
            max_points,
            python_compat,
        };

        // Release the GIL and call the shared processing function
        let process_result = py.allow_threads(|| {
            run_all_samples(
                Path::new(genbank_path),
                Path::new(bam_dir),
                Path::new(output_db),
                &modules,
                annotation_tool,
                &config,
            )
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{}", e)))
        })?;

        let result = PyDict::new(py);
        result.set_item("samples_processed", process_result.samples_processed)?;
        result.set_item("samples_failed", process_result.samples_failed)?;
        result.set_item("total_time", process_result.total_time_secs)?;
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
