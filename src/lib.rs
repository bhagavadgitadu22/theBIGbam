//! MGFeatureViewer Rust Library
//!
//! This module provides both a Rust library and Python bindings (via PyO3).
//!
//! ## Python Usage
//!
//! ```python
//! import mgfeatureviewer_rust as mgfv
//!
//! # Detect sequencing type from BAM
//! seq_type = mgfv.detect_sequencing_type("sample.bam")  # "long", "short-paired", or "short-single"
//!
//! # Process entire BAM file - returns dict of compressed features
//! result = mgfv.process_bam_file(
//!     bam_path="sample.bam",
//!     genbank_path="annotation.gbk",
//!     modules=["coverage", "phagetermini", "assemblycheck"],
//!     annotation_tool="pharokka",
//!     min_coverage=50.0,
//!     step=50,
//!     z_thresh=3.0,
//!     deriv_thresh=3.0,
//!     max_points=10000,
//! )
//! # result = {
//! #     "features": [(contig, feature, position, value), ...],
//! #     "presences": [(contig, coverage_pct), ...],
//! #     "sample_name": "sample",
//! # }
//!
//! # Compress a signal (works with numpy arrays)
//! import numpy as np
//! values = np.array([1.0, 2.0, 3.0, ...])
//! x, y = mgfv.compress_signal(values, "curve", step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000)
//! ```

pub mod bam_reader;
pub mod compress;
pub mod db;
pub mod features;
pub mod genbank;
pub mod parquet;
pub mod processing;
pub mod types;

// Re-export main types for library users
pub use bam_reader::{detect_sequencing_type, process_reads_for_contig};
pub use compress::compress_signal;
pub use features::{calculate_assemblycheck, calculate_coverage, calculate_phagetermini};
pub use genbank::parse_genbank;
pub use processing::{run_all_samples, process_sample, ProcessConfig, ProcessResult};
pub use types::{ContigInfo, FeaturePoint, PlotType, PresenceData, ReadData, SequencingType};

// ============================================================================
// PyO3 Python Bindings
// ============================================================================

#[cfg(feature = "python")]
mod python {
    use super::*;
    use numpy::{PyArray1, PyArrayMethods};
    use pyo3::prelude::*;
    use pyo3::types::{PyDict, PyList};
    use rust_htslib::bam::{IndexedReader, Read as BamRead};
    use std::path::Path;

    /// Detect sequencing type from a BAM file.
    ///
    /// Args:
    ///     bam_path: Path to the BAM file
    ///
    /// Returns:
    ///     str: "long", "short-paired", or "short-single"
    #[pyfunction]
    fn detect_sequencing_type_py(bam_path: &str) -> PyResult<String> {
        let seq_type = detect_sequencing_type(Path::new(bam_path))?;
        Ok(match seq_type {
            SequencingType::Long => "long".to_string(),
            SequencingType::ShortPaired => "short-paired".to_string(),
            SequencingType::ShortSingle => "short-single".to_string(),
        })
    }

    /// Compress a signal for efficient storage.
    ///
    /// Args:
    ///     values: numpy array of feature values
    ///     plot_type: "curve" or "bars"
    ///     step: Keep every Nth point (default 50)
    ///     z_thresh: Z-score threshold for outliers (default 3.0)
    ///     deriv_thresh: Derivative threshold for outliers (default 3.0)
    ///     max_points: Maximum points to keep (default 10000)
    ///
    /// Returns:
    ///     tuple: (x_positions, y_values) as numpy arrays
    #[pyfunction]
    #[pyo3(signature = (values, plot_type, step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000))]
    fn compress_signal_py<'py>(
        py: Python<'py>,
        values: &Bound<'py, PyArray1<f64>>,
        plot_type: &str,
        step: usize,
        z_thresh: f64,
        deriv_thresh: f64,
        max_points: usize,
    ) -> PyResult<(Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<f32>>)> {
        let values_vec: Vec<f64> = unsafe { values.as_slice()?.to_vec() };

        let pt = match plot_type {
            "curve" => PlotType::Curve,
            "bars" => PlotType::Bars,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "plot_type must be 'curve' or 'bars'",
                ))
            }
        };

        let (xs, ys) = compress_signal(&values_vec, pt, step, z_thresh, deriv_thresh, max_points);

        let x_arr = PyArray1::from_vec(py, xs);
        let y_arr = PyArray1::from_vec(py, ys);

        Ok((x_arr, y_arr))
    }

    /// Calculate coverage from preprocessed read data.
    ///
    /// This is a lower-level function. For most use cases, prefer `process_bam_file`.
    ///
    /// Args:
    ///     ref_starts: numpy array of read start positions
    ///     ref_ends: numpy array of read end positions
    ///     ref_length: Length of the reference contig
    ///
    /// Returns:
    ///     numpy array of coverage values per position
    #[pyfunction]
    fn calculate_coverage_py<'py>(
        py: Python<'py>,
        ref_starts: &Bound<'py, PyArray1<i64>>,
        ref_ends: &Bound<'py, PyArray1<i64>>,
        ref_length: usize,
    ) -> PyResult<Bound<'py, PyArray1<u64>>> {
        let starts: Vec<i64> = unsafe { ref_starts.as_slice()?.to_vec() };
        let ends: Vec<i64> = unsafe { ref_ends.as_slice()?.to_vec() };

        // Build minimal ReadData structs for coverage calculation
        let reads: Vec<ReadData> = starts
            .into_iter()
            .zip(ends)
            .map(|(s, e)| ReadData {
                ref_start: s,
                ref_end: e,
                query_length: 0,
                template_length: 0,
                is_read1: false,
                is_proper_pair: false,
                is_reverse: false,
                cigar: Vec::new(),
                md_tag: None,
            })
            .collect();

        let coverage = calculate_coverage(&reads, ref_length);
        Ok(PyArray1::from_vec(py, coverage))
    }

    /// Process a single BAM file and return all computed features.
    ///
    /// This is the main high-level API for feature computation.
    ///
    /// Args:
    ///     bam_path: Path to the BAM file (must be indexed)
    ///     genbank_path: Path to the GenBank annotation file
    ///     modules: List of modules to compute: "coverage", "phagetermini", "assemblycheck"
    ///     annotation_tool: Annotation tool name (e.g., "pharokka")
    ///     min_coverage: Minimum coverage percentage for contig inclusion (default 50.0)
    ///     step: Compression step size (default 50)
    ///     z_thresh: Z-score threshold for outliers (default 3.0)
    ///     deriv_thresh: Derivative threshold for outliers (default 3.0)
    ///     max_points: Maximum points per feature (default 10000)
    ///
    /// Returns:
    ///     dict with keys:
    ///         - "features": list of (contig_name, feature_name, position, value) tuples
    ///         - "presences": list of (contig_name, coverage_pct) tuples
    ///         - "sample_name": str
    ///         - "contigs": list of {"name": str, "length": int} dicts
    #[pyfunction]
    #[pyo3(signature = (bam_path, genbank_path, modules, annotation_tool="", min_coverage=50.0, step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000))]
    fn process_bam_file<'py>(
        py: Python<'py>,
        bam_path: &str,
        genbank_path: &str,
        modules: Vec<String>,
        annotation_tool: &str,
        min_coverage: f64,
        step: usize,
        z_thresh: f64,
        deriv_thresh: f64,
        max_points: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        // Parse GenBank file
        let (contigs, _annotations) = parse_genbank(Path::new(genbank_path), annotation_tool)?;

        // Extract sample name from BAM path
        let bam_path_obj = Path::new(bam_path);
        let sample_name = bam_path_obj
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .replace("_with_MD", "");

        // Detect sequencing type
        let seq_type = detect_sequencing_type(bam_path_obj)?;

        // Open indexed BAM
        let mut bam = IndexedReader::from_path(bam_path_obj)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to open BAM file: {}", e)))?;
        bam.set_threads(2)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to set threads: {}", e)))?;

        let mut all_features: Vec<(String, String, i32, f32)> = Vec::new();
        let mut all_presences: Vec<(String, f32)> = Vec::new();

        // Process each contig
        for contig in &contigs {
            // Get reference length from BAM header
            let tid = bam.header().tid(contig.name.as_bytes());
            let ref_length = if let Some(tid) = tid {
                let bam_length = bam.header().target_len(tid).unwrap_or(contig.length as u64) as usize;
                bam_length / 2 // Circular genome handling
            } else {
                contig.length
            };

            // Fetch reads for this contig
            let reads = process_reads_for_contig(&mut bam, &contig.name, ref_length, &modules, seq_type)?;

            if reads.is_empty() {
                continue;
            }

            // Calculate coverage percentage
            let mut covered = vec![false; ref_length];
            for read in &reads {
                let start = (read.ref_start.max(0) as usize).min(ref_length);
                let end = (read.ref_end as usize).min(ref_length);
                for i in start..end {
                    covered[i] = true;
                }
            }
            let covered_bp: usize = covered.iter().filter(|&&x| x).count();
            let coverage_pct = (covered_bp as f64 / ref_length as f64) * 100.0;

            if coverage_pct < min_coverage {
                continue;
            }

            all_presences.push((contig.name.clone(), coverage_pct as f32));

            // Calculate features based on modules
            if modules.contains(&"coverage".to_string()) || modules.contains(&"assemblycheck".to_string()) {
                let coverage = calculate_coverage(&reads, ref_length);
                let values: Vec<f64> = coverage.iter().map(|&x| x as f64).collect();
                let (xs, ys) = compress_signal(&values, PlotType::Curve, step, z_thresh, deriv_thresh, max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push((contig.name.clone(), "coverage".to_string(), x, y));
                }
            }

            if modules.contains(&"phagetermini".to_string()) {
                let pt_features = calculate_phagetermini(&reads, ref_length, seq_type);

                for feature_name in ["coverage_reduced", "reads_starts", "reads_ends"] {
                    if let Some(data) = pt_features.get(feature_name) {
                        let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                        let pt = if feature_name == "coverage_reduced" {
                            PlotType::Curve
                        } else {
                            PlotType::Bars
                        };
                        let (xs, ys) = compress_signal(&values, pt, step, z_thresh, deriv_thresh, max_points);
                        for (x, y) in xs.into_iter().zip(ys) {
                            all_features.push((contig.name.clone(), feature_name.to_string(), x, y));
                        }
                    }
                }

                // Calculate tau
                if let (Some(cov_red), Some(starts), Some(ends)) = (
                    pt_features.get("coverage_reduced"),
                    pt_features.get("reads_starts"),
                    pt_features.get("reads_ends"),
                ) {
                    let tau: Vec<f64> = cov_red
                        .iter()
                        .zip(starts)
                        .zip(ends)
                        .map(|((&c, &s), &e)| {
                            if c > 0 {
                                (s + e) as f64 / c as f64
                            } else {
                                0.0
                            }
                        })
                        .collect();
                    let (xs, ys) = compress_signal(&tau, PlotType::Bars, step, z_thresh, deriv_thresh, max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push((contig.name.clone(), "tau".to_string(), x, y));
                    }
                }
            }

            if modules.contains(&"assemblycheck".to_string()) {
                let ac_features = calculate_assemblycheck(&reads, ref_length, seq_type);

                for feature_name in [
                    "left_clippings",
                    "right_clippings",
                    "insertions",
                    "deletions",
                    "mismatches",
                    "bad_orientations",
                ] {
                    if let Some(data) = ac_features.get(feature_name) {
                        let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                        let (xs, ys) = compress_signal(&values, PlotType::Bars, step, z_thresh, deriv_thresh, max_points);
                        for (x, y) in xs.into_iter().zip(ys) {
                            all_features.push((contig.name.clone(), feature_name.to_string(), x, y));
                        }
                    }
                }

                // Read lengths (long reads)
                if let (Some(sum_rl), Some(count_rl)) = (
                    ac_features.get("sum_read_lengths"),
                    ac_features.get("count_read_lengths"),
                ) {
                    let values: Vec<f64> = sum_rl
                        .iter()
                        .zip(count_rl)
                        .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                        .collect();
                    let (xs, ys) = compress_signal(&values, PlotType::Curve, step, z_thresh, deriv_thresh, max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push((contig.name.clone(), "read_lengths".to_string(), x, y));
                    }
                }

                // Insert sizes (short-paired)
                if let (Some(sum_is), Some(count_is)) = (
                    ac_features.get("sum_insert_sizes"),
                    ac_features.get("count_insert_sizes"),
                ) {
                    let values: Vec<f64> = sum_is
                        .iter()
                        .zip(count_is)
                        .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                        .collect();
                    let (xs, ys) = compress_signal(&values, PlotType::Curve, step, z_thresh, deriv_thresh, max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push((contig.name.clone(), "insert_sizes".to_string(), x, y));
                    }
                }
            }
        }

        // Build result dict
        let result = PyDict::new(py);

        // Features as list of tuples
        let features_list = PyList::new(
            py,
            all_features
                .iter()
                .map(|(c, f, p, v)| (c.as_str(), f.as_str(), *p, *v)),
        )?;
        result.set_item("features", features_list)?;

        // Presences as list of tuples
        let presences_list = PyList::new(
            py,
            all_presences.iter().map(|(c, p)| (c.as_str(), *p)),
        )?;
        result.set_item("presences", presences_list)?;

        result.set_item("sample_name", &sample_name)?;

        // Contigs info
        let contigs_list = PyList::new(
            py,
            contigs.iter().map(|c| {
                let d = PyDict::new(py);
                d.set_item("name", &c.name).ok();
                d.set_item("length", c.length).ok();
                d
            }),
        )?;
        result.set_item("contigs", contigs_list)?;

        Ok(result)
    }

    /// Parse a GenBank file and return contig information.
    ///
    /// Args:
    ///     genbank_path: Path to the GenBank file
    ///     annotation_tool: Annotation tool name (e.g., "pharokka")
    ///
    /// Returns:
    ///     list of dicts with keys: "name", "length", "annotation_tool"
    #[pyfunction]
    #[pyo3(signature = (genbank_path, annotation_tool=""))]
    fn parse_genbank_py<'py>(
        py: Python<'py>,
        genbank_path: &str,
        annotation_tool: &str,
    ) -> PyResult<Bound<'py, PyList>> {
        let (contigs, _) = parse_genbank(Path::new(genbank_path), annotation_tool)?;

        let result = PyList::new(
            py,
            contigs.iter().map(|c| {
                let d = PyDict::new(py);
                d.set_item("name", &c.name).ok();
                d.set_item("length", c.length).ok();
                d.set_item("annotation_tool", &c.annotation_tool).ok();
                d
            }),
        )?;

        Ok(result)
    }

    /// Process all BAM files in parallel with rayon and progress bar.
    ///
    /// This is the main high-level API that matches the native Rust CLI performance.
    /// It parses GenBank once, then processes all BAM files in parallel using rayon.
    ///
    /// Args:
    ///     genbank_path: Path to the GenBank annotation file
    ///     bam_dir: Path to directory containing BAM files (or single BAM file)
    ///     output_dir: Output directory for results
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
    #[pyo3(signature = (genbank_path, bam_dir, output_dir, modules, threads, annotation_tool="", min_coverage=50.0, step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000))]
    fn process_all_samples<'py>(
        py: Python<'py>,
        genbank_path: &str,
        bam_dir: &str,
        output_dir: &str,
        modules: Vec<String>,
        threads: usize,
        annotation_tool: &str,
        min_coverage: f64,
        step: usize,
        z_thresh: f64,
        deriv_thresh: f64,
        max_points: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        use crate::processing::{run_all_samples, ProcessConfig};

        let config = ProcessConfig {
            threads,
            min_coverage,
            step,
            z_thresh,
            deriv_thresh,
            max_points,
        };

        // Release the GIL and call the shared processing function
        let process_result = py.allow_threads(|| {
            run_all_samples(
                Path::new(genbank_path),
                Path::new(bam_dir),
                Path::new(output_dir),
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
    ///
    /// This module provides fast BAM processing and feature calculation
    /// implemented in Rust for use from Python.
    #[pymodule]
    fn mgfeatureviewer_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(detect_sequencing_type_py, m)?)?;
        m.add_function(wrap_pyfunction!(compress_signal_py, m)?)?;
        m.add_function(wrap_pyfunction!(calculate_coverage_py, m)?)?;
        m.add_function(wrap_pyfunction!(process_bam_file, m)?)?;
        m.add_function(wrap_pyfunction!(parse_genbank_py, m)?)?;
        m.add_function(wrap_pyfunction!(process_all_samples, m)?)?;

        // Add module docstring
        m.add("__doc__", "MGFeatureViewer Rust bindings - fast BAM processing and feature calculation")?;

        Ok(())
    }
}
