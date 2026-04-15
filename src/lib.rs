//! theBIGbam Rust Library
//!
//! This module provides both a Rust library and Python bindings (via PyO3).
//!
//! ## Python Usage
//!
//! ```python
//! import thebigbam_rs as mgfv
//!
//! # Process all BAM files in parallel (main API)
//! result = mgfv.process_all_samples(
//!     genbank_path="annotation.gbk",
//!     bam_dir="/path/to/bams",
//!     output_dir="/path/to/output",
//!     modules=["Coverage", "Misalignment", "Long-reads", "Paired-reads", "Phage termini"],
//!     threads=8,
//! )
//! # result = {"samples_processed": 10, "samples_failed": 0, "total_time": 123.4}
//! ```

pub mod bam_reader;
pub mod blob;
pub mod cigar;
pub mod circular;
pub mod compress;
pub mod db;
pub mod features;
pub mod gc_content;
pub mod mag_blob;
pub mod parser;
pub mod processing;
pub mod processing_completeness;
pub mod processing_phage_packaging;
pub mod types;

// Re-export main types for library users
pub use bam_reader::{detect_sequencing_type, process_contig_streaming};
pub use cigar::{CigarOp, MdTag};
pub use circular::{increment_circular, increment_range};
pub use compress::compress_signal_with_reference;
pub use features::{process_read, FeatureArrays, ModuleFlags};
pub use parser::{parse_annotations, parse_fasta, parse_genbank, parse_gff3};
pub use processing::{run_all_samples, ProcessConfig};
pub use processing_phage_packaging::PhageTerminiConfig;
pub use types::{
    mean_std, ContigInfo, PlotType, PresenceData, SequencingType,
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
    ///     modules: List of modules to compute: "Coverage", "Misalignment", "RNA", "Long-reads", "Paired-reads", "Phage termini"
    ///     threads: Number of threads to use
    ///     min_aligned_fraction: Minimum alignment-length coverage fraction for contig inclusion (default 50.0)
    ///     min_coverage_depth: Minimum mean coverage depth for contig inclusion (default 0.0 = disabled)
    ///     compress_ratio: Compression ratio threshold (default 10.0)
    ///     create_indexes: Whether to create database indexes (default True)
    ///     assembly_path: Path to assembly FASTA file for autoblast (default "")
    ///
    /// Returns:
    ///     dict with keys:
    ///         - "samples_processed": int
    ///         - "samples_failed": int
    ///         - "total_time": float (seconds)
    #[pyfunction]
    #[pyo3(signature = (genbank_path, bam_files, output_db, modules, threads, sequencing_type=None, min_aligned_fraction=50.0, min_coverage_depth=0.0, curve_ratio=10.0, bar_ratio=10.0, create_indexes=true, assembly_path="", extend_db="", min_occurrences=2, enable_timing=false, view="contig", mag_manifest=vec![]))]
    fn process_all_samples<'py>(
        py: Python<'py>,
        genbank_path: &str,
        bam_files: Vec<String>,
        output_db: &str,
        modules: Vec<String>,
        threads: usize,
        sequencing_type: Option<&str>,
        min_aligned_fraction: f64,
        min_coverage_depth: f64,
        curve_ratio: f64,
        bar_ratio: f64,
        create_indexes: bool,
        assembly_path: &str,
        extend_db: &str,
        min_occurrences: u32,
        enable_timing: bool,
        view: &str,
        mag_manifest: Vec<(String, String, String)>,
    ) -> PyResult<Bound<'py, PyDict>> {
        use crate::gc_content::GCParams;
        use crate::processing::{run_all_samples, MagInput, ProcessConfig, ViewMode};
        use crate::processing_phage_packaging::PhageTerminiConfig;
        use std::path::PathBuf;

        // Parse sequencing type: Some(type) if user provided valid value, None for per-sample auto-detection
        let seq_type = sequencing_type.and_then(|s| ProcessConfig::parse_sequencing_type(s));

        let view_mode = ViewMode::from_str(view);
        let mag_manifest: Vec<MagInput> = mag_manifest
            .into_iter()
            .map(|(name, gb, asm)| MagInput {
                name,
                gb_path: if gb.is_empty() { None } else { Some(PathBuf::from(gb)) },
                asm_path: if asm.is_empty() { None } else { Some(PathBuf::from(asm)) },
            })
            .collect();

        if view_mode == ViewMode::Mag && mag_manifest.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "view='mag' requires a non-empty mag_manifest",
            ));
        }

        let config = ProcessConfig {
            threads,
            min_aligned_fraction,
            min_coverage_depth,
            curve_ratio,
            bar_ratio,
            sequencing_type: seq_type,
            phagetermini_config: PhageTerminiConfig::default(),
            gc_params: GCParams::default(),
            min_occurrences,
            enable_timing,
            view_mode,
            mag_manifest,
            contig_drops_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
        };

        // Convert string paths to PathBuf
        let bam_paths: Vec<PathBuf> = bam_files.iter().map(|s| PathBuf::from(s)).collect();

        // Release the GIL and call the shared processing function
        let process_result = py.allow_threads(|| {
            run_all_samples(
                Path::new(genbank_path),
                Path::new(assembly_path),
                &bam_paths,
                Path::new(output_db),
                &modules,
                &config,
                create_indexes,
                Path::new(extend_db),
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

    // ========================================================================
    // Fast chunk decoders for Python (replaces pure-Python _varint_decode)
    // ========================================================================

    /// Varint (LEB128) decode a byte slice into unsigned 64-bit values.
    fn varint_decode_slice(data: &[u8]) -> Vec<u64> {
        let mut values = Vec::with_capacity(data.len() / 2);
        let mut i = 0;
        while i < data.len() {
            let mut val: u64 = 0;
            let mut shift = 0u32;
            loop {
                if i >= data.len() { break; }
                let byte = data[i];
                i += 1;
                val |= ((byte & 0x7F) as u64) << shift;
                shift += 7;
                if byte & 0x80 == 0 { break; }
            }
            values.push(val);
        }
        values
    }

    /// Decode a dense chunk: zstd decompress → varint → zigzag → delta accumulate.
    ///
    /// Input: raw compressed bytes from Feature_blob_chunk.Data / Contig_blob_chunk.Data
    /// Returns: list[int] — cumulative delta-decoded signed values
    #[pyfunction]
    fn decode_dense_chunk(compressed: &[u8]) -> PyResult<Vec<i64>> {
        let decompressed = zstd::stream::decode_all(std::io::Cursor::new(compressed))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(
                format!("zstd decompression failed: {e}")))?;

        let unsigned_vals = varint_decode_slice(&decompressed);

        // Zigzag decode + delta accumulate in one pass
        let mut result = Vec::with_capacity(unsigned_vals.len());
        let mut acc: i64 = 0;
        for v in unsigned_vals {
            let signed = ((v >> 1) as i64) ^ -((v & 1) as i64);
            acc += signed;
            result.push(acc);
        }
        Ok(result)
    }

    /// Decode the positions and values from a sparse chunk.
    ///
    /// Sparse chunk format (uncompressed outer):
    ///   [event_count: u32 LE]
    ///   [pos_compressed_size: u32 LE][pos_compressed: zstd(varint(zigzag(delta)))]
    ///   [val_compressed_size: u32 LE][val_compressed: zstd(varint(zigzag))]
    ///   ... metadata follows (parsed in Python)
    ///
    /// Returns: (positions: list[int], values: list[int])
    /// Positions are delta-decoded, values are zigzag-only (no delta).
    /// Both truncated to event_count.
    #[pyfunction]
    fn decode_sparse_chunk(raw_bytes: &[u8]) -> PyResult<(Vec<i64>, Vec<i64>)> {
        if raw_bytes.len() < 4 {
            return Ok((vec![], vec![]));
        }

        let event_count = u32::from_le_bytes(
            [raw_bytes[0], raw_bytes[1], raw_bytes[2], raw_bytes[3]]) as usize;
        let mut off = 4usize;

        if event_count == 0 {
            return Ok((vec![], vec![]));
        }

        // Positions: size + zstd-compressed varint data
        if off + 4 > raw_bytes.len() { return Ok((vec![], vec![])); }
        let pos_size = u32::from_le_bytes(
            [raw_bytes[off], raw_bytes[off+1], raw_bytes[off+2], raw_bytes[off+3]]) as usize;
        off += 4;
        if off + pos_size > raw_bytes.len() { return Ok((vec![], vec![])); }

        let pos_decompressed = zstd::stream::decode_all(
            std::io::Cursor::new(&raw_bytes[off..off + pos_size]))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(
                format!("zstd pos decompression failed: {e}")))?;
        off += pos_size;

        let pos_unsigned = varint_decode_slice(&pos_decompressed);
        // Zigzag + delta decode positions
        let mut positions = Vec::with_capacity(event_count);
        let mut acc: i64 = 0;
        for v in &pos_unsigned {
            let signed = ((*v >> 1) as i64) ^ -((*v & 1) as i64);
            acc += signed;
            positions.push(acc);
        }
        positions.truncate(event_count);

        // Values: size + zstd-compressed varint data
        if off + 4 > raw_bytes.len() { return Ok((positions, vec![])); }
        let val_size = u32::from_le_bytes(
            [raw_bytes[off], raw_bytes[off+1], raw_bytes[off+2], raw_bytes[off+3]]) as usize;
        off += 4;
        if off + val_size > raw_bytes.len() { return Ok((positions, vec![])); }

        let val_decompressed = zstd::stream::decode_all(
            std::io::Cursor::new(&raw_bytes[off..off + val_size]))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(
                format!("zstd val decompression failed: {e}")))?;

        let val_unsigned = varint_decode_slice(&val_decompressed);
        // Zigzag decode only (no delta for sparse values)
        let mut values: Vec<i64> = val_unsigned.iter()
            .map(|&v| ((v >> 1) as i64) ^ -((v & 1) as i64))
            .collect();
        values.truncate(event_count);

        Ok((positions, values))
    }

    /// theBIGbam Rust bindings for Python.
    #[pymodule]
    fn thebigbam_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(process_all_samples, m)?)?;
        m.add_function(wrap_pyfunction!(decode_dense_chunk, m)?)?;
        m.add_function(wrap_pyfunction!(decode_sparse_chunk, m)?)?;
        m.add("__doc__", "theBIGbam Rust bindings - fast BAM processing")?;
        Ok(())
    }
}
