//! Type definitions for MGFeatureViewer Rust Calculator.
//!
//! This module contains all the struct and enum definitions used throughout the application.

use std::collections::HashMap;

/// Sequencing type detected from BAM file.
/// Python equivalent: return values from `find_sequencing_type_from_bam()` in calculating_data.py:55-78
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SequencingType {
    Long,
    ShortPaired,
    ShortSingle,
}

/// Plot type for feature visualization.
/// Python equivalent: values in `FEATURE_TYPES` dict in calculating_data.py:17-31
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum PlotType {
    Curve,
    Bars,
}

/// Configuration for a feature type.
pub struct FeatureConfig {
    pub plot_type: PlotType,
}

/// Get feature configuration (plot type).
/// Python equivalent: `FEATURE_TYPES` dict in calculating_data.py:17-31
pub fn get_feature_config(feature: &str) -> FeatureConfig {
    // calculating_data.py:17-31 - FEATURE_TYPES dict
    // Maps feature names to "curve" or "bars" plot type
    let plot_type = match feature {
        "coverage" | "coverage_reduced" | "read_lengths" | "insert_sizes" => PlotType::Curve,  // py:18-19, 24-25
        _ => PlotType::Bars,  // py:20-23, 26-31 (reads_starts, reads_ends, tau, bad_orientations, clippings, indels, mismatches)
    };
    FeatureConfig { plot_type }
}

/// Information about a contig from the GenBank file.
#[derive(Clone, Debug)]
pub struct ContigInfo {
    pub name: String,
    pub length: usize,
    pub annotation_tool: String,
}

/// Feature annotation from the GenBank file.
#[derive(Clone, Debug)]
pub struct FeatureAnnotation {
    pub contig_id: i64,
    pub start: i64,
    pub end: i64,
    pub strand: i64,
    pub feature_type: String,
    pub product: Option<String>,
    pub function: Option<String>,
    pub phrog: Option<String>,
}

/// Read preprocessed data for a single read.
/// Python equivalent: data extracted in `preprocess_reads()` in calculating_data.py:497-602
#[derive(Clone, Debug)]
pub struct ReadData {
    pub ref_start: i64,
    pub ref_end: i64,
    pub query_length: i32,
    pub template_length: i32,
    pub is_read1: bool,
    pub is_proper_pair: bool,
    pub is_reverse: bool,
    pub cigar: Vec<(u32, u32)>, // (op, len)
    pub md_tag: Option<Vec<u8>>,
}

/// Feature data point for output.
#[derive(Clone, Debug)]
pub struct FeaturePoint {
    pub contig_name: String,
    pub feature: String,
    pub position: i32,
    pub value: f32,
}

/// Presence data for a contig in a sample.
#[derive(Clone, Debug)]
pub struct PresenceData {
    pub contig_name: String,
    pub coverage_pct: f32,
}

/// Result of feature calculations, keyed by feature name.
pub type FeatureMap = HashMap<String, Vec<u64>>;
