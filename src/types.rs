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

/// Variable configuration for database and plotting.
pub struct VariableConfig {
    pub name: &'static str,
    pub subplot: &'static str,
    pub module: &'static str,
    pub plot_type: PlotType,
    pub color: &'static str,
    pub alpha: f64,
    pub fill_alpha: f64,
    pub size: f64,
    pub title: &'static str,
}

/// All variable configurations - single source of truth.
pub const VARIABLES: &[VariableConfig] = &[
    VariableConfig { name: "coverage", subplot: "Coverage", module: "Coverage", plot_type: PlotType::Curve, color: "#333333", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Coverage depth" },
    VariableConfig { name: "coverage_reduced", subplot: "Coverage reduced", module: "Phage termini", plot_type: PlotType::Curve, color: "#00c53b", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Coverage reduced" },
    VariableConfig { name: "reads_starts", subplot: "Reads termini", module: "Phage termini", plot_type: PlotType::Bars, color: "#215732", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read Starts" },
    VariableConfig { name: "reads_ends", subplot: "Reads termini", module: "Phage termini", plot_type: PlotType::Bars, color: "#6cc24a", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read Ends" },
    VariableConfig { name: "tau", subplot: "Tau", module: "Phage termini", plot_type: PlotType::Bars, color: "#44883e", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Tau" },
    VariableConfig { name: "read_lengths", subplot: "Read lengths", module: "Assembly check", plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Read Lengths" },
    VariableConfig { name: "insert_sizes", subplot: "Insert sizes", module: "Assembly check", plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Insert Sizes" },
    VariableConfig { name: "bad_orientations", subplot: "Bad orientations", module: "Assembly check", plot_type: PlotType::Bars, color: "#c94009", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Bad Orientations" },
    VariableConfig { name: "left_clippings", subplot: "Clippings", module: "Assembly check", plot_type: PlotType::Bars, color: "#7f0091", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Left Clippings" },
    VariableConfig { name: "right_clippings", subplot: "Clippings", module: "Assembly check", plot_type: PlotType::Bars, color: "#8e43e7", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Right Clippings" },
    VariableConfig { name: "insertions", subplot: "Indels", module: "Assembly check", plot_type: PlotType::Bars, color: "#e50001", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Insertions" },
    VariableConfig { name: "deletions", subplot: "Indels", module: "Assembly check", plot_type: PlotType::Bars, color: "#97011a", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Deletions" },
    VariableConfig { name: "mismatches", subplot: "Mismatches", module: "Assembly check", plot_type: PlotType::Bars, color: "#5a0f0b", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mismatches" },
];

/// Get plot type for a feature.
pub fn get_plot_type(feature: &str) -> PlotType {
    VARIABLES.iter()
        .find(|v| v.name == feature)
        .map(|v| v.plot_type)
        .unwrap_or(PlotType::Bars)
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
