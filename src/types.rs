//! Type definitions for MGFeatureViewer Rust Calculator.
//!
//! This module contains all the struct and enum definitions used throughout the application,
//! as well as the single source of truth for feature/variable configuration.
//!
//! # For Python developers new to Rust:
//!
//! ## Key Rust concepts used here:
//!
//! - **`enum`**: Like Python's Enum class, but more powerful. Each variant can hold data.
//! - **`struct`**: Similar to a Python dataclass or namedtuple - a container for related fields.
//! - **`impl`**: Adds methods to a struct or enum (like defining methods inside a Python class).
//! - **`pub`**: Makes something public (accessible from other modules). Without it, it's private.
//! - **`&str`** vs **`String`**: `&str` is a borrowed reference to text (like a view), `String` owns the text.
//! - **`&'static str`**: A string that lives for the entire program (like string literals).
//! - **`Option<T>`**: Rust's way of handling "maybe there's a value" - like Python's `Optional[T]` or `None`.

use std::collections::HashMap;

// ============================================================================
// Sequencing Types
// ============================================================================

/// Represents the type of sequencing technology detected from the BAM file.
///
/// This is determined by examining the first few reads:
/// - Long reads (>1000bp) → `Long` (e.g., PacBio, Nanopore)
/// - Paired-end short reads → `ShortPaired` (e.g., Illumina paired-end)
/// - Single-end short reads → `ShortSingle` (e.g., Illumina single-end)
///
/// # Python equivalent:
/// ```python
/// # Return values from find_sequencing_type_from_bam() in calculating_data.py:55-78
/// sequencing_type = "long"  # or "short_paired" or "short_single"
/// ```
///
/// # Rust concept - `#[derive(...)]`:
/// This automatically generates common functionality:
/// - `Clone`: Can create copies (like Python's copy.copy)
/// - `Copy`: Can be copied implicitly (only for small, simple types)
/// - `PartialEq, Eq`: Can compare with == and !=
/// - `Debug`: Can be printed with {:?} for debugging
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SequencingType {
    /// Long-read sequencing (PacBio, Nanopore) - reads typically >1000bp
    Long,
    /// Paired-end short-read sequencing (Illumina) - two reads per fragment
    ShortPaired,
    /// Single-end short-read sequencing - one read per fragment
    ShortSingle,
}

/// # Rust concept - `impl` block:
/// This is where we add methods to our enum, similar to defining methods in a Python class.
/// In Python you'd write: `class SequencingType: def is_long(self): ...`
impl SequencingType {
    /// Check if this is a long-read sequencing type.
    ///
    /// Used to determine which features to calculate (e.g., read_lengths only for long reads).
    ///
    /// # Rust concept - `#[inline]`:
    /// Suggests to the compiler to insert this function's code directly where it's called,
    /// avoiding the overhead of a function call. Good for tiny functions called many times.
    ///
    /// # Rust concept - `matches!` macro:
    /// A concise way to check if a value matches a pattern. Equivalent to:
    /// ```rust
    /// match self {
    ///     Self::Long => true,
    ///     _ => false,
    /// }
    /// ```
    #[inline]
    pub fn is_long(&self) -> bool {
        matches!(self, Self::Long)
    }

    /// Check if this is paired-end short-read sequencing.
    ///
    /// Used to determine which features to calculate (e.g., insert_sizes only for paired reads).
    #[inline]
    pub fn is_short_paired(&self) -> bool {
        matches!(self, Self::ShortPaired)
    }
}

// ============================================================================
// Plot Types
// ============================================================================

/// How a feature should be visualized in the viewer.
///
/// - `Curve`: Continuous line plot (good for coverage, smoothly varying data)
/// - `Bars`: Bar/stem plot (good for discrete events like read starts, clippings)
///
/// # Python equivalent:
/// ```python
/// # Values in FEATURE_TYPES dict in calculating_data.py:17-31
/// plot_type = "curve"  # or "bars"
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PlotType {
    /// Continuous line plot - for smoothly varying data like coverage
    Curve,
    /// Discrete bar/stem plot - for event counts like read starts, clippings
    Bars,
}

impl PlotType {
    /// Convert to string representation for database storage.
    ///
    /// # Rust concept - `&'static str`:
    /// Returns a reference to a string that exists for the entire program lifetime.
    /// These strings ("curve", "bars") are compiled into the binary itself.
    #[inline]
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Curve => "curve",
            Self::Bars => "bars",
        }
    }
}

// ============================================================================
// Variable Configuration - Single Source of Truth
// ============================================================================

/// Configuration for each variable/feature we calculate.
///
/// This struct holds all the metadata needed for:
/// - Database storage (name, module)
/// - Visualization (subplot, plot_type, color, line_alpha, etc.)
///
/// # Rust concept - `&'static str` fields:
/// Using `&'static str` instead of `String` means these are references to
/// string literals compiled into the program. This is more efficient because:
/// 1. No heap allocation needed
/// 2. No need to free memory
/// 3. Can be used in `const` contexts
#[derive(Clone, Copy, Debug)]
pub struct VariableConfig {
    /// Internal name used in database tables (e.g., "coverage", "reads_starts")
    pub name: &'static str,
    /// Subplot grouping for visualization (e.g., "Coverage", "Reads termini")
    pub subplot: &'static str,
    /// Module this variable belongs to (e.g., "Coverage", "Phage termini", "Assembly check")
    pub module: &'static str,
    /// How to visualize this variable
    pub plot_type: PlotType,
    /// Hex color code for plotting (e.g., "#333333")
    pub color: &'static str,
    /// Line/marker opacity (0.0 = transparent, 1.0 = opaque)
    pub alpha: f64,
    /// Fill opacity for area under curves
    pub fill_alpha: f64,
    /// Line width or marker size
    pub size: f64,
    /// Human-readable title for legends
    pub title: &'static str,
    /// Optional help text for tooltips
    pub help: Option<&'static str>,
}

/// All variable configurations - single source of truth for the entire application.
///
/// # Why a single source of truth?
/// By defining all variables here, we ensure consistency between:
/// - Feature calculation code
/// - Database schema
/// - Visualization settings
///
/// # Rust concept - `const`:
/// `const` values are computed at compile time and embedded in the binary.
/// They cannot be changed at runtime. This is different from `let` which
/// creates a runtime variable.
///
/// # Rust concept - `&[VariableConfig]`:
/// This is a "slice" - a view into a contiguous sequence. The `&` means
/// we're borrowing it (not owning it). Think of it like a Python list view.
pub const VARIABLES: &[VariableConfig] = &[
    // Coverage module - basic read depth
    VariableConfig { name: "primary_reads", subplot: "Primary alignments", module: "Coverage", plot_type: PlotType::Curve, color: "#333333", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Primary reads", help: None },
    VariableConfig { name: "secondary_reads", subplot: "Other alignments", module: "Coverage", plot_type: PlotType::Curve, color: "#1f77b4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Secondary reads", help: None },
    VariableConfig { name: "supplementary_reads", subplot: "Other alignments", module: "Coverage", plot_type: PlotType::Curve, color: "#B45C1F", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Supplementary reads", help: None },

    // Phage termini module - for detecting phage DNA packaging sites
    VariableConfig { name: "coverage_reduced", subplot: "Coverage reduced", module: "Phage termini", plot_type: PlotType::Curve, color: "#00c53b", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Coverage reduced", help: Some("Retains primary mappings starting with an exact match (and ending with an exact match for long reads)") },
    VariableConfig { name: "reads_starts", subplot: "Reads termini", module: "Phage termini", plot_type: PlotType::Bars, color: "#215732", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read Starts", help: None },
    VariableConfig { name: "reads_ends", subplot: "Reads termini", module: "Phage termini", plot_type: PlotType::Bars, color: "#572146", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read Ends", help: None },
    VariableConfig { name: "tau", subplot: "Tau", module: "Phage termini", plot_type: PlotType::Bars, color: "#44883e", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Tau", help: None },

    // Assembly check module - for detecting assembly errors
    VariableConfig { name: "read_lengths", subplot: "Read lengths", module: "Assembly check", plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Read Lengths", help: None },
    VariableConfig { name: "insert_sizes", subplot: "Insert sizes", module: "Assembly check", plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Insert Sizes", help: None },
    VariableConfig { name: "bad_orientations", subplot: "Bad orientations", module: "Assembly check", plot_type: PlotType::Bars, color: "#c94009", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Bad Orientations", help: None },
    VariableConfig { name: "left_clippings", subplot: "Clippings", module: "Assembly check", plot_type: PlotType::Bars, color: "#8e43e7", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Left Clippings", help: None },
    VariableConfig { name: "right_clippings", subplot: "Clippings", module: "Assembly check", plot_type: PlotType::Bars, color: "#9CE743", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Right Clippings", help: None },
    VariableConfig { name: "insertions", subplot: "Indels", module: "Assembly check", plot_type: PlotType::Bars, color: "#e50001", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Insertions", help: Some("Extra bases are present in the read but not in the reference") },
    VariableConfig { name: "deletions", subplot: "Indels", module: "Assembly check", plot_type: PlotType::Bars, color: "#00E5E4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Deletions", help: Some("A stretch of the reference has no corresponding bases in the read") },
    VariableConfig { name: "mismatches", subplot: "Mismatches", module: "Assembly check", plot_type: PlotType::Bars, color: "#5a0f0b", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mismatches", help: None },
];

// ============================================================================
// Feature Name Constants - Quick Lookup Lists
// ============================================================================

/// Feature names calculated by the phagetermini module.
///
/// These are used for phage terminus detection (DNA packaging sites).
/// - `coverage_reduced`: Coverage counting only "clean" reads (no clipping/mismatch at ends)
/// - `reads_starts`: Count of read 5' ends at each position
/// - `reads_ends`: Count of read 3' ends at each position
pub const PHAGETERMINI_FEATURES: &[&str] = &["coverage_reduced", "reads_starts", "reads_ends"];

/// Feature names calculated by the assemblycheck module.
///
/// These help identify potential assembly errors:
/// - `left_clippings`: Soft/hard clips at read starts (potential breakpoints)
/// - `right_clippings`: Soft/hard clips at read ends
/// - `insertions`: Insertion events from CIGAR
/// - `deletions`: Deletion events from CIGAR
/// - `mismatches`: Base mismatches from MD tag
/// - `bad_orientations`: Reads not in proper pair orientation
pub const ASSEMBLYCHECK_FEATURES: &[&str] = &[
    "left_clippings",
    "right_clippings",
    "insertions",
    "deletions",
    "mismatches",
    "bad_orientations",
];

// ============================================================================
// Helper Functions for Variables
// ============================================================================

/// Get the plot type for a feature by its name.
///
/// # Arguments
/// * `feature` - The feature name (e.g., "coverage", "reads_starts")
///
/// # Returns
/// The PlotType for visualization, defaulting to Bars if not found.
///
/// # Rust concept - Iterator chain:
/// ```rust
/// VARIABLES.iter()           // Create an iterator over the slice
///     .find(|v| ...)         // Find first matching element (returns Option)
///     .map(|v| v.plot_type)  // Transform the found value (if any)
///     .unwrap_or(...)        // Provide default if None
/// ```
/// This is similar to Python's: `next((v.plot_type for v in VARIABLES if v.name == feature), PlotType.Bars)`
#[inline]
pub fn get_plot_type(feature: &str) -> PlotType {
    VARIABLES
        .iter()
        .find(|v| v.name == feature)
        .map(|v| v.plot_type)
        .unwrap_or(PlotType::Bars)
}

/// Get the full variable configuration by name.
///
/// # Returns
/// `Some(&VariableConfig)` if found, `None` if not found.
///
/// # Rust concept - `Option<&T>`:
/// Instead of returning `None` or raising an exception like Python,
/// Rust uses `Option` to explicitly handle "might not exist" cases.
/// The caller must handle both `Some(value)` and `None`.
#[inline]
pub fn get_variable(name: &str) -> Option<&'static VariableConfig> {
    VARIABLES.iter().find(|v| v.name == name)
}

/// Get an iterator over all variables belonging to a specific module.
///
/// # Example
/// ```rust
/// for var in variables_for_module("Phage termini") {
///     println!("{}: {}", var.name, var.title);
/// }
/// ```
///
/// # Rust concept - `impl Iterator`:
/// Instead of returning a concrete type, we return "something that implements Iterator".
/// This is like Python's type hint `Iterator[VariableConfig]`.
///
/// # Rust concept - `+ '_`:
/// The `'_` is a lifetime annotation. It tells Rust that the returned iterator
/// borrows from `module` and cannot outlive it. Don't worry too much about this -
/// the compiler usually figures it out, but sometimes needs hints.
pub fn variables_for_module(module: &str) -> impl Iterator<Item = &'static VariableConfig> + '_ {
    VARIABLES.iter().filter(move |v| v.module == module)
}

/// Generate the database table name for a feature.
///
/// Convention: `Feature_<name>` (e.g., "Feature_coverage", "Feature_reads_starts")
#[inline]
pub fn feature_table_name(feature: &str) -> String {
    format!("Feature_{}", feature)
}

// ============================================================================
// Data Structures for Processing Pipeline
// ============================================================================

/// Information about a contig (chromosome/sequence) from the GenBank file.
///
/// # Bioinformatics context:
/// A "contig" is a contiguous sequence - could be a chromosome, plasmid,
/// or assembled sequence fragment. In phage genomics, each phage genome
/// is typically one contig.
#[derive(Clone, Debug)]
pub struct ContigInfo {
    /// Unique identifier for the contig (e.g., "phage_lambda", "contig_1")
    pub name: String,
    /// Length in base pairs
    pub length: usize,
    /// Tool used for annotation (e.g., "pharokka", "prokka")
    pub annotation_tool: String,
}

/// A gene/feature annotation from the GenBank file.
///
/// # Bioinformatics context:
/// Represents a CDS (coding sequence), gene, or other annotated feature
/// on a contig. Used for displaying gene tracks in the viewer.
#[derive(Clone, Debug)]
pub struct FeatureAnnotation {
    /// Database ID of the parent contig
    pub contig_id: i64,
    /// Start position (1-based, as in GenBank)
    pub start: i64,
    /// End position (1-based, inclusive)
    pub end: i64,
    /// Strand: 1 for forward (+), -1 for reverse (-)
    pub strand: i64,
    /// Feature type (e.g., "CDS", "tRNA", "gene")
    pub feature_type: String,
    /// Product description (e.g., "terminase large subunit")
    /// Option<String> means this field might be empty/missing
    pub product: Option<String>,
    /// Functional annotation
    pub function: Option<String>,
    /// PHROG database annotation (phage-specific)
    pub phrog: Option<String>,
}

/// A single data point for a calculated feature.
///
/// After calculating features (coverage, reads_starts, etc.) and compressing
/// the signal, we store individual (position, value) points.
///
/// # Why compress?
/// A 50kb phage genome at single-base resolution = 50,000 points per feature.
/// With 10 features × 100 samples = 50 million points! Compression reduces
/// A single feature data point (run-length encoded).
///
/// Each FeaturePoint represents a run of consecutive positions with the same value.
/// For singleton points (spikes), start_pos == end_pos.
///
/// Before compression, genomic signals have millions of (position, value) pairs.
/// Run-length encoding with adaptive thresholding reduces storage by encoding
/// constant regions as single runs while preserving rapid changes as short runs.
#[derive(Clone, Debug)]
pub struct FeaturePoint {
    /// Which contig this point belongs to
    pub contig_name: String,
    /// Which feature (e.g., "coverage", "reads_starts")
    pub feature: String,
    /// First position in the run (1-indexed, inclusive)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// The value for this run
    pub value: f32,
}

/// Records whether a contig was detected in a sample.
///
/// A contig is "present" in a sample if it has sufficient coverage.
/// This is used for presence/absence matrices across samples.
#[derive(Clone, Debug)]
pub struct PresenceData {
    /// The contig name
    pub contig_name: String,
    /// Percentage of bases with at least 1x coverage (0-100)
    pub coverage_pct: f32,
}

/// Type alias for feature calculation results.
///
/// Maps feature names to their per-position values.
///
/// # Rust concept - `type` alias:
/// This creates a shorthand name for a complex type. Instead of writing
/// `HashMap<String, Vec<u64>>` everywhere, we can write `FeatureMap`.
///
/// # Example:
/// ```rust
/// let mut results: FeatureMap = HashMap::new();
/// results.insert("coverage".to_string(), vec![10, 15, 20, 18, ...]);
/// results.insert("reads_starts".to_string(), vec![0, 1, 0, 0, 5, ...]);
/// ```
pub type FeatureMap = HashMap<String, Vec<u64>>;

// ============================================================================
// Statistics Helpers
// ============================================================================

/// Calculate mean and standard deviation of a slice of values.
///
/// # Arguments
/// * `values` - A slice of f64 values (like a numpy array view)
///
/// # Returns
/// A tuple `(mean, std)` where std is at least 1e-9 to avoid division by zero.
///
/// # Rust concept - Tuple return:
/// Rust functions can return multiple values using tuples, similar to Python.
/// The caller destructures with: `let (mean, std) = mean_std(&values);`
///
/// # Rust concept - Iterator methods:
/// Instead of explicit loops, Rust idiomatically uses iterator adapters:
/// - `.iter()` - create an iterator
/// - `.sum()` - add up all values
/// - `.map(|x| ...)` - transform each value (like Python's map() or list comprehension)
#[inline]
pub fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len();

    // Handle empty input
    if n == 0 {
        return (0.0, 1e-9);
    }

    // Calculate mean: sum / count
    // The `::<f64>` is a "turbofish" - it tells .sum() what type to produce
    let mean = values.iter().sum::<f64>() / n as f64;

    // Calculate variance: mean of squared differences from mean
    // .map() transforms each value, .sum() adds them up
    let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;

    // Standard deviation is sqrt of variance
    // .max(1e-9) ensures we never return 0 (avoids division by zero later)
    let std = variance.sqrt().max(1e-9);

    (mean, std)
}
