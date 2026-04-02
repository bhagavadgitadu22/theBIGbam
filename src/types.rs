//! Type definitions for theBIGbam Rust Calculator.
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

    /// Convert to database-friendly string representation.
    #[inline]
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Long => "long",
            Self::ShortPaired => "paired-short",
            Self::ShortSingle => "single-short",
        }
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
    /// Order within the module for display (1-based)
    pub module_order: i32,
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
    // Genome module - genomic properties
    // Repeat count subplot (button 1)
    VariableConfig { name: "direct_repeat_count", subplot: "Repeat count", module: "Genome", module_order: 1, plot_type: PlotType::Curve, color: "#c1121f", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Direct repeat count", help: Some("Direct repeats detected by self-BLAST (e.g., terminal repeats)") },
    VariableConfig { name: "inverted_repeat_count", subplot: "Repeat count", module: "Genome", module_order: 1, plot_type: PlotType::Curve, color: "#12C1B4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Inverted repeat count", help: Some("Inverted repeats detected by self-BLAST (e.g., terminal repeats)") },
    // Max repeat identity subplot (button 2)
    VariableConfig { name: "direct_repeat_identity", subplot: "Max repeat identity", module: "Genome", module_order: 2, plot_type: PlotType::Curve, color: "#c1121f", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Max direct repeat identity", help: Some("Direct repeats detected by self-BLAST (e.g., terminal repeats)") },
    VariableConfig { name: "inverted_repeat_identity", subplot: "Max repeat identity", module: "Genome", module_order: 2, plot_type: PlotType::Curve, color: "#12C1B4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Max inverted repeat identity", help: Some("Inverted repeats detected by self-BLAST (e.g., terminal repeats)") },
    VariableConfig { name: "gc_content", subplot: "GC content", module: "Genome", module_order: 3, plot_type: PlotType::Curve, color: "#693efe", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "GC content", help: None },
    VariableConfig { name: "gc_skew", subplot: "GC skew", module: "Genome", module_order: 4, plot_type: PlotType::Curve, color: "#C8A2C8", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "GC skew", help: None },   
    
    // Coverage module - basic read depth
    VariableConfig { name: "primary_reads", subplot: "Primary alignments", module: "Coverage", module_order: 1, plot_type: PlotType::Curve, color: "#333333", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Primary reads", help: None },
    VariableConfig { name: "primary_reads_plus_only", subplot: "Alignments by strand", module: "Coverage", module_order: 2, plot_type: PlotType::Curve, color: "#a1665e", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Primary reads (+ only)", help: None },
    VariableConfig { name: "primary_reads_minus_only", subplot: "Alignments by strand", module: "Coverage", module_order: 2, plot_type: PlotType::Curve, color: "#5E99A1", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Primary reads (- only)", help: None },
    VariableConfig { name: "secondary_reads", subplot: "Other alignments", module: "Coverage", module_order: 3, plot_type: PlotType::Curve, color: "#1f77b4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Secondary reads", help: None },
    VariableConfig { name: "supplementary_reads", subplot: "Other alignments", module: "Coverage", module_order: 3, plot_type: PlotType::Curve, color: "#B45C1F", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Supplementary reads", help: None },
    VariableConfig { name: "mapq", subplot: "MAPQ", module: "Coverage", module_order: 4, plot_type: PlotType::Curve, color: "#77dd77", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mapping quality", help: None },
        
    // Per position errors from reads (Assembly check)
    VariableConfig { name: "left_clippings", subplot: "Clippings", module: "Misalignment", module_order: 1, plot_type: PlotType::Bars, color: "#8e43e7", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Left clippings", help: Some("Extra bases are missing on the left of the contig") },
    VariableConfig { name: "right_clippings", subplot: "Clippings", module: "Misalignment", module_order: 1, plot_type: PlotType::Bars, color: "#9CE743", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Right clippings", help: Some("Extra bases are missing on the right of the contig") },
    VariableConfig { name: "insertions", subplot: "Indels", module: "Misalignment", module_order: 2, plot_type: PlotType::Bars, color: "#e50001", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Insertions", help: Some("Extra bases are present in the read but not in the contig") },
    VariableConfig { name: "deletions", subplot: "Indels", module: "Misalignment", module_order: 2, plot_type: PlotType::Bars, color: "#00E5E4", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Deletions", help: Some("A stretch of the contig has no corresponding bases in the read") },
    VariableConfig { name: "mismatches", subplot: "Mismatches", module: "Misalignment", module_order: 3, plot_type: PlotType::Bars, color: "#5a0f0b", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mismatches", help: None },
    
    // Per read metrics (long reads)
    VariableConfig { name: "read_lengths", subplot: "Read lengths", module: "Long-reads", module_order: 1, plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Read lengths", help: None },

    // Per read metrics (paired-reads)  
    VariableConfig { name: "insert_sizes", subplot: "Insert sizes", module: "Paired-reads", module_order: 1, plot_type: PlotType::Curve, color: "#ed8b00", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Insert sizes", help: None },
    VariableConfig { name: "non_inward_pairs", subplot: "Non-inward pairs", module: "Paired-reads", module_order: 2, plot_type: PlotType::Curve, color: "#c94009", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Non-inward pairs", help: None },
    VariableConfig { name: "mate_not_mapped", subplot: "Missing mates", module: "Paired-reads", module_order: 3, plot_type: PlotType::Curve, color: "#302DD2", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mate unmapped", help: None },
    VariableConfig { name: "mate_on_another_contig", subplot: "Missing mates", module: "Paired-reads", module_order: 4, plot_type: PlotType::Curve, color: "#CFD22D", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Mate mapped on another contig", help: None },

    // Phage termini module - for detecting phage DNA packaging sites
    VariableConfig { name: "coverage_reduced", subplot: "Coverage reduced", module: "Phage termini", module_order: 1, plot_type: PlotType::Curve, color: "#00c53b", alpha: 0.8, fill_alpha: 0.4, size: 1.0, title: "Coverage reduced", help: Some("Retains primary mappings starting with an exact match (and ending with an exact match for long reads)") },
    VariableConfig { name: "reads_starts", subplot: "Reads termini", module: "Phage termini", module_order: 2, plot_type: PlotType::Bars, color: "#215732", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read starts", help: None },
    VariableConfig { name: "reads_ends", subplot: "Reads termini", module: "Phage termini", module_order: 3, plot_type: PlotType::Bars, color: "#572146", alpha: 0.6, fill_alpha: 0.4, size: 1.0, title: "Read ends", help: None },
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
/// - `left_clippings`: Count of soft/hard clips at read starts (with Mean/Median/Std columns)
/// - `right_clippings`: Count of soft/hard clips at read ends (with Mean/Median/Std columns)
/// - `insertions`: Count of insertion events from CIGAR (with Mean/Median/Std columns)
/// - `deletions`: Deletion events from CIGAR
/// - `mismatches`: Base mismatches from MD tag
/// - `non_inward_pairs`: Mates on same contig but not properly oriented (for paired reads)
/// - `mate_not_mapped`: Mate is unmapped (for paired reads)
/// - `mate_on_another_contig`: Mate mapped to different contig (for paired reads)
pub const ASSEMBLYCHECK_FEATURES: &[&str] = &[
    "left_clippings",
    "right_clippings",
    "insertions",
    "deletions",
    "mismatches",
    "non_inward_pairs",
    "mate_not_mapped",
    "mate_on_another_contig",
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

/// Map feature name to a 1-based feature_id for the Feature_blob table.
/// Returns None if the feature is not found in VARIABLES.
#[inline]
pub fn feature_name_to_id(name: &str) -> Option<i16> {
    VARIABLES.iter().position(|v| v.name == name).map(|i| (i + 1) as i16)
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
    /// Sequence data for GC content computation (optional, may be large)
    pub sequence: Option<Vec<u8>>,
}

/// A gene/feature annotation from a GenBank or GFF3 file.
///
/// # Bioinformatics context:
/// Represents a CDS (coding sequence), gene, or other annotated feature
/// on a contig. Structural fields (position, strand, type) are stored directly;
/// all qualifiers (product, function, locus_tag, etc.) go into a generic HashMap
/// so that every qualifier from the source file is preserved.
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
    /// All qualifiers from the source file (key-value pairs).
    /// GenBank: product, function, locus_tag, phrog, protein_id, gene, note, db_xref, ...
    /// GFF3: product, Name, ID, Parent, locus_tag, ...
    pub qualifiers: HashMap<String, String>,
    /// Nucleotide sequence for CDS features (extracted from contig sequence)
    pub nucleotide_sequence: Option<String>,
    /// Protein sequence for CDS features (translated from nucleotide_sequence)
    pub protein_sequence: Option<String>,
}

/// Coverage data for a contig in a sample.
///
/// Stores coverage metrics used to assess read depth and abundance.
#[derive(Clone, Debug)]
pub struct PresenceData {
    /// The contig name
    pub contig_name: String,
    /// Percentage of bases with at least 1x coverage (0-100)
    pub coverage_pct: f32,
    /// Whether aligned fraction is above expected: exact_af >= (1 - e^(-0.883 × coverage_mean)) × 100
    pub above_expected_aligned_fraction: bool,
    /// Number of primary reads mapped to this contig
    pub read_count: u64,
    /// Coverage relative coverage roughness (Fano factor style): 1/(n-1) * Σ(cov(i+1) - cov(i))² / mean_cov
    /// Measures coverage smoothness normalized by mean (low = uniform, high = variable)
    pub coverage_relative_coverage_roughness: f32,
    /// Coverage coefficient of variation: Coefficient of Variation (CV) = std_dev / mean
    /// Measures overall spread from the mean, normalized to remove correlation with coverage depth
    /// Scaled by 1,000,000 for integer storage (same as coverage_relative_coverage_roughness)
    pub coverage_coefficient_of_variation: f32,
    /// Mean coverage depth across all positions
    pub coverage_mean: f32,
    /// Median coverage depth across all positions
    pub coverage_median: f32,
    /// Trimmed mean coverage (5% trim both sides)
    pub coverage_trimmed_mean: f32,
}

/// A terminus area with full metadata for database storage.
#[derive(Clone, Debug)]
pub struct TerminusArea {
    /// First position of the area (1-indexed)
    pub start_pos: i32,
    /// Last position of the area (1-indexed)
    pub end_pos: i32,
    /// Position with highest SPC (1-indexed)
    pub center_pos: i32,
    /// Total SPC (sum of read starts/ends in area)
    pub total_spc: u32,
    /// Coverage at center position
    pub coverage: u32,
    /// Tau value for the area
    pub tau: f64,
    /// Number of positions merged into this area
    pub number_peaks: u32,
    /// Poisson test p-value (peak significance)
    pub pvalue: f64,
    /// Bonferroni-adjusted p-value
    pub adjusted_pvalue: f64,
    /// Whether this terminus passed the Poisson significance test
    pub passed_poisson_test: bool,
    /// Whether this terminus passed the clipping excess test
    pub passed_clipping_test: bool,
    /// Sum of pre-filtered clippings used for filtering
    pub sum_clippings: u64,
    /// Global clipped ratio for this contig/sample
    pub clipped_ratio: f64,
    /// Expected clippings (threshold from statistical test)
    pub expected_clippings: f64,
    /// Expected SPC from Poisson test (lambda_w)
    pub expected_spc: f64,
    /// DTR-aware effective size of the area
    pub size: i32,
    /// Median clipping length within this area (for diagnostics)
    pub median_clippings: f64,
}

/// Phage packaging mechanism data for a contig in a sample.
///
/// Stores the detected packaging mechanism and terminus areas.
/// This is stored separately from PresenceData.
#[derive(Clone, Debug)]
pub struct PackagingData {
    /// The contig name
    pub contig_name: String,
    /// Phage packaging mechanism (e.g., "PAC", "COS", "DTR_short_5'", etc.)
    pub mechanism: String,
    /// Left terminus areas - start peaks with full metadata
    pub left_termini: Vec<TerminusArea>,
    /// Right terminus areas - end peaks with full metadata
    pub right_termini: Vec<TerminusArea>,
    /// Whether termini are in a duplication: Some(true) = DTR, Some(false) = ITR, None = no/mixed
    pub duplication: Option<bool>,
    /// Distance between start and end peak centers (when exactly 1 of each), None otherwise
    pub repeat_length: Option<i32>,
    /// Median clipping lengths for kept left termini (one per kept terminus)
    pub median_left_termini_clippings: Vec<f64>,
    /// Median clipping lengths for kept right termini (one per kept terminus)
    pub median_right_termini_clippings: Vec<f64>,
}

/// Type alias for feature calculation results.
///
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
