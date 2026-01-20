//! Annotation file parsing functions.
//!
//! # Overview
//!
//! This module parses annotation files in multiple formats to extract:
//! - Contig information (name, length)
//! - Gene annotations (start, end, strand, product, function, phrog)
//!
//! # Supported Formats
//!
//! - **GenBank** (`.gbk`, `.gbff`): Standard NCBI flat-file format
//! - **GFF3** (`.gff`, `.gff3`): Generic Feature Format version 3
//!
//! # Usage
//!
//! Use the main entry point `parse_annotations()` which auto-detects format:
//! ```rust
//! let (contigs, annotations) = parse_annotations(path, annotation_tool)?;
//! ```

use anyhow::{anyhow, Context, Result};
use gb_io::reader::SeqReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::types::{ContigInfo, FeatureAnnotation};

// ============================================================================
// Format Detection
// ============================================================================

/// Supported annotation file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    /// GenBank format (.gbk, .gbff, .gb)
    GenBank,
    /// GFF3 format (.gff, .gff3)
    Gff3,
}

/// Detect the annotation format from file extension.
fn detect_format(path: &Path) -> Result<AnnotationFormat> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "gbk" | "gbff" | "gb" | "genbank" => Ok(AnnotationFormat::GenBank),
        "gff" | "gff3" => Ok(AnnotationFormat::Gff3),
        _ => Err(anyhow!(
            "Unknown annotation file format: '{}'. Supported extensions: .gbk, .gbff, .gb, .gff, .gff3",
            ext
        )),
    }
}

// ============================================================================
// Main Entry Point
// ============================================================================

/// Parse annotation file and extract contigs and annotations.
///
/// This is the main entry point that auto-detects file format based on extension.
///
/// # Arguments
/// * `path` - Path to the annotation file (.gbk, .gbff, .gff, .gff3)
/// * `annotation_tool` - Name of the annotation tool (e.g., "pharokka", "bakta")
///
/// # Returns
/// A tuple of (contigs, annotations) where:
/// - contigs: Vector of ContigInfo with name, length, and annotation_tool
/// - annotations: Vector of FeatureAnnotation with feature details
pub fn parse_annotations(
    path: &Path,
    annotation_tool: &str,
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
    let format = detect_format(path)?;

    match format {
        AnnotationFormat::GenBank => parse_genbank(path, annotation_tool),
        AnnotationFormat::Gff3 => parse_gff3(path, annotation_tool),
    }
}

// ============================================================================
// GenBank Parser
// ============================================================================

/// Check if location is complemented (reverse strand).
fn is_complement(loc: &gb_io::seq::Location) -> bool {
    matches!(loc, gb_io::seq::Location::Complement(_))
}

/// Print context lines around a problematic area in the file
fn print_genbank_context(path: &Path, record_num: usize) {
    // Try to find the start of the problematic record by counting LOCUS lines
    if let Ok(file) = File::open(path) {
        let reader = BufReader::new(file);
        let mut locus_count = 0;
        let mut line_num = 0;
        let mut record_start_line = 0;

        for line in reader.lines() {
            line_num += 1;
            if let Ok(ref l) = line {
                if l.starts_with("LOCUS") {
                    locus_count += 1;
                    if locus_count == record_num {
                        record_start_line = line_num;
                    }
                }
            }
        }

        // Now print lines around the problematic record
        if record_start_line > 0 {
            eprintln!(
                "=== Context around record {} (starting at line {}) ===",
                record_num, record_start_line
            );
            if let Ok(file2) = File::open(path) {
                let reader2 = BufReader::new(file2);
                for (i, line) in reader2.lines().enumerate() {
                    let ln = i + 1;
                    // Print 50 lines from the start of the problematic record
                    if ln >= record_start_line && ln < record_start_line + 50 {
                        if let Ok(l) = line {
                            eprintln!("{:6}: {}", ln, l);
                        }
                    }
                }
            }
            eprintln!("=== End context ===");
        }
    }
}

/// Parse GenBank file and extract contigs and annotations.
///
/// Handles `.gbk`, `.gbff`, and `.gb` files using the gb_io crate.
pub fn parse_genbank(
    path: &Path,
    annotation_tool: &str,
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
    let file = File::open(path).context("Failed to open GenBank file")?;
    let reader = SeqReader::new(file);

    let mut contigs = Vec::new();
    let mut annotations = Vec::new();
    let mut contig_id = 1i64;
    let mut record_num = 0usize;

    for result in reader {
        record_num += 1;
        let seq = match result {
            Ok(s) => s,
            Err(e) => {
                eprintln!(
                    "ERROR parsing GenBank record #{} in file: {}",
                    record_num,
                    path.display()
                );
                eprintln!("gb_io error: {:?}", e);
                print_genbank_context(path, record_num);
                return Err(anyhow!(
                    "Failed to parse GenBank record #{}: {}",
                    record_num,
                    e
                ));
            }
        };

        let name = seq
            .name
            .clone()
            .unwrap_or_else(|| format!("contig_{}", contig_id));
        let length = seq.seq.len();

        contigs.push(ContigInfo {
            name: name.clone(),
            length,
            annotation_tool: annotation_tool.to_string(),
        });

        // Extract features
        for feature in &seq.features {
            let feature_type = feature.kind.to_string();
            if feature_type == "source" || feature_type == "gene" {
                continue;
            }

            // Use find_bounds() for start/end (returns 0-based exclusive range)
            let (start, end) = match feature.location.find_bounds() {
                Ok((s, e)) => (s + 1, e), // Convert to 1-based
                Err(_) => continue,
            };

            // Check if complement (reverse strand)
            let strand = if is_complement(&feature.location) {
                -1
            } else {
                1
            };

            let product = feature
                .qualifier_values("product")
                .next()
                .map(|s| s.to_string());
            let function = feature
                .qualifier_values("function")
                .next()
                .map(|s| s.to_string());
            let phrog = feature
                .qualifier_values("phrog")
                .next()
                .and_then(|s| s.parse::<i32>().ok());

            annotations.push(FeatureAnnotation {
                contig_id,
                start,
                end,
                strand,
                feature_type,
                product,
                function,
                phrog,
            });
        }

        contig_id += 1;
    }

    Ok((contigs, annotations))
}

// ============================================================================
// GFF3 Parser
// ============================================================================

/// Parse GFF3 attributes column (column 9) into key-value pairs.
///
/// GFF3 attributes are semicolon-separated, with key=value format.
/// Values may be URL-encoded.
fn parse_gff3_attributes(attrs_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();

    for attr in attrs_str.split(';') {
        let attr = attr.trim();
        if attr.is_empty() {
            continue;
        }

        if let Some((key, value)) = attr.split_once('=') {
            // URL-decode common escapes
            let value = value
                .replace("%3B", ";")
                .replace("%3D", "=")
                .replace("%26", "&")
                .replace("%2C", ",")
                .replace("%09", "\t")
                .replace("%0A", "\n")
                .replace("%25", "%");
            attrs.insert(key.to_string(), value);
        }
    }

    attrs
}

/// Parse GFF3 file and extract contigs and annotations.
///
/// GFF3 format has 9 tab-separated columns:
/// 1. seqid - contig/chromosome name
/// 2. source - annotation source
/// 3. type - feature type (CDS, gene, etc.)
/// 4. start - 1-based start position
/// 5. end - 1-based end position (inclusive)
/// 6. score - score (often ".")
/// 7. strand - + or - (or . for unknown)
/// 8. phase - 0, 1, 2 for CDS (or .)
/// 9. attributes - key=value pairs separated by semicolons
///
/// Contig lengths are obtained from:
/// 1. `##sequence-region` directives (preferred)
/// 2. Maximum feature end position (fallback)
pub fn parse_gff3(
    path: &Path,
    annotation_tool: &str,
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
    let file = File::open(path).context("Failed to open GFF file")?;
    let reader = BufReader::new(file);

    // First pass: collect sequence regions and track max positions per contig
    let mut sequence_regions: HashMap<String, usize> = HashMap::new();
    let mut max_positions: HashMap<String, usize> = HashMap::new();
    let mut features_data: Vec<(String, String, i64, i64, i64, HashMap<String, String>)> =
        Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.context(format!("Failed to read line {}", line_num + 1))?;
        let line = line.trim();

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        // Handle directives
        if line.starts_with("##") {
            // Parse ##sequence-region seqid start end
            if line.starts_with("##sequence-region") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 4 {
                    let seqid = parts[1].to_string();
                    if let Ok(end) = parts[3].parse::<usize>() {
                        sequence_regions.insert(seqid, end);
                    }
                }
            }
            continue;
        }

        // Skip comment lines
        if line.starts_with('#') {
            continue;
        }

        // Parse feature line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue; // Invalid line, skip
        }

        let seqid = fields[0].to_string();
        let feature_type = fields[2].to_string();

        // Skip "source", "gene", and "region" features (like GenBank parser)
        if feature_type == "source" || feature_type == "gene" || feature_type == "region" {
            // But still track position for contig length
            if let (Ok(start), Ok(end)) = (fields[3].parse::<usize>(), fields[4].parse::<usize>()) {
                let current_max = max_positions.entry(seqid.clone()).or_insert(0);
                *current_max = (*current_max).max(end).max(start);
            }
            continue;
        }

        // Parse positions
        let start: i64 = fields[3]
            .parse()
            .context(format!("Invalid start position on line {}", line_num + 1))?;
        let end: i64 = fields[4]
            .parse()
            .context(format!("Invalid end position on line {}", line_num + 1))?;

        // Track max position for contig length fallback
        let current_max = max_positions.entry(seqid.clone()).or_insert(0);
        *current_max = (*current_max).max(end as usize).max(start as usize);

        // Parse strand
        let strand = match fields[6] {
            "+" => 1i64,
            "-" => -1i64,
            _ => 1i64, // Default to forward strand
        };

        // Parse attributes
        let attrs = parse_gff3_attributes(fields[8]);

        features_data.push((seqid, feature_type, start, end, strand, attrs));
    }

    // Build contig list with proper ordering
    let mut contig_names: Vec<String> = sequence_regions
        .keys()
        .chain(max_positions.keys())
        .cloned()
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    contig_names.sort(); // Consistent ordering

    // Create contig ID mapping
    let mut contig_id_map: HashMap<String, i64> = HashMap::new();
    let mut contigs = Vec::new();

    for (idx, name) in contig_names.iter().enumerate() {
        let contig_id = (idx + 1) as i64;
        contig_id_map.insert(name.clone(), contig_id);

        // Get length from sequence-region directive, or fallback to max position
        let length = sequence_regions
            .get(name)
            .copied()
            .or_else(|| max_positions.get(name).copied())
            .unwrap_or(0);

        contigs.push(ContigInfo {
            name: name.clone(),
            length,
            annotation_tool: annotation_tool.to_string(),
        });
    }

    // Build annotations
    let mut annotations = Vec::new();

    for (seqid, feature_type, start, end, strand, attrs) in features_data {
        let contig_id = *contig_id_map.get(&seqid).unwrap_or(&1);

        // Extract product - try multiple attribute names used by different tools
        let product = attrs
            .get("product")
            .or_else(|| attrs.get("Product"))
            .or_else(|| attrs.get("Name"))
            .or_else(|| attrs.get("name"))
            .or_else(|| attrs.get("description"))
            .cloned();

        // Extract function
        let function = attrs
            .get("function")
            .or_else(|| attrs.get("Function"))
            .or_else(|| attrs.get("note"))
            .or_else(|| attrs.get("Note"))
            .cloned();

        // Extract phrog ID (pharokka-specific)
        let phrog = attrs
            .get("phrog")
            .or_else(|| attrs.get("PHROG"))
            .and_then(|s| s.parse::<i32>().ok());

        annotations.push(FeatureAnnotation {
            contig_id,
            start,
            end,
            strand,
            feature_type,
            product,
            function,
            phrog,
        });
    }

    if contigs.is_empty() {
        return Err(anyhow!(
            "No contigs found in GFF file: {}. Make sure the file contains valid GFF3 data.",
            path.display()
        ));
    }

    Ok((contigs, annotations))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gff3_attributes() {
        let attrs = parse_gff3_attributes("ID=cds1;Name=test;product=hypothetical protein");
        assert_eq!(attrs.get("ID"), Some(&"cds1".to_string()));
        assert_eq!(attrs.get("Name"), Some(&"test".to_string()));
        assert_eq!(
            attrs.get("product"),
            Some(&"hypothetical protein".to_string())
        );
    }

    #[test]
    fn test_parse_gff3_attributes_with_escapes() {
        let attrs = parse_gff3_attributes("Name=test%3Bvalue;product=a%3Db");
        assert_eq!(attrs.get("Name"), Some(&"test;value".to_string()));
        assert_eq!(attrs.get("product"), Some(&"a=b".to_string()));
    }

    #[test]
    fn test_detect_format() {
        assert_eq!(
            detect_format(Path::new("test.gbk")).unwrap(),
            AnnotationFormat::GenBank
        );
        assert_eq!(
            detect_format(Path::new("test.gbff")).unwrap(),
            AnnotationFormat::GenBank
        );
        assert_eq!(
            detect_format(Path::new("test.gff")).unwrap(),
            AnnotationFormat::Gff3
        );
        assert_eq!(
            detect_format(Path::new("test.gff3")).unwrap(),
            AnnotationFormat::Gff3
        );
        assert!(detect_format(Path::new("test.txt")).is_err());
    }
}
