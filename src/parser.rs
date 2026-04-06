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
//! let (contigs, annotations) = parse_annotations(path)?;
//! ```
//!
//! Pharokka-specific fields (`phrog`, `function` categories) are extracted
//! automatically when present — no annotation tool flag needed.

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

/// Parse annotation file and extract contigs, annotations, and contig-level qualifiers.
///
/// This is the main entry point that auto-detects file format based on extension.
///
/// # Arguments
/// * `path` - Path to the annotation file (.gbk, .gbff, .gff, .gff3)
///
/// # Returns
/// A tuple of (contigs, annotations, contig_qualifiers) where:
/// - contigs: Vector of ContigInfo with name and length
/// - annotations: Vector of FeatureAnnotation with feature details
/// - contig_qualifiers: Vector of (contig_id, HashMap) from source features
pub fn parse_annotations(
    path: &Path,
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>, Vec<(i64, HashMap<String, String>)>)> {
    let format = detect_format(path)?;

    match format {
        AnnotationFormat::GenBank => parse_genbank(path),
        AnnotationFormat::Gff3 => parse_gff3(path),
    }
}

// ============================================================================
// GenBank Parser
// ============================================================================

/// Check if location is complemented (reverse strand).
fn is_complement(loc: &gb_io::seq::Location) -> bool {
    matches!(loc, gb_io::seq::Location::Complement(_))
}

/// Recursively walk a gb-io Location, flattening Join/Order parts and
/// stripping Complement wrappers, to produce a list of 1-based inclusive
/// (start, end) sub-intervals in the order they appear in the source file.
fn collect_ranges(loc: &gb_io::seq::Location, out: &mut Vec<(i64, i64)>) {
    use gb_io::seq::Location::*;
    match loc {
        // gb-io stores ranges 0-based half-open; convert to 1-based inclusive.
        Range((a, _), (b, _)) => out.push((*a + 1, *b)),
        Complement(inner) => collect_ranges(inner, out),
        Join(parts) | Order(parts) => {
            for p in parts {
                collect_ranges(p, out);
            }
        }
        // Bond, OneOf, External, Gap, etc. — not representable as a genomic
        // interval on this contig; silently skipped.
        _ => {}
    }
}

/// Extract (bounding_box, segments, strand) from a gb-io Location.
///
/// - `bounding_box` is `(min_start, max_end)` over all sub-intervals, 1-based inclusive.
/// - `segments` is the full list of sub-intervals when the location has more
///   than one part; empty for a single Range (no need for Annotation_segments rows).
/// - `strand` is -1 when the outermost wrapper is `Complement`, else 1.
fn extract_location(
    loc: &gb_io::seq::Location,
) -> Option<((i64, i64), Vec<(i64, i64)>, i64)> {
    let strand = if is_complement(loc) { -1i64 } else { 1i64 };
    let mut parts: Vec<(i64, i64)> = Vec::new();
    collect_ranges(loc, &mut parts);
    if parts.is_empty() {
        return None;
    }
    let bbox_start = parts.iter().map(|&(s, _)| s).min().unwrap();
    let bbox_end = parts.iter().map(|&(_, e)| e).max().unwrap();
    let segments = if parts.len() >= 2 { parts } else { Vec::new() };
    Some(((bbox_start, bbox_end), segments, strand))
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
///
/// Returns (contigs, annotations, contig_qualifiers) where contig_qualifiers
/// are (contig_id, HashMap) pairs extracted from "source" features.
pub fn parse_genbank(
    path: &Path,
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>, Vec<(i64, HashMap<String, String>)>)> {
    let file = File::open(path).context("Failed to open GenBank file")?;
    let reader = SeqReader::new(file);

    let mut contigs = Vec::new();
    let mut annotations = Vec::new();
    let mut contig_qualifiers: Vec<(i64, HashMap<String, String>)> = Vec::new();
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
            .version
            .clone()
            .or_else(|| seq.name.clone())
            .unwrap_or_else(|| format!("contig_{}", contig_id));
        let length = seq.seq.len();

        // Extract sequence bytes for GC content computation
        let sequence = if seq.seq.is_empty() {
            None
        } else {
            Some(seq.seq.clone())
        };

        contigs.push(ContigInfo {
            name: name.clone(),
            length,
            sequence,
        });

        // Extract features
        let mut mrna_counter: usize = 0;
        let contig_feature_start = annotations.len();
        for feature in &seq.features {
            let feature_type = feature.kind.to_string();

            // Skip "source" — its qualifiers are stored as contig-level metadata
            if feature_type == "source" {
                // Collect source qualifiers for this contig
                let source_quals = collect_qualifiers(&feature.qualifiers);
                if !source_quals.is_empty() {
                    contig_qualifiers.push((contig_id, source_quals));
                }
                continue;
            }

            // Walk the Location tree: bbox + per-interval segments + strand.
            // Single-Range features yield empty `segments`; Join/Order features
            // carry the sub-intervals so we can store them in Annotation_segments.
            let ((start, end), segments, strand) = match extract_location(&feature.location) {
                Some(v) => v,
                None => continue,
            };

            // Collect all qualifiers from the feature
            let qualifiers = collect_qualifiers(&feature.qualifiers);

            // Assign a synthetic self_key for mRNA features so that CDS
            // containment matching (resolve_genbank_parents) can point at them.
            let self_key = if feature_type == "mRNA" {
                mrna_counter += 1;
                Some(format!("gbk_mrna_{}_{}", contig_id, mrna_counter))
            } else {
                None
            };

            annotations.push(FeatureAnnotation {
                contig_id,
                start,
                end,
                strand,
                feature_type,
                qualifiers,
                nucleotide_sequence: None,
                protein_sequence: None,
                segments,
                parent_key: None,
                self_key,
            });
        }

        // For GenBank: resolve CDS → mRNA relationships within this contig by
        // segment containment + shared locus_tag/gene. Runs once per contig so
        // the matching window is small.
        resolve_genbank_parents(&mut annotations[contig_feature_start..]);

        contig_id += 1;
    }

    Ok((contigs, annotations, contig_qualifiers))
}

/// Resolve CDS → mRNA parent links within a single GenBank contig's features.
///
/// GenBank has no explicit `Parent=` attribute, so we match each CDS against
/// the mRNAs sharing the same locus_tag (or gene as fallback) on the same
/// strand, keeping only the mRNA whose segments fully cover the CDS's
/// segments. If exactly one mRNA matches, we set the CDS's `parent_key` to
/// that mRNA's `self_key`; zero or multiple matches leave `parent_key = None`
/// and the CDS falls through to the per-locus first-occurrence rule.
fn resolve_genbank_parents(features: &mut [FeatureAnnotation]) {
    // Collect an immutable snapshot of (index, locus_key, self_key, strand, segments or bbox) for mRNAs.
    let mrnas: Vec<(usize, String, String, i64, Vec<(i64, i64)>)> = features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.feature_type == "mRNA")
        .filter_map(|(i, f)| {
            let key = f
                .qualifiers
                .get("locus_tag")
                .or_else(|| f.qualifiers.get("gene"))
                .cloned()?;
            let sk = f.self_key.clone()?;
            // Use segments when present; otherwise fall back to a single bbox interval
            // (unspliced mRNA — rare but handled).
            let segs = if f.segments.is_empty() {
                vec![(f.start, f.end)]
            } else {
                f.segments.clone()
            };
            Some((i, key, sk, f.strand, segs))
        })
        .collect();

    if mrnas.is_empty() {
        return;
    }

    for cds in features.iter_mut() {
        if cds.feature_type != "CDS" {
            continue;
        }
        let cds_key = match cds
            .qualifiers
            .get("locus_tag")
            .or_else(|| cds.qualifiers.get("gene"))
        {
            Some(v) => v.clone(),
            None => continue,
        };
        let cds_segs: Vec<(i64, i64)> = if cds.segments.is_empty() {
            vec![(cds.start, cds.end)]
        } else {
            cds.segments.clone()
        };

        // A match requires: same locus_key, same strand, and every CDS segment
        // is fully contained in some mRNA segment.
        let mut matched: Option<&str> = None;
        let mut ambiguous = false;
        for (_, mrna_key, mrna_sk, mrna_strand, mrna_segs) in &mrnas {
            if *mrna_key != cds_key || *mrna_strand != cds.strand {
                continue;
            }
            let all_contained = cds_segs.iter().all(|&(cs, ce)| {
                mrna_segs.iter().any(|&(ms, me)| ms <= cs && me >= ce)
            });
            if !all_contained {
                continue;
            }
            if matched.is_some() {
                ambiguous = true;
                break;
            }
            matched = Some(mrna_sk.as_str());
        }
        if !ambiguous {
            if let Some(sk) = matched {
                cds.parent_key = Some(sk.to_string());
            }
        }
    }
}

/// Collect all qualifiers from a GenBank feature into a HashMap.
///
/// When a qualifier key appears multiple times (e.g., multiple db_xref entries),
/// values are joined with "; ".
fn collect_qualifiers(qualifiers: &[(std::borrow::Cow<'static, str>, Option<String>)]) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for (key, value) in qualifiers {
        if let Some(val) = value {
            let key_str = key.to_string();
            map.entry(key_str)
                .and_modify(|existing: &mut String| {
                    existing.push_str("; ");
                    existing.push_str(val);
                })
                .or_insert_with(|| val.clone());
        }
    }
    map
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

/// Normalize GFF3 attribute aliases into canonical GenBank qualifier names.
/// Ensures the Contig_annotation view works consistently regardless of format.
fn normalize_gff_qualifiers(qualifiers: &mut HashMap<String, String>) {
    // product: GFF3 tools use Name, name, description, or product
    if !qualifiers.contains_key("product") {
        if let Some(val) = qualifiers
            .get("Product")
            .or_else(|| qualifiers.get("Name"))
            .or_else(|| qualifiers.get("name"))
            .or_else(|| qualifiers.get("description"))
            .cloned()
        {
            qualifiers.insert("product".to_string(), val);
        }
    }
    // function: GFF3 tools use note/Note as a fallback
    if !qualifiers.contains_key("function") {
        if let Some(val) = qualifiers
            .get("Function")
            .or_else(|| qualifiers.get("note"))
            .or_else(|| qualifiers.get("Note"))
            .cloned()
        {
            qualifiers.insert("function".to_string(), val);
        }
    }
    // phrog: normalize casing
    if !qualifiers.contains_key("phrog") {
        if let Some(val) = qualifiers.get("PHROG").cloned() {
            qualifiers.insert("phrog".to_string(), val);
        }
    }
    // locus_tag: normalize hyphen variant
    if !qualifiers.contains_key("locus_tag") {
        if let Some(val) = qualifiers.get("locus-tag").cloned() {
            qualifiers.insert("locus_tag".to_string(), val);
        }
    }
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
) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>, Vec<(i64, HashMap<String, String>)>)> {
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
        // "source": Describes the entire sequence/contig—metadata rather than an actual coding feature
        // "gene": Often a parent/container feature in GFF3 that wraps CDS/mRNA children. You typically want just the actual coding regions
        // "region": A generic feature type often used for metadata or structural regions, not functional genes
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
            sequence: None, // GFF3 doesn't contain sequence data
        });
    }

    // Normalize qualifiers on every row so ID/Parent lookups and downstream
    // queries all see canonical keys (product, function, locus_tag, ...).
    for row in features_data.iter_mut() {
        normalize_gff_qualifiers(&mut row.5);
    }

    // Build the ID → index map for rows that expose an `ID=` attribute.
    let mut id_to_index: HashMap<String, usize> = HashMap::new();
    for (idx, row) in features_data.iter().enumerate() {
        if let Some(id) = row.5.get("ID") {
            id_to_index.entry(id.clone()).or_insert(idx);
        }
    }

    // Group child rows (exon / CDS) under their transcript-like parent and
    // mark them as absorbed so they aren't emitted as standalone annotations.
    let mut exon_segments_by_parent: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
    let mut cds_segments_by_parent: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
    // For the collapsed CDS row we need a strand (first child's) and a
    // qualifier overlay (first child's — typically has protein_id / phase).
    let mut cds_overlay_by_parent: HashMap<usize, (i64, HashMap<String, String>)> =
        HashMap::new();
    let mut absorbed: Vec<bool> = vec![false; features_data.len()];

    for idx in 0..features_data.len() {
        let child_type = features_data[idx].1.clone();
        if child_type != "exon" && child_type != "CDS" {
            continue;
        }
        let parent_raw = match features_data[idx].5.get("Parent").cloned() {
            Some(p) => p,
            None => continue,
        };
        // Parent= can be a comma-separated list; pick the first resolvable
        // transcript-like ancestor.
        let parent_idx_opt = parent_raw
            .split(',')
            .filter_map(|p| id_to_index.get(p.trim()).copied())
            .find(|&pi| {
                matches!(
                    features_data[pi].1.as_str(),
                    "mRNA" | "transcript" | "primary_transcript"
                )
            });
        let parent_idx = match parent_idx_opt {
            Some(i) => i,
            None => continue,
        };
        let (_, _, c_start, c_end, c_strand, c_attrs) = features_data[idx].clone();
        if child_type == "exon" {
            exon_segments_by_parent
                .entry(parent_idx)
                .or_default()
                .push((c_start, c_end));
        } else {
            cds_segments_by_parent
                .entry(parent_idx)
                .or_default()
                .push((c_start, c_end));
            cds_overlay_by_parent
                .entry(parent_idx)
                .or_insert((c_strand, c_attrs));
        }
        absorbed[idx] = true;
    }

    // Sort collected segments in genomic order once.
    for segs in exon_segments_by_parent.values_mut() {
        segs.sort_by_key(|&(s, _)| s);
    }
    for segs in cds_segments_by_parent.values_mut() {
        segs.sort_by_key(|&(s, _)| s);
    }

    // Build annotations: walk features_data in file order, skipping absorbed
    // rows. For each transcript-like parent that has exon/CDS children, fold
    // them into the transcript and synthesise a collapsed CDS sibling.
    let mut annotations = Vec::new();

    for (idx, (seqid, feature_type, start, end, strand, attrs)) in
        features_data.into_iter().enumerate()
    {
        if absorbed[idx] {
            continue;
        }
        let contig_id = *contig_id_map.get(&seqid).unwrap_or(&1);
        let qualifiers = attrs;

        let is_transcript = matches!(
            feature_type.as_str(),
            "mRNA" | "transcript" | "primary_transcript"
        );
        let self_id = qualifiers.get("ID").cloned();

        // Exon children → segments of the transcript itself.
        let exon_segs = if is_transcript {
            exon_segments_by_parent.remove(&idx).unwrap_or_default()
        } else {
            Vec::new()
        };
        // CDS children → a separate collapsed CDS annotation that references
        // this transcript as its parent.
        let cds_segs = if is_transcript {
            cds_segments_by_parent.remove(&idx).unwrap_or_default()
        } else {
            Vec::new()
        };
        let cds_overlay = if is_transcript {
            cds_overlay_by_parent.remove(&idx)
        } else {
            None
        };

        // Emit the transcript (or non-transcript) row itself. Its bounding box
        // is the row's own start/end; if exon children were present their
        // min/max matches the transcript's declared extent anyway.
        let transcript_segments = if exon_segs.len() >= 2 {
            exon_segs.clone()
        } else {
            Vec::new()
        };
        annotations.push(FeatureAnnotation {
            contig_id,
            start,
            end,
            strand,
            feature_type,
            qualifiers: qualifiers.clone(),
            nucleotide_sequence: None,
            protein_sequence: None,
            segments: transcript_segments,
            parent_key: None,
            self_key: if is_transcript { self_id.clone() } else { None },
        });

        // Synthesise the collapsed CDS sibling when the transcript carried
        // CDS children.
        if !cds_segs.is_empty() {
            let (cds_strand, cds_child_attrs) = cds_overlay.unwrap_or((strand, HashMap::new()));
            // Start from the transcript's qualifiers (they carry product,
            // gene, locus_tag) and overlay the first CDS child's attributes
            // (protein_id, phase, etc.).
            let mut cds_qualifiers = qualifiers.clone();
            for (k, v) in cds_child_attrs {
                cds_qualifiers.entry(k).or_insert(v);
            }
            let cds_start = cds_segs.iter().map(|&(s, _)| s).min().unwrap();
            let cds_end = cds_segs.iter().map(|&(_, e)| e).max().unwrap();
            let cds_segments_vec = if cds_segs.len() >= 2 {
                cds_segs
            } else {
                Vec::new()
            };
            annotations.push(FeatureAnnotation {
                contig_id,
                start: cds_start,
                end: cds_end,
                strand: cds_strand,
                feature_type: "CDS".to_string(),
                qualifiers: cds_qualifiers,
                nucleotide_sequence: None,
                protein_sequence: None,
                segments: cds_segments_vec,
                parent_key: self_id,
                self_key: None,
            });
        }
    }

    if contigs.is_empty() {
        return Err(anyhow!(
            "No contigs found in GFF file: {}. Make sure the file contains valid GFF3 data.",
            path.display()
        ));
    }

    // GFF3 has no "source" features, so contig_qualifiers is empty
    Ok((contigs, annotations, Vec::new()))
}

// ============================================================================
// FASTA Parser
// ============================================================================

/// Parse a FASTA file and return (name, sequence) pairs.
///
/// Reads `>name ...` headers (takes first whitespace-delimited word as name)
/// and concatenates sequence lines as bytes.
pub fn parse_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(path).with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.context("Failed to read FASTA line")?;
        let line = line.trim_end();

        if line.starts_with('>') {
            // Save previous record
            if let Some(name) = current_name.take() {
                if !current_seq.is_empty() {
                    records.push((name, current_seq.clone()));
                }
            }
            // Parse header: take first whitespace-delimited word after '>'
            let header = &line[1..];
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }

    // Save last record
    if let Some(name) = current_name {
        if !current_seq.is_empty() {
            records.push((name, current_seq));
        }
    }

    Ok(records)
}

// ============================================================================
// Sequence Translation
// ============================================================================

/// Compute the reverse complement of a DNA sequence.
fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Translate a DNA sequence to protein using the standard genetic code.
/// Handles partial codons at the end by ignoring them.
fn translate_dna(dna: &[u8]) -> String {
    let mut protein = String::with_capacity(dna.len() / 3 + 1);
    for codon in dna.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        let aa = match &[codon[0].to_ascii_uppercase(), codon[1].to_ascii_uppercase(), codon[2].to_ascii_uppercase()] {
            b"TTT" | b"TTC" => 'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
            b"ATT" | b"ATC" | b"ATA" => 'I',
            b"ATG" => 'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
            b"TAT" | b"TAC" => 'Y',
            b"TAA" | b"TAG" | b"TGA" => '*',
            b"CAT" | b"CAC" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"AAT" | b"AAC" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"GAT" | b"GAC" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"TGT" | b"TGC" => 'C',
            b"TGG" => 'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
            _ => 'X',
        };
        protein.push(aa);
    }
    protein
}

/// Compute nucleotide and protein sequences for CDS annotations.
///
/// For each CDS annotation, extracts the subsequence from the contig sequence,
/// handles strand (reverse complement for -1), and translates to protein.
/// This must be called AFTER parse_annotations() AND merge_sequences_from_fasta()
/// so that contig sequences are available.
pub fn compute_annotation_sequences(
    contigs: &[ContigInfo],
    annotations: &mut [FeatureAnnotation],
) {
    // Build contig_id -> sequence lookup
    let seq_by_id: HashMap<i64, &[u8]> = contigs
        .iter()
        .enumerate()
        .filter_map(|(i, c)| c.sequence.as_ref().map(|s| ((i + 1) as i64, s.as_slice())))
        .collect();

    let mut count = 0;
    for ann in annotations.iter_mut() {
        if ann.feature_type != "CDS" {
            continue;
        }
        let Some(seq) = seq_by_id.get(&ann.contig_id) else {
            continue;
        };

        // Spliced CDS: concatenate exon segments in genomic order, then apply
        // reverse complement if on the minus strand. For an unspliced CDS the
        // segments list is empty and we fall back to the bounding box slice.
        let nuc: Vec<u8> = if !ann.segments.is_empty() {
            let mut buf = Vec::with_capacity(
                ann.segments.iter().map(|&(s, e)| (e - s + 1) as usize).sum::<usize>(),
            );
            let mut ok = true;
            for &(s, e) in &ann.segments {
                let s0 = (s - 1) as usize;
                let e0 = e as usize;
                if e0 > seq.len() || s0 >= e0 {
                    ok = false;
                    break;
                }
                buf.extend_from_slice(&seq[s0..e0]);
            }
            if !ok {
                continue;
            }
            buf
        } else {
            let start = (ann.start - 1) as usize;
            let end = ann.end as usize;
            if end > seq.len() || start >= end {
                continue;
            }
            seq[start..end].to_vec()
        };

        if nuc.is_empty() {
            continue;
        }
        let nuc = if ann.strand == -1 {
            reverse_complement(&nuc)
        } else {
            nuc
        };

        let nuc_str = String::from_utf8_lossy(&nuc).into_owned();
        let protein = translate_dna(&nuc);

        ann.nucleotide_sequence = Some(nuc_str);
        ann.protein_sequence = Some(protein);
        count += 1;
    }
    eprintln!("Computed sequences for {count} CDS annotations");
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
