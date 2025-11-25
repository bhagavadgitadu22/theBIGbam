//! GenBank file parsing functions.
//!
//! This module handles parsing GenBank files to extract contig and annotation information.

use anyhow::{Context, Result};
use gb_io::reader::SeqReader;
use std::fs::File;
use std::path::Path;

use crate::types::{ContigInfo, FeatureAnnotation};

/// Check if location is complemented (reverse strand).
fn is_complement(loc: &gb_io::seq::Location) -> bool {
    matches!(loc, gb_io::seq::Location::Complement(_))
}

/// Parse GenBank file and extract contigs and annotations.
pub fn parse_genbank(path: &Path, annotation_tool: &str) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
    let file = File::open(path).context("Failed to open GenBank file")?;
    let reader = SeqReader::new(file);

    let mut contigs = Vec::new();
    let mut annotations = Vec::new();
    let mut contig_id = 1i64;

    for result in reader {
        let seq = result.context("Failed to parse GenBank record")?;

        let name = seq.name.clone().unwrap_or_else(|| format!("contig_{}", contig_id));
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
                Ok((s, e)) => (s + 1, e),  // Convert to 1-based
                Err(_) => continue,
            };

            // Check if complement (reverse strand)
            let strand = if is_complement(&feature.location) { -1 } else { 1 };

            let product = feature.qualifier_values("product").next().map(|s| s.to_string());
            let function = feature.qualifier_values("function").next().map(|s| s.to_string());
            let phrog = feature.qualifier_values("phrog").next().map(|s| s.to_string());

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
