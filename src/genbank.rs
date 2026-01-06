//! GenBank file parsing functions.
//!
//! # Overview
//!
//! This module parses GenBank files to extract:
//! - Contig information (name, length)
//! - Gene annotations (start, end, strand, product, function, phrog)
//!
//! GenBank is a standard flat-file format for annotated nucleotide sequences.
//! Each record contains a sequence with associated features like CDS, tRNA, etc.
//!
//! # Python Equivalent
//!
//! This module corresponds to the GenBank parsing in `calculating_data.py`:
//! ```python
//! # Python: Uses BioPython's SeqIO
//! from Bio import SeqIO
//!
//! for rec in SeqIO.parse(args.genbank, "genbank"):
//!     contig_name = rec.name
//!     contig_length = len(rec.seq)
//!
//!     for f in rec.features:
//!         start = int(f.location.start) + 1
//!         end = int(f.location.end)
//!         strand_val = f.location.strand
//!         ftype = f.type
//!         # ... extract qualifiers
//! ```
//!
//! # Rust Implementation
//!
//! Uses the `gb_io` crate for GenBank parsing, which is pure Rust and
//! faster than BioPython for large files.

use anyhow::{Context, Result, anyhow};
use gb_io::reader::SeqReader;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::types::{ContigInfo, FeatureAnnotation};

/// Check if location is complemented (reverse strand).
fn is_complement(loc: &gb_io::seq::Location) -> bool {
    matches!(loc, gb_io::seq::Location::Complement(_))
}

/// Print context lines around a problematic area in the file
fn print_file_context(path: &Path, record_num: usize) {
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
            eprintln!("=== Context around record {} (starting at line {}) ===", record_num, record_start_line);
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
pub fn parse_genbank(path: &Path, annotation_tool: &str) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
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
                eprintln!("ERROR parsing GenBank record #{} in file: {}", record_num, path.display());
                eprintln!("gb_io error: {:?}", e);
                print_file_context(path, record_num);
                return Err(anyhow!("Failed to parse GenBank record #{}: {}", record_num, e));
            }
        };

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
            let phrog = feature.qualifier_values("phrog").next().and_then(|s| s.parse::<i32>().ok());

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
