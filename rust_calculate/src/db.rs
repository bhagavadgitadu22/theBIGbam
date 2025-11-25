//! SQLite database operations.
//!
//! This module handles creating and updating the metadata SQLite database.

use anyhow::Result;
use rusqlite::Connection;
use std::collections::HashMap;
use std::path::Path;

use crate::types::{ContigInfo, FeatureAnnotation, PresenceData};

/// Create metadata SQLite database with schema and initial data.
pub fn create_metadata_db(
    db_path: &Path,
    contigs: &[ContigInfo],
    annotations: &[FeatureAnnotation],
) -> Result<()> {
    let conn = Connection::open(db_path)?;

    // Create tables
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Annotation_tool TEXT
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Sample_name TEXT UNIQUE
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Presences (
            Presence_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Sample_id INTEGER,
            Coverage_percentage REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Contig_annotation (
            Contig_annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Start INTEGER,
            End INTEGER,
            Strand INTEGER,
            Type TEXT,
            Product TEXT,
            Function TEXT,
            Phrog INTEGER,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id)
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Variable (
            Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Variable_name TEXT UNIQUE,
            Subplot TEXT,
            Module TEXT,
            Type TEXT,
            Color TEXT,
            Alpha REAL,
            Fill_alpha REAL,
            Size REAL,
            Title TEXT,
            Help TEXT,
            Feature_table_name TEXT
        )",
        [],
    )?;

    // Insert contigs
    for contig in contigs {
        conn.execute(
            "INSERT INTO Contig (Contig_name, Contig_length, Annotation_tool) VALUES (?1, ?2, ?3)",
            [&contig.name, &contig.length.to_string(), &contig.annotation_tool],
        )?;
    }

    // Insert annotations in a transaction for speed
    conn.execute("BEGIN TRANSACTION", [])?;
    for ann in annotations {
        conn.execute(
            "INSERT INTO Contig_annotation (Contig_id, Start, End, Strand, Type, Product, Function, Phrog)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            rusqlite::params![
                ann.contig_id,
                ann.start,
                ann.end,
                ann.strand,
                ann.feature_type,
                ann.product,
                ann.function,
                ann.phrog,
            ],
        )?;
    }
    conn.execute("COMMIT", [])?;

    // Insert default variables (matching Python)
    let variables = [
        ("coverage", "Coverage", "Coverage", "curve", "#333333", 0.8, 0.4, 1.0, "Coverage depth", ""),
        ("coverage_reduced", "Coverage reduced", "Phage termini", "curve", "#00c53b", 0.8, 0.4, 1.0, "Coverage reduced", ""),
        ("reads_starts", "Reads termini", "Phage termini", "bars", "#215732", 0.6, 0.4, 1.0, "Read Starts", ""),
        ("reads_ends", "Reads termini", "Phage termini", "bars", "#6cc24a", 0.6, 0.4, 1.0, "Read Ends", ""),
        ("tau", "Tau", "Phage termini", "bars", "#44883e", 0.6, 0.4, 1.0, "Tau", ""),
        ("read_lengths", "Read lengths", "Assembly check", "curve", "#ed8b00", 0.8, 0.4, 1.0, "Read Lengths", ""),
        ("insert_sizes", "Insert sizes", "Assembly check", "curve", "#ed8b00", 0.8, 0.4, 1.0, "Insert Sizes", ""),
        ("bad_orientations", "Bad orientations", "Assembly check", "bars", "#c94009", 0.6, 0.4, 1.0, "Bad Orientations", ""),
        ("left_clippings", "Clippings", "Assembly check", "bars", "#7f0091", 0.6, 0.4, 1.0, "Left Clippings", ""),
        ("right_clippings", "Clippings", "Assembly check", "bars", "#8e43e7", 0.6, 0.4, 1.0, "Right Clippings", ""),
        ("insertions", "Indels", "Assembly check", "bars", "#e50001", 0.6, 0.4, 1.0, "Insertions", ""),
        ("deletions", "Indels", "Assembly check", "bars", "#97011a", 0.6, 0.4, 1.0, "Deletions", ""),
        ("mismatches", "Mismatches", "Assembly check", "bars", "#5a0f0b", 0.6, 0.4, 1.0, "Mismatches", ""),
    ];

    for (name, subplot, module, vtype, color, alpha, fill_alpha, size, title, help) in variables {
        let table_name = format!("Feature_{}", name);
        conn.execute(
            "INSERT INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            rusqlite::params![name, subplot, module, vtype, color, alpha, fill_alpha, size, title, help, table_name],
        )?;
    }

    Ok(())
}

/// Update Sample and Presences tables in SQLite database.
pub fn update_sample_presences(
    db_path: &Path,
    results: &[Result<(String, Vec<PresenceData>, f64)>],
    contigs: &[ContigInfo],
) -> Result<()> {
    let conn = Connection::open(db_path)?;

    // Build contig name -> id mapping
    let contig_name_to_id: HashMap<String, i64> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
        .collect();

    // Insert samples and presences
    conn.execute("BEGIN TRANSACTION", [])?;

    let mut sample_id = 1i64;
    for result in results {
        if let Ok((sample_name, presences, _)) = result {
            // Insert sample
            conn.execute(
                "INSERT INTO Sample (Sample_name) VALUES (?1)",
                [sample_name],
            )?;

            // Insert presences for this sample
            for presence in presences {
                if let Some(&contig_id) = contig_name_to_id.get(&presence.contig_name) {
                    conn.execute(
                        "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?1, ?2, ?3)",
                        rusqlite::params![contig_id, sample_id, presence.coverage_pct],
                    )?;
                }
            }
            sample_id += 1;
        }
    }

    conn.execute("COMMIT", [])?;
    Ok(())
}
