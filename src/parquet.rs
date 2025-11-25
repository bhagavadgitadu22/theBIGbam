//! Parquet file writing functions.
//!
//! This module handles writing feature and presence data to Parquet files.

use anyhow::Result;
use arrow::array::{ArrayRef, Float32Array, Int32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use crate::types::{FeaturePoint, PresenceData};

/// Write features to Parquet file.
pub fn write_features_parquet(features: &[FeaturePoint], path: &Path) -> Result<()> {
    let schema = Schema::new(vec![
        Field::new("contig_name", DataType::Utf8, false),
        Field::new("feature", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("value", DataType::Float32, false),
    ]);

    let contig_names: Vec<&str> = features.iter().map(|f| f.contig_name.as_str()).collect();
    let feature_names: Vec<&str> = features.iter().map(|f| f.feature.as_str()).collect();
    let positions: Vec<i32> = features.iter().map(|f| f.position).collect();
    let values: Vec<f32> = features.iter().map(|f| f.value).collect();

    let batch = RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(StringArray::from(contig_names)) as ArrayRef,
            Arc::new(StringArray::from(feature_names)) as ArrayRef,
            Arc::new(Int32Array::from(positions)) as ArrayRef,
            Arc::new(Float32Array::from(values)) as ArrayRef,
        ],
    )?;

    let file = File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

/// Write presences to Parquet file.
pub fn write_presences_parquet(presences: &[PresenceData], path: &Path) -> Result<()> {
    let schema = Schema::new(vec![
        Field::new("contig_name", DataType::Utf8, false),
        Field::new("coverage_pct", DataType::Float32, false),
    ]);

    let contig_names: Vec<&str> = presences.iter().map(|p| p.contig_name.as_str()).collect();
    let coverage_pcts: Vec<f32> = presences.iter().map(|p| p.coverage_pct).collect();

    let batch = RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(StringArray::from(contig_names)) as ArrayRef,
            Arc::new(Float32Array::from(coverage_pcts)) as ArrayRef,
        ],
    )?;

    let file = File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}
