# Database Architecture Update: Feature_blob Migration

This document explains the database storage architecture after the migration from
per-feature tables (`Feature_*`) to a single compressed BLOB table (`Feature_blob`).

---

## Overview

Previously, each feature had its own DuckDB table (e.g., `Feature_mismatches`,
`Feature_left_clippings`, `Feature_primary_reads_plus_only`), each storing
run-length encoded (RLE) rows with `First_position`, `Last_position`, `Value`
columns. This created ~20 separate tables with ~1.7M rows per sample.

Now, all per-position feature data is stored in a single `Feature_blob` table.
Each row holds one feature for one contig/sample pair as a compressed binary BLOB.
The only exception is:

- **Contig-level tables** (`Contig_GCContent`, `Contig_GCSkew`, repeat tables) —
  these are contig-wide (no sample dimension) and remain as separate RLE tables.

---

## Database Schema

### Core Tables

| Table          | Description                                                                                                                                                                                                                                                       |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Contig`       | One row per contig. `Contig_id`, `Contig_name`, `Contig_length`, GC stats, etc.                                                                                                                                                                                   |
| `Sample`       | One row per BAM file. `Sample_id`, `Sample_name`, `Sequencing_type` (long/paired-short/single-short), read counts, circular mapping flag.                                                                                                                         |
| `Coverage`     | One row per (contig, sample). Mean/median/trimmed coverage, CV, roughness, aligned fraction. All integer-scaled (×100 or ×1,000,000). Primary key: `(Contig_id, Sample_id)`.                                                                                      |
| `Variable`     | Feature metadata for the plotting layer. Maps each feature to its display properties (color, plot type, subplot grouping) and its `Feature_table_name` (always `"Feature_blob"` for sample-level features, or a `Contig_*` table name for contig-level features). |
| `Feature_blob` | **The main data table.** One compressed BLOB per (contig, sample, feature). See below.                                                                                                                                                                            |

### Feature_blob Table

```sql
CREATE TABLE Feature_blob (
    Contig_id   INTEGER NOT NULL REFERENCES Contig,
    Sample_id   INTEGER NOT NULL REFERENCES Sample,
    Feature_id  SMALLINT NOT NULL,  -- 1-based index into VARIABLES array
    Data        BLOB NOT NULL,      -- compressed feature data
    PRIMARY KEY (Contig_id, Sample_id, Feature_id)
);
```

Each row stores all per-position data for one feature on one contig in one sample.
The `Feature_id` is a 1-based index into the compile-time `VARIABLES` array (defined
in `src/types.rs`), which is the single source of truth for all feature definitions.

### Contig-Level Tables

These have no `Sample_id` — they describe properties of the contig itself:

| Table                             | Content                                                                                                                    |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `Contig_GCContent`                | GC percentage per window. RLE: `(Contig_id, First_position, Last_position, Value)`                                         |
| `Contig_GCSkew`                   | GC skew per window (×100). Same RLE schema.                                                                                |
| `Contig_direct_repeat_count`      | Number of overlapping direct repeats at each position. Materialized from sweep-line algorithm over `Contig_directRepeats`. |
| `Contig_inverted_repeat_count`    | Same for inverted repeats.                                                                                                 |
| `Contig_direct_repeat_identity`   | Max percent identity of overlapping direct repeats.                                                                        |
| `Contig_inverted_repeat_identity` | Same for inverted repeats.                                                                                                 |
| `Contig_directRepeats`            | Raw BLAST repeat hits: `(Position1, Position2, Position1prime, Position2prime, Pident)`                                    |
| `Contig_invertedRepeats`          | Same for inverted repeats.                                                                                                 |

### Other Tables

| Table                    | Content                                                                                      |
| ------------------------ | -------------------------------------------------------------------------------------------- |
| `Contig_annotation_core` | Gene annotations (start, end, strand, type, product, locus_tag, etc.)                        |
| `Annotation_sequence`    | Nucleotide/protein sequences for annotations (separate to avoid bloating annotation queries) |
| `Phage_mechanisms`       | Detected phage packaging mechanisms per (contig, sample)                                     |
| `Phage_termini`          | Individual terminus peaks with statistics (SPC, tau, p-values)                               |
| `Misassembly`            | Per-contig/sample counts of assembly errors at ≥50% prevalence                               |
| `Microdiversity`         | Same at ≥10% prevalence                                                                      |
| `Side_misassembly`       | Clipping events at contig edges                                                              |
| `Topology`               | Circularisation metrics (circularising reads, insert size deviation)                         |
| `Codon_table`            | Lookup: 64 codons (id 0–63)                                                                  |
| `Contig_sequence`        | Full contig nucleotide sequence                                                              |
| `Column_scales`          | Scaling factors for all stored values                                                        |
| `Constants_for_plotting` | Metadata flags about database content                                                        |

---

## Feature List

There are 25 features defined in `VARIABLES` (in `src/types.rs`), grouped into modules.
The `Feature_id` in `Feature_blob` is the 1-based position in this array.

### Genome Module (contig-level, NOT in Feature_blob)

These are stored in dedicated `Contig_*` tables, not BLOBs, because they have no
sample dimension.

| #   | Feature                    | Subplot             | Plot Type |
| --- | -------------------------- | ------------------- | --------- |
| 1   | `direct_repeat_count`      | Repeat count        | Curve     |
| 2   | `inverted_repeat_count`    | Repeat count        | Curve     |
| 3   | `direct_repeat_identity`   | Max repeat identity | Curve     |
| 4   | `inverted_repeat_identity` | Max repeat identity | Curve     |
| 5   | `gc_content`               | GC content          | Curve     |
| 6   | `gc_skew`                  | GC skew             | Curve     |

### Coverage Module (dense BLOBs)

All values stored at full base-pair resolution. Every position has a value.

| #   | Feature                    | Subplot              | Scale | Description                              |
| --- | -------------------------- | -------------------- | ----- | ---------------------------------------- |
| 7   | `primary_reads`            | Primary alignments   | Raw   | Primary alignment depth at each position |
| 8   | `primary_reads_plus_only`  | Alignments by strand | Raw   | Forward-strand primary reads             |
| 9   | `primary_reads_minus_only` | Alignments by strand | Raw   | Reverse-strand primary reads             |
| 10  | `secondary_reads`          | Other alignments     | Raw   | Secondary alignment depth                |
| 11  | `supplementary_reads`      | Other alignments     | Raw   | Supplementary alignment depth            |
| 12  | `mapq`                     | MAPQ                 | ×100  | Mean mapping quality per position        |

### Misalignment Module (sparse BLOBs, filtered)

These are **event-based**: only positions with significant signal are stored.
A position is kept only if **both** conditions are met:

1. `event_count > coverage[position] × bar_ratio` (default bar_ratio = 10%)
2. `event_count > min_occurrences` (default = 2)

This filters out noise at high-coverage regions while retaining meaningful events
at low-coverage regions.

| #   | Feature           | Subplot    | Metadata          | Description                                                                                                                           |
| --- | ----------------- | ---------- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| 13  | `left_clippings`  | Clippings  | stats + sequence  | Soft/hard clips at read starts. Stats = mean/median/std of clip length. Sequence = dominant clipped bases.                            |
| 14  | `right_clippings` | Clippings  | stats + sequence  | Same at read ends.                                                                                                                    |
| 15  | `insertions`      | Indels     | stats + sequence  | Insertion events from CIGAR. Stats = mean/median/std of insertion length. Sequence = dominant inserted sequence.                      |
| 16  | `deletions`       | Indels     | (none)            | Deletion events from CIGAR. Count only, no metadata.                                                                                  |
| 17  | `mismatches`      | Mismatches | sequence + codons | Base mismatches from MD tag. Sequence = dominant alternate base. Codons = synonymous/non-synonymous classification + codon/AA change. |

### Long-Reads Module (dense BLOB)

| #   | Feature        | Subplot      | Scale      | Description                                      |
| --- | -------------- | ------------ | ---------- | ------------------------------------------------ |
| 18  | `read_lengths` | Read lengths | Dense, ×10 | Mean read length of reads covering each position |

### Paired-Reads Module (mixed encoding)

| #   | Feature                  | Subplot          | Encoding         | Description                                      |
| --- | ------------------------ | ---------------- | ---------------- | ------------------------------------------------ |
| 19  | `insert_sizes`           | Insert sizes     | Dense, ×10       | Mean insert size of pairs covering each position |
| 20  | `non_inward_pairs`       | Non-inward pairs | Sparse, filtered | Pairs on same contig but wrong orientation       |
| 21  | `mate_not_mapped`        | Missing mates    | Sparse, filtered | Reads whose mate is unmapped                     |
| 22  | `mate_on_another_contig` | Missing mates    | Sparse, filtered | Reads whose mate maps to a different contig      |

Features 20–22 use the same `bar_ratio + min_occurrences` filtering as the
Misalignment module.

### Phage Termini Module (mixed encoding)

| #   | Feature            | Subplot          | Encoding                     | Description                                                                    |
| --- | ------------------ | ---------------- | ---------------------------- | ------------------------------------------------------------------------------ |
| 23  | `coverage_reduced` | Coverage reduced | Dense, Raw                   | Coverage counting only "clean" reads (no clipping/mismatch at alignment edges) |
| 24  | `reads_starts`     | Reads termini    | Sparse, filtered, + sequence | Count of read 5' ends. Sequence = dominant clipped bases at that terminus.     |
| 25  | `reads_ends`       | Reads termini    | Sparse, filtered, + sequence | Count of read 3' ends.                                                         |

---

## BLOB Binary Format

Every BLOB starts with a 32-byte header:

```
Offset  Size  Field
------  ----  -----
0       4     Magic: b"TBB\x01"
4       1     Format version (1)
5       1     Metadata flags (bitfield: 0x01=sparse, 0x02=has_stats, 0x04=has_sequence, 0x08=has_codons)
6       1     Scale code (0=Raw, 1=×100, 2=×1000, 3=×10)
7       1     Number of zoom levels (2)
8       4     Contig length (u32 LE)
12      4     Base block offset (u32 LE)
16      4     Base block compressed size (u32 LE)
20      4     Zoom index offset (u32 LE)
24      4     Sparse metadata offset (u32 LE, 0 if none)
28      4     Sparse metadata compressed size (u32 LE, 0 if none)
```

### Dense Encoding (continuous features)

Used for coverage, mapq, read_lengths, insert_sizes, coverage_reduced.

**Compression pipeline:** raw i32 values → delta encoding → zigzag encoding →
variable-length integer (LEB128/varint) → Zstd compression (level 3).

Values are split into 65,536-position chunks for random-access decompression.

### Sparse Encoding (event-based features)

Used for mismatches, insertions, deletions, clippings, read termini, paired-read events.

**Stored data:**

- **Positions:** delta-encoded u32 → zigzag → varint → Zstd
- **Values:** zigzag-encoded i32 → varint → Zstd (no delta — values aren't correlated)
- **Event count:** u32

**Metadata block (optional, after positions+values):**

Stored in columnar layout for better compression:

1. **Stats** (if `has_stats` flag): arrays of `mean[n]`, `median[n]`, `std[n]` (i32, ×100)
2. **Sequences** (if `has_sequence` flag): `[length: u8, bytes...][n]` + `prevalence[n]` (i16, ×10 = percentage)
3. **Codons** (if `has_codons` flag): `category[n]`, `codon_id[n]`, `aa_id[n]` (u8 arrays)

### Zoom Levels

Every BLOB contains 3 pre-computed zoom levels (100bp, 1000bp, 10000bp bins)
for fast visualization at different scales. Contig_blob has 1 level (10000bp).

**Dense zoom bins**: mean (i32) per bin, one entry per bin. Zstd-compressed.

**Sparse zoom bins**: only nonzero bins are stored, matching base-resolution
sparse philosophy. Format per level (inside Zstd envelope):

- `nonzero_count` (u32)
- `bin_indices[n]`: delta + zigzag + varint encoded
- `max_values[n]`: zigzag + varint encoded

This avoids storing zeros for empty bins and keeps zoom output consistent
with base resolution (only positions/bins with events are returned).

---

## Value Scales

Values are stored as integers with a scale factor to preserve precision without
floating-point overhead:

| Scale Code | Factor | Used By                                                               |
| ---------- | ------ | --------------------------------------------------------------------- |
| 0 (Raw)    | 1      | primary_reads, secondary_reads, supplementary_reads, coverage_reduced |
| 1 (×100)   | 100    | mapq                                                                  |
| 2 (×1000)  | 1000   | All sparse/filtered features (relative to coverage)                   |
| 3 (×10)    | 10     | read_lengths, insert_sizes                                            |

When decoding: `float_value = stored_i32 / scale_factor`.

Sparse features use ×1000 scale because their values represent the fraction of
coverage at that position (e.g., 150 mismatches at 1000× coverage = 0.150 = stored as 150).

---

## Filtering: bar_ratio and min_occurrences

Sparse features (bars and sparse curves) do not store every position — only positions
where the signal is significant relative to local coverage. Two thresholds must both
be satisfied:

### bar_ratio (default: 10.0, meaning 10%)

A position is kept if:

```
event_count > coverage[position] × (bar_ratio / 100)
```

This makes the filter **adaptive to local coverage depth**. At 1000× coverage,
a mismatch needs >100 occurrences. At 10× coverage, it needs >1 occurrence.

### min_occurrences (default: 2)

A position is kept if:

```
event_count > min_occurrences
```

This prevents isolated single events (likely sequencing errors) from being stored.

### Combined filter

```
kept = (event_count > coverage × bar_ratio/100) AND (event_count > min_occurrences)
```

**Example:**

```
Position:      1      2      3      4      5
Coverage:    1000   1000     50     50   1000
Mismatches:     5      5      5      5      5
bar_ratio = 10.0, min_occurrences = 2

Position 1: 5 > 1000×0.10 = 100? NO  → filtered out
Position 2: 5 > 1000×0.10 = 100? NO  → filtered out
Position 3: 5 > 50×0.10 = 5?     NO  → filtered out
Position 4: 5 > 50×0.10 = 5?     NO  → filtered out
Position 5: 5 > 1000×0.10 = 100? NO  → filtered out

With bar_ratio = 5.0:
Position 3: 5 > 50×0.05 = 2.5? YES, 5 > 2? YES → KEPT
Position 4: 5 > 50×0.05 = 2.5? YES, 5 > 2? YES → KEPT
```

### Features affected by filtering

All sparse features are filtered: `left_clippings`, `right_clippings`, `insertions`,
`deletions`, `mismatches`, `non_inward_pairs`, `mate_not_mapped`,
`mate_on_another_contig`, `reads_starts`, `reads_ends`.

Dense features (coverage, mapq, read_lengths, insert_sizes, coverage_reduced) store
**every position** — no filtering.

---

## What Changed (Legacy Removal)

### Rust Side

- **Removed `FeaturePoint` struct** and `write_features()` function from `db.rs`.
  Previously, each feature position was written as a row to a dedicated `Feature_*`
  table. Now, only the BLOB path exists.

- **Removed `create_feature_table_if_needed()`** which dynamically created
  `Feature_*` tables with varying column schemas (some with `Mean`/`Median`/`Std`,
  some with `Sequence`/`Sequence_prevalence`, some with `Codon_*` columns).

- **Removed dual-write** from `processing.rs`. Previously, `add_features_from_arrays()`
  both populated a `Vec<FeaturePoint>` (for Feature_* tables) and encoded BLOBs.
  Now it only encodes BLOBs.

- **Removed `feature_table_name()`** from `types.rs` which mapped feature names to
  table names like `"Feature_mismatches"`, `"Feature_left_clippings"`, etc.

### Python Side

- **Removed legacy SQL fallback** in `plotting_data_per_sample.py`. When the BLOB
  path failed to find data, it previously fell through to queries against `Feature_*`
  tables. Now it skips the feature with `continue`.

- **Removed Feature_* query branches** for downsampling and full-resolution paths.
  Only `Feature_primary_reads` (a real table) and `Contig_*` tables are queried
  via SQL. All other features come from BLOB decoding.

- **Removed legacy constants**: `RELATIVE_SCALED_FEATURES`, `_FEATURES_WITH_STATS`,
  `_FEATURES_WITH_SEQUENCES` — these described which Feature_* tables had which
  columns. No longer needed since BLOB metadata flags encode this information.

- **Simplified CSV export** in `downloading_data.py`. The UNION ALL query builder
  now only handles `Contig_*` tables. All sample-level features are exported by
  decoding BLOBs in Python.

### Backward Compatibility

Old databases with `Feature_*` tables are **not supported** after this change.
The Python plotting layer will skip any feature whose BLOB is not found, rather
than falling back to legacy table queries.
