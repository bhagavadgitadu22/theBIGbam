# Binning and Downsampling Strategy

theBIGbam uses intelligent binning to render plots efficiently for large genomic regions. This document explains the binning strategies for different feature types and the hard limits for certain tracks.

## Overview

When displaying a large genomic window (e.g., viewing an entire 5 Mbp chromosome), rendering every data point would overwhelm both the browser and the user. theBIGbam applies SQL-side binning to reduce data volume while preserving biologically meaningful signals.

**Key principle**: Always preserve the **maximum value** within each bin to avoid hiding important spikes (e.g., clipping peaks, coverage drops, terminus signals).

---

## Hard Limits (No Binning — Track Disabled Beyond Threshold)

Some tracks cannot be meaningfully binned and are simply disabled when the viewing window exceeds a configurable threshold.

### Gene Map

| Parameter       | Default                 | Configurable                  |
| --------------- | ----------------------- | ----------------------------- |
| Max window size | **100,000 bp** (100 kb) | Yes — "Gene map (bp)" spinner |

**Rationale**: Gene annotations require rendering individual arrows, labels, and strand indicators. Binning would make the track unreadable. For large windows, hide the gene map entirely and let users zoom in.

**Message shown**: `Gene map not plotted: window > {threshold} bp`

### Nucleotide Sequence

| Parameter       | Default             | Configurable                        |
| --------------- | ------------------- | ----------------------------------- |
| Max window size | **1,000 bp** (1 kb) | Yes — "Sequence plots (bp)" spinner |

**Rationale**: The nucleotide sequence track displays individual A/T/G/C letters. Beyond ~1,000 bp, characters become illegible. Binning letters makes no biological sense.

**Message shown**: `Sequence not plotted: window > {threshold} bp`

### Translated Sequence (Amino Acids)

| Parameter       | Default             | Configurable                        |
| --------------- | ------------------- | ----------------------------------- |
| Max window size | **1,000 bp** (1 kb) | Same spinner as nucleotide sequence |

**Rationale**: Same as nucleotide sequence — individual amino acid letters need to be readable.

**Message shown**: `Translated sequence not plotted: window > {threshold} bp`

---

## Binning Strategy for Feature Plots

For curve and bar plots, binning is applied when the viewing window exceeds the **downsample threshold**.

### Downsample Threshold

| Parameter         | Value                           |
| ----------------- | ------------------------------- |
| Default threshold | **100,000 bp** (100 kb)         |
| Number of bins    | **1,000**                       |
| Bin width         | `window_size / 1000` (adaptive) |

When `(xend - xstart) > threshold`, SQL-side binning is triggered. Otherwise, full-resolution data is returned.

---

## Curve Features (Continuous Lines)

**Examples**: Coverage, GC content, GC skew, insert sizes, read lengths, repeat count, repeat identity, MAPQ, non-inward pairs, mate unmapped, mate on another contig

### Binning Method: MAX with Position Preservation

```
┌─────────────────────────────────────────────────────────────┐
│  Bin 1    │  Bin 2    │  Bin 3    │  Bin 4    │  Bin 5    │
│  MAX=45   │  MAX=67   │  MAX=23   │  MAX=89   │  MAX=56   │
│  pos=150  │  pos=380  │  pos=520  │  pos=712  │  pos=950  │
└─────────────────────────────────────────────────────────────┘
```

**Algorithm**:

1. Each RLE (run-length encoded) database row is expanded to all bins it overlaps
2. For each bin, find the **maximum value** (`MAX(Value)`)
3. Track the **mid-position** of the max-value run, clamped to bin boundaries
4. **Zero-fill** gaps between occupied bins to prevent interpolation artifacts in `varea` rendering

**Why MAX instead of MEAN?**

- Preserves spikes (e.g., coverage peaks at terminus areas)
- Highlights anomalies (e.g., sudden coverage drops indicating misassembly)
- Biologically meaningful: a single outlier position may be the signal of interest

**SQL pattern** (simplified):

```sql
SELECT bin_idx,
       MAX(Value) AS max_value,
       ARG_MAX(mid_pos, Value) AS position
FROM (
    SELECT Value, (First_position + Last_position) / 2 AS mid_pos,
           UNNEST(generate_series(bin_start_idx, bin_end_idx)) AS bin_idx
    FROM Feature_table
    WHERE ...
)
GROUP BY bin_idx
```

---

## Bar Features (Discrete Spikes)

**Examples**: Clippings (left/right), insertions, deletions, mismatches, read starts, read ends

### Binning Method: MAX with Full Position Preservation

Bar features represent discrete events at specific positions. The binning strategy differs from curves:

```
┌─────────────────────────────────────────────────────────────┐
│  Bin 1         │  Bin 2         │  Bin 3         │  Bin 4  │
│  ▮ pos=125     │  ▮ pos=380     │                │  ▮▮pos=712-715
│  MAX=45        │  MAX=67        │  (empty)       │  MAX=89 │
└─────────────────────────────────────────────────────────────┘
```

**Algorithm**:

1. For each bin, find the row with the **maximum value**
2. Preserve the **exact First_position and Last_position** of that row (not clamped)
3. **Deduplicate**: If the same run is the max in multiple adjacent bins, emit it only once
4. **No zero-fill**: Empty bins are not rendered (bars should be sparse)

**Why preserve exact positions?**

- Bar features mark specific genomic positions (e.g., "clipping peak at position 4521")
- Users need to identify the exact coordinate for further investigation
- Metadata (sequence, prevalence) is tied to the specific position

**Output includes**:

- `x_coords`: Midpoint of the bar
- `y_coords`: Value (count or relative frequency)
- `width_coords`: Width of the bar (`Last_position - First_position + 1`)
- `first_pos_coords`, `last_pos_coords`: Original coordinates for tooltip display
- `seq_coords`, `prev_coords`: (if applicable) Dominant sequence and prevalence

---

## Contig-Level vs Sample-Level Features

### Contig-Level Features (Sample-Independent)

These features have the same value regardless of which sample is selected:

| Feature                  | Table/View                               | Notes                    |
| ------------------------ | ---------------------------------------- | ------------------------ |
| GC content               | `Contig_GCContent`                       | Stored as RLE            |
| GC skew                  | `Contig_GCSkew`                          | Stored as INTEGER × 100  |
| Direct repeat count      | `Contig_direct_repeat_count` (view)      | Sweep-line aggregation   |
| Inverted repeat count    | `Contig_inverted_repeat_count` (view)    | Sweep-line aggregation   |
| Direct repeat identity   | `Contig_direct_repeat_identity` (view)   | MAX identity per segment |
| Inverted repeat identity | `Contig_inverted_repeat_identity` (view) | MAX identity per segment |

**Binning**: Uses the same `_rle_weighted_bin_sql()` function as sample-level features. The query omits `Sample_id` but otherwise follows identical MAX binning logic.

### Sample-Level Features

All other features (coverage, clippings, termini, etc.) are stored per-sample in `Feature_*` tables and require a `Sample_id` filter.

---

## Repeat Feature Views

Repeat features are special because the raw data is stored as interval pairs, not RLE values:

```
Contig_directRepeats:
  Position1, Position2, Position1prime, Position2prime, Pident
```

SQL views perform **sweep-line aggregation** to convert intervals to RLE-like segments:

1. Collect all interval boundaries (start and end+1 of each repeat region)
2. For each segment between consecutive boundaries:
   - **Count view**: COUNT of overlapping intervals
   - **Identity view**: MAX(Pident) ÷ 100 of overlapping intervals
3. Output format matches standard Feature tables: `(Contig_id, First_position, Last_position, Value)`

This allows repeat features to use the standard binning pipeline, fixing slow rendering on large contigs.

---

## Configuration Summary

| Track Type          | Default Limit | When Exceeded           |
| ------------------- | ------------- | ----------------------- |
| Gene map            | 100,000 bp    | Track hidden            |
| Nucleotide sequence | 1,000 bp      | Track hidden            |
| Translated sequence | 1,000 bp      | Track hidden            |
| Curve features      | 100,000 bp    | MAX binning (1000 bins) |
| Bar features        | 100,000 bp    | MAX binning with dedup  |

---

## User Interface Controls

Located in the **Plotting parameters → Max window size for plotting** section:

- **Gene map (bp)**: Spinner to adjust gene map threshold (default: 100,000)
- **Sequence plots (bp)**: Spinner to adjust sequence tracks threshold (default: 1,000)

The downsampling threshold for curve/bar features is currently not exposed in the UI (fixed at 100 kb).

Coverage normalization

For features that depend on local coverage depth (clippings, indels, mismatches, reads starts/ends, paired-read anomalies), an optional **"Plot relative to local coverage"** checkbox normalizes values to facilitate comparison across regions with varying coverage:

- **Enabled**: Y-axis shows the **ratio** of events to local coverage (value ÷ coverage), scaled between 0 and 1. For example, a value of 0.10 means 10% of reads at that position have the feature (e.g., 10 clippings per 100× coverage, or 50 clippings per 500× coverage).

- **Disabled** (default): Y-axis shows **absolute counts** (e.g., number of clippings, insertions, mismatches).

This normalization reveals whether anomalies are proportional to coverage (expected sequencing noise) or represent true biological signal or assembly errors that persist regardless of depth.
