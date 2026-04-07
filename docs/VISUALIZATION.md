# Binning and Downsampling Strategy

theBIGbam uses pre-computed zoom levels embedded in BLOB storage to render plots efficiently for large genomic regions. This document explains the binning levels, how zoom selection works, and the hard limits for certain tracks.

## Overview

When displaying a large genomic window (e.g., viewing an entire 5 Mbp chromosome), rendering every data point would overwhelm both the browser and the user. Each feature BLOB contains per-position base data plus pre-computed zoom levels at multiple bin sizes. The viewer selects the appropriate level based on the current window size.

**Key principle**: Bars use the **maximum value** within each bin to preserve spikes; curves use the **mean** for smooth profiles.

---

## BLOB Zoom Levels

### Sample-Level Features

Features stored per sample (coverage, clippings, insertions, deletions, mismatches, splicings, read termini, etc.) have **3 zoom levels**:

| Level | Bin Size  | Use Case                                         |
| ----- | --------- | ------------------------------------------------ |
| 1     | 100 bp    | 10 kb - 1 Mb windows (100-10,000 points)         |
| 2     | 1,000 bp  | 1 Mb - 10 Mb windows (1000-10,000 points)        |
| 3     | 10,000 bp | > 10 Mb windows ((1000-genome_size/10000 points) |

Each zoom bin stores: `bin_start`, `bin_end`, `max`, `mean`.

### Contig-Level Features

Features independent of sample (GC content, GC skew, repeat count, repeat identity) have **1 zoom level**:

| Level | Bin Size  | Use Case           |
| ----- | --------- | ------------------ |
| 1     | 10,000 bp | Large windows only |

GC content (500 bp window) and GC skew (1000 bp window) are already pre-windowed, so base resolution is used up to 10 Mbp windows before falling back to the zoom level.

---

## Zoom Level Selection

When plotting, the viewer picks the resolution based on the current window size and the **base resolution threshold** (configurable via the "Feature plots without binning (bp)" spinner, default **10,000 bp**):

1. If `window <= threshold`: use **base resolution** (per-position data from BLOB)
2. If `window > threshold`: try zoom levels in order [100, 1000, 10000], picking the **smallest bin size** where `window / bin_size <= 10,000` data points
3. If no level fits: fallback to **10,000 bp** bins (coarsest zoom)

### Rendering by Plot Type

- **Bars** (clippings, insertions, deletions, mismatches, read termini): use `max` from zoom bins to preserve spikes
- **Curves** (coverage, GC content, insert sizes, read lengths, splicings, etc.): use `mean` from zoom bins for smooth profiles

---

## Hard Limits (Track Disabled Beyond Threshold)

Some tracks cannot be meaningfully binned and are hidden when the viewing window exceeds a configurable threshold.

### Gene Map

| Parameter       | Default                 | Configurable                  |
| --------------- | ----------------------- | ----------------------------- |
| Max window size | **100,000 bp** (100 kb) | Yes — "Gene map (bp)" spinner |

Gene annotations require rendering individual arrows, labels, and strand indicators. For large windows, the gene map is hidden.

### Nucleotide Sequence

| Parameter       | Default             | Configurable                        |
| --------------- | ------------------- | ----------------------------------- |
| Max window size | **1,000 bp** (1 kb) | Yes — "Sequence plots (bp)" spinner |

The nucleotide sequence track displays individual A/T/G/C letters. Beyond ~1,000 bp, characters become illegible.

### Translated Sequence (Amino Acids)

| Parameter       | Default             | Configurable                        |
| --------------- | ------------------- | ----------------------------------- |
| Max window size | **1,000 bp** (1 kb) | Same spinner as nucleotide sequence |

Same as nucleotide sequence — individual amino acid letters need to be readable.

---

## Configuration Summary

| Track Type          | Default Limit | When Exceeded                    |
| ------------------- | ------------- | -------------------------------- |
| Gene map            | 100,000 bp    | Track hidden                     |
| Nucleotide sequence | 1,000 bp      | Track hidden                     |
| Translated sequence | 1,000 bp      | Track hidden                     |
| Feature plots       | 10,000 bp     | Zoom levels (100/1k/10k bp bins) |

---

## User Interface Controls

Located in the **Plotting parameters** section:

### Max window size for plotting

- **Gene map (bp)**: Spinner to adjust gene map threshold (default: 100,000)
- **Sequence plots (bp)**: Spinner to adjust sequence tracks threshold (default: 1,000)
- **Feature plots without binning (bp)**: Spinner to adjust base resolution threshold (default: 10,000)

### Other plotting parameters

- **Minimum frequency for coverage-related features**: Filters out positions where the value is below this fraction of the maximum (default: 0.0 = no filtering). Applies only to sample-level coverage features, not contig-level features (GC content, GC skew, repeats).
