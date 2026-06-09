# Database Structure

theBIGbam stores all computed data in a single DuckDB file. This page describes the table layout and the BLOB encoding used for per-position data. This information is primarily relevant for users querying the database directly with SQL or developers interested in helping improve theBIGbam. 

---

## Table overview

Tables are grouped by what they describe. Not all tables are present in every database — tables marked *(BAM)* require BAM input, and tables marked *(MAG)* are only created in MAG mode. Annotation and sequence tables (e.g. Contig_annotation_core, Contig_sequence) are always created but remain empty if no annotation or assembly files are provided.

### Core entities

| Table                               | Key                   | Description                                                                                                                     |
| ----------------------------------- | --------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| **Contig**                          | `Contig_id`           | One row per contig. Stores name, length, GC statistics, duplication percentage. In contig mode, also stores `Number_of_samples` |
| **Sample** *(BAM)*                  | `Sample_id`           | One row per BAM file. Stores name, sequencing type, total/mapped read counts, circular mapping flag                             |
| **MAG** *(MAG)*                     | `MAG_id`              | One row per MAG bin. Stores name, total length, N50, number of contigs, GC statistics                                           |
| **MAG_contigs_association** *(MAG)* | `(MAG_id, Contig_id)` | Maps contigs to their MAG. `Offset_in_MAG` gives the cumulative bp offset under longest-first ordering                          |

### Per-contig-per-sample metrics

| Table                        | Key                      | Description                                                                                                                                      |
| ---------------------------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Coverage** *(BAM)*         | `(Contig_id, Sample_id)` | Contig-level coverage summaries: mean, median, trimmed mean, coefficient of variation, relative coverage roughness, aligned fraction, read count |
| **MAG_coverage** *(MAG)*     | `(MAG_id, Sample_id)`    | Same metrics as Coverage but aggregated at MAG level                                                                                             |
| **Misassembly** *(BAM)*      | `(Contig_id, Sample_id)` | Event counts at >=50% prevalence threshold: mismatches, deletions, insertions, clippings, collapse/expansion bp                                  |
| **Microdiversity** *(BAM)*   | `(Contig_id, Sample_id)` | Event counts at >=10% prevalence threshold, plus microdiverse bp on reference and reads                                                          |
| **Side_misassembly** *(BAM)* | `(Contig_id, Sample_id)` | Left/right contig-end clipping events: collapse/expansion at each end, misjoint mates                                                            |
| **Topology** *(BAM)*         | `(Contig_id, Sample_id)` | Circularisation metrics: circularising reads/inserts, median circularising length                                                                |

### Per-position data (BLOB storage)

| Table                          | Key                                             | Description                                                                     |
| ------------------------------ | ----------------------------------------------- | ------------------------------------------------------------------------------- |
| **Feature_blob** *(BAM)*       | `(Contig_id, Sample_id, Feature_id)`            | Zoom-level summary BLOB per mapping feature (coverage, MAPQ, mismatches, etc.)  |
| **Feature_blob_chunk** *(BAM)* | `(Contig_id, Sample_id, Feature_id, Chunk_idx)` | Base-resolution data in 65,536 bp chunks per mapping feature                    |
| **Contig_blob**                | `(Contig_id, Feature_id)`                       | Zoom-level summary BLOB per contig-level feature (GC content, GC skew, repeats) |
| **Contig_blob_chunk**          | `(Contig_id, Feature_id, Chunk_idx)`            | Base-resolution chunks for contig-level features                                |
| **MAG_blob** *(MAG)*           | `(MAG_id, Sample_id, Feature_id)`               | Zoom-level summary BLOB for MAG-scale mapping features                          |
| **MAG_contig_blob** *(MAG)*    | `(MAG_id, Feature_id)`                          | Zoom-level summary BLOB for MAG-scale contig features                           |

### Phage packaging

| Table                        | Key            | Description                                                                                        |
| ---------------------------- | -------------- | -------------------------------------------------------------------------------------------------- |
| **Phage_mechanisms** *(BAM)* | `Packaging_id` | One row per detected packaging mechanism per contig per sample                                     |
| **Phage_termini** *(BAM)*    | `Terminus_id`  | Individual terminus areas linked to `Phage_mechanisms`, with Poisson and clipping test diagnostics |

### Annotations

| Table                      | Key                              | Description                                                                                                                       |
| -------------------------- | -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| **Contig_annotation_core** | `Annotation_id`                  | Structural annotation data: contig, start, end, strand, type, main isoform flag, parent annotation link                           |
| **Annotation_segments**    | `(Annotation_id, Segment_index)` | Sub-intervals for spliced features (GenBank `join(...)`, GFF3 exon/CDS children). Only populated when a feature has >= 2 segments |
| **Annotation_qualifier**   | `(Annotation_id, Key)`           | Key-value store for all annotation qualifiers (product, locus_tag, gene, note, db_xref, ...)                                      |
| **Annotation_sequence**    | `Annotation_id`                  | Nucleotide and protein sequences, plus S-sites and N-sites counts for dN/dS                                                       |
| **Contig_qualifier**       | `(Contig_id, Key)`               | Key-value store for contig-level qualifiers from GenBank `source` features (organism, host, strain, ...)                          |
| **Contig_sequence**        | `Contig_id`                      | Full contig reference sequence                                                                                                    |
| **Codon_table**            | `Codon`                          | Standard genetic code lookup (64 codons)                                                                                          |
| **Annotated_types**        | `Type_id`                        | Distinct annotation type values ordered by frequency                                                                              |

### Repeats and BLAST hits

| Table                         | Key                                     | Description                                                             |
| ----------------------------- | --------------------------------------- | ----------------------------------------------------------------------- |
| **Contig_directRepeats**      | `(Contig_id, positions)`                | Direct repeats detected by nucmer. Pident stored as integer x100        |
| **Contig_invertedRepeats**    | `(Contig_id, positions)`                | Inverted repeats detected by nucmer. Pident stored as integer x100      |
| **Contig_blast_hits** *(MAG)* | `(Contig_id_1, Contig_id_2, positions)` | Inter-contig BLAST hits, stored canonically (id_1 <= id_2). Pident x100 |

### Visualization metadata

| Table                      | Key                           | Description                                                                                                       |
| -------------------------- | ----------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| **Variable**               | `Variable_id`                 | Feature display properties: subplot assignment, module, color, alpha, size, encoding type, help text              |
| **Column_scales**          | `(Feature_name, Column_name)` | Documents the integer scaling factor for each stored column so Python can decode back to real values              |
| **Constants_for_plotting** | `Constant`                    | Boolean metadata flags about database content (e.g., whether modules are present)                                 |
| **Color_templates**        | `Template_id`                 | Named color schemes (e.g., pharokka, generic)                                                                     |
| **Color_rules**            | `Rule_id`                     | Individual coloring rules belonging to a template, matching qualifier values to colors                            |
| **Database_metadata**      | `Key`                         | Key-value store for how the database was created/modified (theBIGbam version, creation date, modules, parameters) |

### Views

Views compute derived metrics on the fly from the tables above. They are not stored — removing a sample or contig automatically updates all view results.

| View                                        | Description                                                                                                             |
| ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| **Contig_annotation**                       | Joins `Contig_annotation_core` + `Annotation_sequence` + `Annotation_segments` into a single convenient row per feature |
| **Explicit_coverage**                       | Human-readable coverage with RPKM and TPM computed on the fly                                                           |
| **Explicit_misassembly**                    | Misassembly counts normalized per 100 kbp                                                                               |
| **Explicit_microdiversity**                 | Microdiversity counts normalized per 100 kbp                                                                            |
| **Explicit_side_misassembly**               | Side misassembly with normalized misjoint mates                                                                         |
| **Explicit_topology**                       | Circularisation metrics with human-readable names                                                                       |
| **Explicit_phage_mechanisms**               | Packaging mechanisms with contig/sample names                                                                           |
| **Explicit_phage_termini**                  | Terminus areas with contig/sample names and all diagnostics                                                             |
| **MAG_annotation_core** *(MAG)*             | Annotations shifted to MAG coordinates via `Offset_in_MAG`                                                              |
| **Contig_blast_hits_symmetric** *(MAG)*     | BLAST hits in both directions (both (A,B) and (B,A))                                                                    |
| **Explicit_coverage_per_MAG** *(MAG)*       | MAG-level coverage with RPKM and TPM                                                                                    |
| **Explicit_misassembly_per_MAG** *(MAG)*    | MAG-level misassembly counts per 100 kbp                                                                                |
| **Explicit_microdiversity_per_MAG** *(MAG)* | MAG-level microdiversity counts per 100 kbp                                                                             |

---

## Integer scaling conventions

All feature values are stored as **INTEGER** in DuckDB for efficient columnar compression. Real-valued quantities are scaled up before storage and scaled back down on read. The `Column_scales` table documents every scaling factor, but the main conventions are:

| Scale      | Used for                                                                | Example                      |
| ---------- | ----------------------------------------------------------------------- | ---------------------------- |
| x10        | Aligned fraction, duplication percentage                                | 95.3% stored as 953          |
| x100       | Coverage mean/median/trimmed mean, MAPQ, repeat Pident, GC sd/amplitude | 42.57x stored as 4257        |
| x1000      | Coverage-relative sparse features, GC content/skew per-position values  | 12.5 per-mille stored as 125 |
| x1,000,000 | Coefficient of variation, relative coverage roughness                   | 0.000832 stored as 832       |

---

## BLOB encoding

Per-position data is too large to store efficiently as one database row per genomic position. Instead, theBIGbam compresses each per-position signal into a binary BLOB. Each row stores a contig identifier, sample identifier, feature identifier, and a BLOB containing the values for all positions along the contig. This representation is **lossless**.

### Encoding pipeline

For dense features (one value per position), the pipeline is:

1. **Delta encoding:** store `[v[0], v[1]-v[0], v[2]-v[1], ...]` instead of raw values. Consecutive positions with similar values produce small deltas.
2. **Zigzag encoding:** map signed integers to unsigned integers. The sign is encoded in the least significant bit, ie positive numbers map to **even** values, negative numbers map to **odd** values. -1 becomes 1 (1 byte) instead of 0xFFFFFFFF (5 bytes)
3. **Varint encoding:** (LEB128) encode each unsigned integer as 1-5 bytes depending on magnitude. Small values (common after delta+zigzag) use 1 byte. 5 bytes is enough to include the full integer 32-bit range.
4. **Zstd compression:** (level 3) final byte-level compression on the varint stream. Exploits repeated byte sequences and statistical redundancies.

Sparse features differ from dense features in that values are present only at a limited number of genomic positions rather than at every position along the contig. As a result, both the positions and their associated values must be stored explicitly. The position array and value array are each compressed using the same delta+zigzag+varint+zstd pipeline and then packed together into a single BLOB. Sparse features may also carry per-event metadata, such as the dominant sequence, length statistics (mean, median, and standard deviation), or codon classification (e.g., synonymous, non-synonymous, or intergenic). This metadata is stored alongside the positions and values within the same BLOB.

### Decoding the BLOB data

When visualising a range of positions from a contig and clicking APPLY, BLOB data must be decoded to fetch and return the values associated to the positions requested. For huge contigs, to avoid decoding the entire length, dense feature BLOBs are split into **65,536 bp chunks** (`Feature_blob_chunk` table). Each chunk is independently compressed. When the user zooms into a genomic window, Python computes which chunk index(es) cover that window and fetches only those rows from Feature_blob_chunk, avoiding the cost of decoding the entire contig's data. Sparse BLOBs are small enough (only non-zero positions) to be stored and fetched whole — they are not chunked.

### Zoom levels

At low zoom levels, users may view very large genomic regions. To avoid the computational cost of decoding multiple chunks simultaneously and to prevent visual clutter in the plots, each feature has a pre-computed **zoom BLOB** (stored in the `Feature_blob` or `Contig_blob` table) containing summaries at three resolutions: 100 bp, 1 kbp, and 10 kbp. Only 1 zoom level (10 kbp) is used for GC content and GC skew as those features are already 500 and 1000 bp summaries. **Dense zoom bins** store the **mean** value across the bin. **Sparse zoom bins** store the **max** value across the bin (so peaks are never hidden when zoomed out).

With `thebigbam serve`, zoom BLOBs are used instead of per-position data when the displayed window exceeds 10 kbp. This threshold can be adjusted through the **Feature plots without binning (bp)** parameter in the **Plotting parameters** section. See [the main page](README.md#adaptive-resolution-rendering) for more details regarding the plotting with zooms.

Conversely, `thebigbam inspect` always returns base-resolution data, never zoom summaries.

### MAG BLOBs

In MAG mode, per-contig BLOBs are still written as usual. In addition, MAG-wide zoom-only BLOBs are built by concatenating per-contig data (shifted by each contig's Offset_in_MAG). MAG_blob stores one zoom BLOB per MAG per sample per mapping feature. MAG_contig_blob stores one zoom BLOB per MAG per contig-level feature (GC, repeats). 

No MAG-level base-resolution chunks are stored, in order to avoid redundant data duplication. For zoomed-in views, the viewer retrieves data directly from the per-contig chunk tables and uses coordinate offsets to map MAG-level positions back to contig coordinates.
