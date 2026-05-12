# Phage Packaging Classification Algorithm

The `processing_phage_packaging.rs` module classifies phage packaging mechanisms through a multi-step filtering and classification pipeline.

## 1. Computing Read Starts and Ends

The phage termini detection relies on **reads_starts** and **reads_ends** arrays that count where reads begin and end along the contig.

### Filtering: Only Matching Termini

A read's terminus is only counted if it starts with a **match** within the first min_clipping_length basepairs (default 5). This is verified using both the CIGAR string and MD tag. Reads starting with clippings or insertions longer than 5bp are excluded because they likely represent sequencing artifacts or misalignments rather than true biological termini. Reads where the terminus is counted are called "clean reads".

### Short vs Long Reads

**Short reads**: For paired-end sequencing, the 5' end of read1 (+ strand) represents the fragment start; the 5' end of read2 (- strand) represents the fragment end. Thus:

- reads_starts = start positions of + strand reads
- reads_ends = start positions of - strand reads (which are fragment 3' ends)

For single-end sequencing, only one fragment terminus is observed. The read is classified as `reads_starts` if it maps to the + strand and as `reads_ends` if it maps to the − strand.

**Long reads**: Each read is **split in half** and each terminus is evaluated independently. This prevents losing ~80% of reads that have clipping at one end but a clean match at the other. Both strands contribute to both arrays:

- reads_starts = 5' ends from both strands (start_plus + start_minus)
- reads_ends = 3' ends from both strands (end_plus + end_minus)

## 2. Contig Selection

Only contigs where ≥ **min_aligned_fraction** (default 90%) of positions have at least one aligned read are processed.

## 3. Terminal repeats (DTR/ITR) are identified:

Repeat detection results (BLAST self-alignment) are filtered to keep real terminal repeats:

- Identity ≥ **min_identity_dtr** (default 90%)
- One region within **max_distance_duplication** (default 100bp) of contig start
- Other region within **max_distance_duplication** of contig end

These regions enable:

- **Both-copies confirmation** for clipping validation (DTR only)
- Deduplication of peaks at equivalent positions in both DTR copies (at classification time)

## 3. Position Filtering

A position must pass two criteria to be considered as potential termini:

- read_starts ≥ **min_frequency** × coverage_reduced (default 10%)
- read_starts ≥ **min_events** (default 10)

Where the raw count used is read_starts (for start positions) or read_ends (for end positions).

## 4. Peak Merging

Positions passing the filter are merged into peak areas:

- Nearby positions within **max_distance_peaks** (default 20bp) merge together
- For circular genomes, wrap-around merging occurs
- Each area stores: start_pos, end_pos, peak_size, center_pos (position with highest signal), total_spc (spc stands for Starting Position Coverage)

## 5. Peak Testing

2 tests are performed to identify real termini. The Poisson test validates that the peak itself is significant (more starts than expected under uniform fragmentation). The subsequent clipping test checks whether a significant peak is an artifact (excess clippings indicating misalignment). These are complementary: the Poisson test filters noise, the clipping test filters artifacts. An area must pass **both** tests (`passed_poisson_test && passed_clipping_test`) to be kept.

First, each peak area is tested for statistical significance using a **true Poisson exact test**. Under uniform random fragmentation, read starts at each position follow a Poisson distribution with rate proportional to local coverage.

**Background rate computation:**

- `pstart = sum(reads_starts) / sum(coverage_reduced)` (computed separately for starts and ends)

**Per-window expected count:**

- The window spans from `start_pos` to `end_pos` of the merged area
- `λ_W = sum(coverage_reduced[i] × pstart)` for all positions i in the window
- This accounts for variable coverage across the window

**Observed count:**

- `X_W = area.total_spc` (total read starts/ends in the window)

**Significance test:**

- Model: `X_W ~ Poisson(λ_W)`
- p-value: `P(X >= X_W | Poisson(λ_W))` (upper-tail, via `statrs` crate)
- Bonferroni correction: `adjusted_pvalue = pvalue × K` where K = number of windows tested (separate K for starts and ends)
- Areas with `adjusted_pvalue > 0.05` pass the test

## 6. Clipping Testing

### 6a. Clipping Pre-filter

Individual clipping positions are pre-filtered before aggregation. A clipping position is significant if:

- clippings ≥ **min_frequency** × primary_reads (default 10%)
- clippings ≥ **min_events** (default 10)

For start areas: evaluate left_clippings at each position
For end areas: evaluate right_clippings at each position

### 6b. DTR Both-Copies Confirmation

For **DTR only** (not ITR): real biological clippings appear at both DTR copies at equivalent positions. Boundary artifacts (reads can't extend past contig edge) appear at only one copy.

Any clipping that lacks a counterpart at the equivalent position in the other DTR copy is zeroed out. This prevents false positives from contig boundary effects.

### 6c. Area Clipping Aggregation

For each peak area, sum all pre-filtered clippings:

- Within the area bounds (start_pos to end_pos)
- OR within **max_distance_peaks** of area edges

### 6d. Statistical Test

Test whether a peak area has **significantly more clippings** than expected globally using a **z-test with a normal approximation to a binomial model**.

**Expected clippings calculation:**

- `clipped_ratio = (primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean`
- `expected_clippings = clipped_ratio × total_spc`

**Decision logic:**

- If `sum_clippings > expected_clippings` significantly → **DISCARD** (`passed_clipping_test = false`)
- Otherwise → **KEEP**

Statistical significance is determined using **clipping_significance** (default 0.05 → z_critical = 1.645). Only peaks that do not have too many clippings within their area are kept.

## 7. DTR Deduplication

Before classification, remaining peak areas are deduplicated:

- Each area's `center_pos` is translated to the first DTR region (canonical position)
- Duplicate areas at the same canonical position (within 20bp tolerance) are removed
- Result: unique peak count for classification

**Important:** Canonical positions ensure correct distance calculation for classification (e.g., T7 phage: start at pos 1, end at pos 160 → DTR distance = 159bp, not ~39000bp if using second copy position).

## 8. Classification

Based on unique peak count after deduplication:

| Start Areas | End Areas | Classification                 |
| ----------- | --------- | ------------------------------ |
| 0           | 0         | No_packaging                   |
| 1           | 0         | PAC                            |
| 0           | 1         | PAC                            |
| 1           | 1         | See distance-based rules below |
| >1 or >1    |           | Unknown_packaging              |

For contigs harboring DTR, 2 peaks that are located at equivalent positions on both DTR count as only 1 peak.

### Distance-based Classification (1 start, 1 end)

For ITR configuration (both peaks in ITR regions, distance ≤20bp):

- ITR length ≤1000bp → **ITR_short_5'/3'**
- ITR length ≤10% genome → **ITR_long_5'/3'**
- ITR length >10% genome → **ITR_outlier_5'/3'**

Otherwise, based on distance between canonical positions:

- <2bp → **COS** (blunt cohesive ends)
- ≤20bp → **COS_5'/3'** (cohesive with overhang)
- ≤1000bp → **DTR_short_5'/3'**
- ≤10% genome → **DTR_long_5'/3'**
- \>10% genome → **DTR_outlier_5'/3'**

The suffix (_5' or _3') indicates overhang orientation based on whether end comes before or after start in genomic coordinates.

## 9. Output

Results are stored as two DuckDB views in the database: `Explicit_phage_mechanisms` (one row per contig×sample) and `Explicit_phage_termini` (one row per terminus area).

### Explicit_phage_mechanisms

| Column                         | Description                                                                                                                                      |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Contig_name                    | Contig identifier                                                                                                                                |
| Sample_name                    | Sample identifier                                                                                                                                |
| Packaging_mechanism            | Classification string (e.g., "DTR_short_5'")                                                                                                     |
| Left_termini                   | Comma-separated center positions of kept start areas                                                                                             |
| Median_left_termini_clippings  | Comma-separated median clipping lengths of clean starts per left terminus (one value per terminus, same order as Left_termini)                   |
| Right_termini                  | Comma-separated center positions of kept end areas                                                                                               |
| Median_right_termini_clippings | Comma-separated median clipping lengths of clean starts per right terminus (one value per terminus, same order as Right_termini)                 |
| Duplication                    | "DTR" if all peaks in DTR, "ITR" if all in ITR, NULL otherwise                                                                                   |
| Total_peaks                    | Count of areas that passed both Poisson and clipping tests                                                                                       |
| Repeat_length                  | Distance between start and end peak centers (if 1 of each)                                                                                       |
| Terminase_distance             | Distance between the two peak centers (if applicable, from terminase large subunit annotation — comes from external metadata, not this pipeline) |
| Terminase_percentage           | `Terminase_distance / Contig_length × 100`                                                                                                       |

### Explicit_phage_termini

Contains all columns from `Explicit_phage_mechanisms` (Contig_name through Terminase_percentage) plus the following per-area columns:

| Column              | Description                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Start               | First position of the merged area                                                                                                                                                                                                                                                                                                                                                                                         |
| End                 | Last position of the merged area                                                                                                                                                                                                                                                                                                                                                                                          |
| Size                | Distance between Start and End + 1                                                                                                                                                                                                                                                                                                                                                                                        |
| Center              | Position with highest SPC in the area                                                                                                                                                                                                                                                                                                                                                                                     |
| Status              | "start" (left terminus) or "end" (right terminus)                                                                                                                                                                                                                                                                                                                                                                         |
| SPC                 | Total starting position coverage in the area                                                                                                                                                                                                                                                                                                                                                                              |
| Median_clippings    | Median clipping length across all clean reads in the area. Clean reads are those with clipping < min_clipping_length (default 5): exact-match reads contribute 0, near-match reads contribute their actual clip length (1–4). Values range from 0 to min_clipping_length − 1                                                                                                                                              |
| Coverage            | Coverage at center position                                                                                                                                                                                                                                                                                                                                                                                               |
| Tau                 | SPC / Coverage × 100 (stored as integer)                                                                                                                                                                                                                                                                                                                                                                                  |
| NumberPeaks         | Number of positions merged into this area                                                                                                                                                                                                                                                                                                                                                                                 |
| Passed_PoissonTest  | "yes" or "no"                                                                                                                                                                                                                                                                                                                                                                                                             |
| Expected_SPC        | Expected SPC from Poisson model (λ_W, rounded)                                                                                                                                                                                                                                                                                                                                                                            |
| Pvalue              | Raw Poisson p-value (compact format)                                                                                                                                                                                                                                                                                                                                                                                      |
| Adjusted_pvalue     | Bonferroni-adjusted p-value (compact format)                                                                                                                                                                                                                                                                                                                                                                              |
| Passed_ClippingTest | "yes" or "no"                                                                                                                                                                                                                                                                                                                                                                                                             |
| Clippings           | Sum of pre-filtered clippings in/near area                                                                                                                                                                                                                                                                                                                                                                                |
| Clipping_excess     | Relative amount of clipped reads to the number of reads starting with a match: `100*(#reads-# reads_starting_with_a_match)/# reads_starting_with_a_match`. `#reads-# reads_starting_with_a_match` corresponds to the number of clipped reads. For long reads, each read is split in two to account for both termini so the formula becomes:  `100*(2*#reads-# reads_starting_with_a_match)/# reads_starting_with_a_match` |
| Expected_clippings  | `ROUND(expected_clippings)` — expected number of clippings under the null hypothesis, computed from the local SPC and the clipping_excess value for this contig/sample                                                                                                                                                                                                                                                    |

### Accessing the results

Both DuckDB views can be explored with a database viewer allowing DuckDB databases like DBeaver or directly in the terminal via:

- duckdb <DATABASE_NAME> "SELECT * FROM Explicit_phage_mechanisms"

- duckdb "SELECT * FROM Explicit_phage_termini"

The `Explicit_phage_mechanisms` results are also integrated in the visualization local webpage.
