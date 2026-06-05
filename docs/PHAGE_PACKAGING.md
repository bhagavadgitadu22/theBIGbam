# Phage packaging classification

The classification is made through a multi-step filtering and classification pipeline. Steps 1-8 identify positions on the contig where significantly more reads than expected start or end their alignments. Those steps are applied independently to read start positions and read end positions. For simplicity, the descriptions below focus on read start positions, but the same procedure is applied identically to read end positions. In the final step (Step 9), information from both sets of positions is combined to assign the contig to a packaging mechanism category.

**Beware:** this works only for libraries prepared without fragmentation or with random fragmentation (i.e. restriction enzyme digestions would not produce enrichment).

## 1. Only keep reliable reads

To identify positions enriched in read starts and ends, we first filter primary alignments. We retain only reads that begin with a match within the first *min_clipping_length* base pairs (default: 5 bp). This condition is validated using both the CIGAR string and the MD tag. Reads starting with soft clipping or insertions longer than 5 bp are excluded, as they are more likely to reflect sequencing or alignment artifacts rather than true biological termini.

For **long reads**, we consider both ends of the alignment. Each read is therefore split into two segments, and the end of the second segment is treated as an additional start position. This allows both alignment boundaries to contribute to terminus detection.

The resulting set of filtered primary reads is used to compute a new coverage profile, referred to as *coverage_reduced*. Only those reads are considered to locate the contig termini.

## 2. Computing read starts and ends

The phage termini detection relies on *read_starts* and *read_ends* arrays that count where reads begin and end along the contig.

### Counting read termini depending on sequencing type

During **paired-end sequencing**, a DNA fragment is sequenced in 2 parts: one associated to the + strand (read 1) and its mate associated to the - strand (read 2). The 5' end of read 1 represents the fragment start (contribute to *read_starts* array) ; the 5' end of read 2 represents the fragment end (contribute to *read_ends* array).

For **single-end sequencing** only one fragment terminus is observed, the other terminus of a fragment is never sequenced. The 5' end of the read is counted in `reads_starts` if it maps to the + strand. Alternatively, the 3' end of the reads is counted in `reads_ends` if it maps to the − strand.

For **long reads**, we consider both ends of each read. For reads mapping on the + strand, the 5' end contribute to *reads_starts* and the 3' end to *reads_ends*. For reads mapping on the - strand, this is the opposite (3' end contribute to *reads_starts* and the 5' end to *reads_ends*).

## 3. Contig selection

Only contigs where ≥ *min_aligned_fraction* (default 90%) of positions have at least one aligned read are processed further.

## 4. Terminal repeats (DTR/ITR) are identified:

Repeat detection results (BLAST self-alignment) are filtered to keep real terminal repeats:

- Identity ≥ *min_identity_dtr* (default 90%)
- One of the repeats within *max_distance_duplication* (default 100 bp) of contig start
- The other repeat within *max_distance_duplication* of contig end

Identification of DTR regions enable the filtering of clippings and termini peaks later in the pipeline.

## 5. Position filtering

A position must pass two criteria to be considered as a potential start termini:

- *read_starts* ≥ *min_frequency* × coverage_reduced (default 10%)
- *read_starts* ≥ *min_events* (default 10)

## 6. Peak merging

Contig termini are not always confined to a single base position. They may be associated with several nearby positions due to biological variability in the DNA packaging and cleavage process. To account for this phenomenon, *read_starts* positions passing the filter are merged into peak areas:

- Nearby positions within *max_distance_peaks* (default 20 bp) are merged together. For genomes with circular mapping, position 1 and the last position on the contig are considered at a distance of 1 bp

Each peak area stores:

- beginning position

- finishing position

- center position (position where the peak with the highest signal is located)

- total SPC (Starting Position Coverage) with the total number of reads starting aligning in the peak area

## 7. Peak testing

2 tests are performed to identify real termini.

- **Poisson test:** validates that the peak area itself is significant. The peak should harbor more starts than expected under uniform fragmentation.

- **Clipping test**:** checks whether a significant peak is likely to be an artifact. High clipping rates near a peak indicate potential misalignment rather than a true terminus. This test requires mappings generated with `thebigbam mapping-per-sample --circular`, as standard linear mappings produce artificial clipping at contig ends, making genuine termini indistinguishable from alignment artifacts.

These are complementary: the Poisson test filters noise, the clipping test filters artifacts. An area must pass both tests to be kept.

### 7a. Details of the Poisson test

First, each peak area is tested for statistical significance using a **true Poisson exact test**. Under uniform random fragmentation, read starting at each position follow a Poisson distribution with rate proportional to local coverage.

**Per-area expected count:**

- Average expected number of starts per position: `pstart = sum(reads_starts) / sum(coverage_reduced)`
- `λ_W = sum(coverage_reduced[i] × pstart)` for all positions i in the peak
- This accounts for variable coverage across the window

**Observed count:**

- `X_W = area.total_spc` (total read starts in the window)

**Significance test:**

- Model: `X_W ~ Poisson(λ_W)`
- p-value: `P(X >= X_W | Poisson(λ_W))`
- Bonferroni correction: `adjusted_pvalue = pvalue × K` where K = number of windows tested
- Areas with `adjusted_pvalue <= 0.05` pass the test

### 7b. Details of the clipping test

#### Clipping pre-filter

Individual clipping positions are pre-filtered before aggregation. A clipping position is significant if:

- clippings ≥ *min_frequency* × primary_reads (default 10%)
- clippings ≥ *min_events* (default 10)

Left clippings are considered for start peak areas. Right clippings are considered for end peak areas.

#### DTR both-copies confirmation

Phages with direct terminal repeats (DTRs) are circular during part of their infection cycle and therefore exist in forms containing only a single copy of the repeat. As a result, sequencing reads (particularly long reads) may span the repeat region and contain sequence extending on both sides of the repeat.

If the assembled contig contains both copies of the DTR, these reads cannot align continuously across the junction between the two repeat copies. Even when using `thebigbam mapping-per-sample --circular`, they will therefore generate split alignments and an excess of clipped reads at the DTR boundaries.

These clipping events differ from clipping caused by misalignments. Misalignment-induced clipping typically occurs at both DTR copies, whereas clipping caused by reads spanning a single-repeat genome form occurs only at the artificial junction between the two DTR copies in the assembly.

To account for this effect, clipping signals are compared between the two DTR copies. Any clipping event that does not have a corresponding signal at the equivalent position in the other repeat copy is set to zero. This correction removes clipping artifacts caused by the assembly structure and prevents false positives during terminus detection.

#### Area clipping aggregation

For each peak area, sum all pre-filtered clippings:

- Located within the area bounds (start_pos to end_pos)
- OR within *max_distance_peaks* of area edges

#### Statistical test

Test whether a peak area has **significantly more clippings** than expected globally using a **z-test with a normal approximation to a binomial model**.

**Expected clippings calculation:**

- We calculate per-read probability of being clipped: `clipped_ratio = (primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean`
- For long reads, each read contributes to both termini so primary_reads is doubled: ``clipped_ratio = (2*primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean`
- If clippings were uniformly distributed, how many clippings do we expect in this area given the number of reads starting there: `expected_clippings = clipped_ratio × total_spc`

**Decision logic:**

Only peaks that do not have too many clippings within their area are kept:

- Areas where `sum_clippings` is significantly higher than `expected_clippings` are discarded
- Statistical significance is determined using **clipping_significance** (default 0.05 → z_critical = 1.645)

## 8. DTR deduplication

Before classification, remaining peak regions are deduplicated. When the same peak is detected in both copies of a direct terminal repeat (DTR), they are merged into the occurence from the first repeat copy.

This deduplication serves two purposes. First, it ensures that each biological terminus is counted only once, allowing the correct number of termini to be inferred. Second, it enables accurate distance calculations between termini during classification. For example, in the case of the T7 phage, the termini are located near positions 1 and 160, corresponding to a DTR length of 159 bp. Without deduplication, one of the termini could be assigned to the second DTR copy near the end of the contig, yielding an incorrect distance of approximately 39 kbp instead of the true DTR length.

## 9. Classification

Based on unique peak count after deduplication:

| Start Areas | End Areas | Classification                 |
| ----------- | --------- | ------------------------------ |
| 0           | 0         | No_packaging                   |
| 1           | 0         | PAC                            |
| 0           | 1         | PAC                            |
| 1           | 1         | See distance-based rules below |
| > 1         | OR > 1    | Unknown_packaging              |

### Distance-based Classification (1 start, 1 end)

For ITR configuration (both peaks in ITR regions, distance ≤ *max_distance_peaks*):

- ITR length ≤ 1000 bp → **ITR_short_5'/3'**
- ITR length ≤ 10% genome → **ITR_long_5'/3'**
- ITR length > 10% genome → **ITR_outlier_5'/3'**

Otherwise, based on distance between peaks:

- < 2 bp → **COS** (blunt cohesive ends)
- ≤ 20 bp → **COS_5'/3'** (cohesive with overhang)
- ≤ 1000 bp → **DTR_short_5'/3'**
- ≤ 10% genome → **DTR_long_5'/3'**
- \> 10% genome → **DTR_outlier_5'/3'**

The suffix (_5' or _3') indicates overhang orientation based on whether end comes before or after start in genomic coordinates.

## 10. Output

If you want to use the phage packaging results outside thebigbam visualization tool, the data structure is described below.

Results are exposed as 2 DuckDB views:

- `Explicit_phage_termini` contains one row per terminus area per contig per sample (results of step 8)

- `Explicit_phage_mechanisms` contains the final classification made per contig sample (results of step 9)

### Explicit_phage_termini

| Column              | Description                                                                                                                                                                                                                                                                                                        |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Start               | First position of the merged area                                                                                                                                                                                                                                                                                  |
| End                 | Last position of the merged area                                                                                                                                                                                                                                                                                   |
| Size                | Distance between Start and End + 1                                                                                                                                                                                                                                                                                 |
| Center              | Position with highest SPC in the area                                                                                                                                                                                                                                                                              |
| Status              | "start" (left terminus) or "end" (right terminus)                                                                                                                                                                                                                                                                  |
| SPC                 | Total starting position coverage in the area                                                                                                                                                                                                                                                                       |
| Median_clippings    | Median clipping length across all clean reads starting in the area. Clean reads are those with clipping < *min_clipping_length* (default 5): exact-match reads contribute 0, near-match reads contribute their actual clip length (1–4).                                                                           |
| Coverage            | Coverage at center position                                                                                                                                                                                                                                                                                        |
| Tau                 | SPC / Coverage × 100 (stored as integer)                                                                                                                                                                                                                                                                           |
| NumberPeaks         | Number of positions merged into this area                                                                                                                                                                                                                                                                          |
| Passed_PoissonTest  | "yes" or "no"                                                                                                                                                                                                                                                                                                      |
| Expected_SPC        | Expected SPC from Poisson model (λ_W, rounded)                                                                                                                                                                                                                                                                     |
| Pvalue              | Raw Poisson p-value (compact format)                                                                                                                                                                                                                                                                               |
| Adjusted_pvalue     | Bonferroni-adjusted p-value (compact format)                                                                                                                                                                                                                                                                       |
| Passed_ClippingTest | "yes" or "no"                                                                                                                                                                                                                                                                                                      |
| Clippings           | Sum of pre-filtered clippings in/near area                                                                                                                                                                                                                                                                         |
| Clipping_excess     | % of clipped reads relative to clean reads (those starting with a match). Formula: `100 × (primary_reads − clean_reads) / clean_reads`. For long reads, each read contributes to both termini so primary_reads is doubled: `100 × (2 × primary_reads − clean_reads) / clean_reads`. Stored as an integer (rounded) |
| Expected_clippings  | Expected number of clippings under the null hypothesis, computed from the local SPC and the clipping_excess value for this contig/sample                                                                                                                                                                           |

### Explicit_phage_mechanisms

| Column                         | Description                                                                                                                                                                                    |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Contig_name                    | Contig identifier                                                                                                                                                                              |
| Sample_name                    | Sample identifier                                                                                                                                                                              |
| Packaging_mechanism            | Classification string (e.g., "DTR_short_5'")                                                                                                                                                   |
| Left_termini                   | Comma-separated center positions of kept start areas                                                                                                                                           |
| Median_left_termini_clippings  | Comma-separated median clipping lengths of clean starts per left terminus (one value per terminus, same order as Left_termini)                                                                 |
| Right_termini                  | Comma-separated center positions of kept end areas                                                                                                                                             |
| Median_right_termini_clippings | Comma-separated median clipping lengths of clean starts per right terminus (one value per terminus, same order as Right_termini)                                                               |
| Duplication                    | "DTR" if all kept peaks are in DTR regions, "ITR" if all in ITR regions. NULL when status is mixed or no peaks were kept                                                                       |
| Total_peaks                    | Total number of terminus areas (both starts and ends) that passed both the Poisson test and the clipping test                                                                                  |
| Repeat_length                  | Genomic distance between the unique start and end peak centers. Only set when there is exactly 1 kept start area and 1 kept end area (after DTR deduplication), NULL otherwise                 |
| Terminase_distance             | Distance between the kept peak center(s) and the nearest terminase annotation (from GFF/GBK `product` qualifier matching "terminase"). NULL when no terminase annotation exists for the contig |
| Terminase_percentage           | `Terminase_distance / Contig_length × 100`                                                                                                                                                     |

### Accessing the results

Both DuckDB views can be explored with a database viewer allowing DuckDB databases like DBeaver or directly in the terminal via:

- duckdb <DB_NAME> "SELECT * FROM Explicit_phage_mechanisms"

- duckdb <DB_NAME> "SELECT * FROM Explicit_phage_termini"

The `Explicit_phage_mechanisms` results are also integrated in the visualization webpage.
