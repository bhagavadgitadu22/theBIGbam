# Filtering logic

The Filtering section lets you build complex queries to narrow down which contig/sample pairs are displayed in your plots. 

Each filter row consists of: 

- A **category**: Contig, Sample, Presences, Completeness, or Termini

- A **metric** within that category: examples "Contig length" for Contig, "Coverage mean" for Presences, "Sample name" for Sample

- A **comparison operator**: "=", "!=", ">" and "<" for numeric columns, "="" and "!=" for text columns

- An **input value**: a numeric value for numeric columns, or a text value for other columns. For text inputs, a list of existing column values is provided as autocomplete suggestions

- A **"-" button** to remove the filter row

Multiple conditions can be added within a section (visually grouped by a green bar on the left) using the green "+ Add AND/OR" button. 

Multiple sections can be added using the pink "+ Add AND/OR" button. 

Clicking “+ Add AND/OR” adds a new **connector** and a new **filter row**. You can choose whether the logical operator is **AND** or **OR**.

<img title="" src="../static/FILTERING.png" alt="image" data-align="center" width="323">

This two-level grouping allows expressive compound queries. For example, the example above is equivalent to the condition:

$(Contig\_length>50kbp AND Sequencing\_type="long") OR (Aligned\_fraction>70% OR Circularising\_reads>10)$

Once your filters are set, the available contigs and samples update accordingly, and only matching pairs will appear in the generated plots. This is particularly useful for focusing your analysis on subsets of interest, such as high-coverage contigs, specific packaging mechanisms, or contigs meeting completeness thresholds.

# Available filters per category

## Contig

Those metrics are available only if a Genbank file was provided (or at least an assembly file for some of them).

| Metric                  | Definition                                                                                                                           |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Contig length           | Length of the contig sequence, in base pairs (bp)                                                                                    |
| GC (mean, median or sd) | Mean, median, or standard deviation of GC content (%) calculated across the contig using 100 bp sliding windows                      |
| Duplication (%)         | Proportion of the contig that is repeated at least once, calculated from an autoblast of the contig                                  |
| Annotation metrics      | Represents the info extracted from the Genbank file: for pharokka for instance, 3 columns Product, Function and Phrog can be queried |

Additional columns might be available depending on the information you added per contig using **the bigbam add-contig-metadata** command.

## Sample

| Metric                 | Definition                                                                           |
| ---------------------- | ------------------------------------------------------------------------------------ |
| Sequencing type        | Type of sequencing technology used for the sample (long, paired-short, single-short) |
| Number of reads        | Total number of sequencing reads generated for the sample                            |
| Number of mapped reads | Number of sequencing reads that successfully map to the contigs                      |

Additional columns might be available depending on the information you added per sample using **the bigbam add-sample-metadata** command.

## Presences

| Metric                                                   | Definition                                                                                                                                                                                                                                                                                                                                                                                      |
| -------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Aligned fraction (%)                                     | Percentage of the contig length covered by at least one mapped read                                                                                                                                                                                                                                                                                                                             |
| Coverage (mean, median, sd)                              | Mean, median, and standard deviation of number of reads mapping per position                                                                                                                                                                                                                                                                                                                    |
| Coverage variation                                       | Measures coverage smoothness: a low value means coverage changes gradually across positions, a high value means there are sharp jumps. A sum of coverage difference between consecutive positions is computed: $coverage\_variation = sqrt( Σ(cov[i+1] - cov[i])² / (n-1) ) / mean\_coverage$                                                                                                   |
| Coverage mean/median corrected by number of reads        | Mean or median coverage normalized to account for differences in sequencing depth between samples. The **normalization ratio by number of reads** corresponds to each sample's total read count divided by the minimum total read count across all samples. This produces a ratio relative to the least-sequenced sample (e.g., a sample with 10k reads when the min is 5k gets a ratio of 2.0) |
| Coverage mean/median corrected by number of mapped reads | Mean or median coverage normalized to account for differences in sequencing depth between samples (accounting for mapped-reads only). The same logic is used to compute the **normalization ratio by number of mapped reads** but using mapped reads instead of total reads                                                                                                                     |

Metrics other than the coverage mean and median are independent of sequencing-depth and do not need normalization.

## Completeness

### Inner completeness metrics

| Metric                             | Definition                                                                                                                                                                                                                                                               |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Mismatch frequency                 | Per-position mismatch counts (from the MD tag) divided by local coverage, summed across all positions                                                                                                                                                                    |
| Deletion frequency                 | Per-position deletion events normalized by local coverage, weighted by the median deletion length at that position                                                                                                                                                       |
| Insertion frequency                | Same logic as deletions but for insertions (bases in reads but not in reference)                                                                                                                                                                                         |
| Read-based clipping frequency      | Sum of soft/hard clipping events at each position, normalized by coverage and weighted by median clip length. Represents sequence present in reads but not aligned to the reference                                                                                      |
| Reference-based clipping frequency | When a right-clipping event is immediately followed by a left-clipping event (sorted by position), the gap between them represents sequence in the reference not covered by reads. The contribution is weighted by the average prevalence (count/coverage) of both clips |
| Whole completeness (%)             | Where:<br> $score\_completeness = mismatch\_frequency + insertion\_frequency + read\_based\_clipping\_frequency$<br>These represent bases missing from the reference but present in reads (i.e. the reference is incomplete)                                             |
| Whole contamination (%)            | Where: <br> $score\_contamination = mismatch\_frequency + deletion\_frequency + reference\_based\_clipping\_frequency$<br>These represent bases present in the reference but not in reads (i.e. the reference has extra sequence)                                        |

### Side completeness metrics

| Metric                             | Definition                                                                                                                                                                                                                                                                                                                                        |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Left/right completeness (%)        | **For left:** Identifies the leftmost significant\* soft-clipping event and checks whether the median clipped length extends beyond the start of the contig (position 1). If so:<br>$Left\_completeness = 100-(clip\_count / coverage × 100)$ <br>**For right:** same logic with the rightmost right-clipping position relative to the contig end |
| Distance contaminated left/right   | The distance from the terminal clipping event to the contig boundary (in bp)                                                                                                                                                                                                                                                                      |
| Min distance left/right            | The median clipped length at the terminal clipping position                                                                                                                                                                                                                                                                                       |
| Circularising reads (number and %) |                                                                                                                                                                                                                                                                                                                                                   |

*A significant clipping event depends on the `--coverage_percentage` threshold. By default, this is 10%, meaning that only positions where more than 10% of locally mapped reads are clipped at that position are considered significant.

Warning on long reads tends to give long clippings
