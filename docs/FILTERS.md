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

*(Contig_length > 50kbp AND Sequencing_type = "long") OR (Aligned_fraction > 70% OR Circularising_reads > 10)*

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
| Coverage variation                                       | Measures coverage smoothness: a low value means coverage changes gradually across positions, a high value means there are sharp jumps. A sum of coverage difference between consecutive positions is computed:<br>*Coverage_variation = sqrt( Σ(Coverage[i+1] - Coverage[i])² / (n-1) ) / Mean_coverage*                                                                                        |
| Coverage mean/median corrected by number of reads        | Mean or median coverage normalized to account for differences in sequencing depth between samples. The **normalization ratio by number of reads** corresponds to each sample's total read count divided by the minimum total read count across all samples. This produces a ratio relative to the least-sequenced sample (e.g., a sample with 10k reads when the min is 5k gets a ratio of 2.0) |
| Coverage mean/median corrected by number of mapped reads | Mean or median coverage normalized to account for differences in sequencing depth between samples (accounting for mapped-reads only). The same logic is used to compute the **normalization ratio by number of mapped reads** but using mapped reads instead of total reads                                                                                                                     |

Metrics other than the coverage mean and median are independent of sequencing-depth and do not need normalization.

## Completeness

Those metrics aim to answer whether a contig is complete or contaminated, and if complete whether it is a linear contig or a circular one.

### Side completeness metrics

Left and right completeness are evaluated based on the presence of significant clipping events at the contig ends. A significant clipping event is determined by the `--coverage_percentage` threshold. By default, this is 10%, meaning that only positions where more than 10% of locally mapped reads are clipped are considered significant.

To assess completeness at the left end of a contig, significant left clipping events **where the median clipped length extends beyond the start of the contig (position 1)** are identified and the prevalence of those event are calculated as:

*Prevalence left clippings = Number of left clippings​ / Local coverage*

The position with the maximal prevalence (called the **terminal left-clipping event**) is used to defined the left completeness metrics.

| Metric                          | Definition                                                                                              |
| ------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Left/right completeness (%)     | *Left completeness = 100 − Prevalence(terminal left-clipping event)*. Same logic for right completeness |
| Left/right contamination length | The distance from the terminal left-clipping event to the contig start (in bp)                          |
| Left/right missing length       | The median clipped length at the terminal left-clipping event                                           |

The same logic applies for the right completeness metrics using right clipping events last position of the position.

**3 important warnings regarding those metrics:** 

- Side completeness metrics tend to **underestimate the incompleteness of the contig**, especially when using short reads. Short reads cannot extend far beyond the last properly assembled base, so Left/Right completeness and Left/Right contamination length are typically slightly lower than the true values. The **Left/Right missing length**, in particular, can be substantially smaller than the actual missing sequence. 

- For circular contigs, completeness values can sometimes **overestimate missing sequence**. This typically occurs when mappings are performed without theBIGbam or without the`--circular` option, **because reads that would normally span the junction between the last and first contig positions are clipped: standard mappers do not support circular alignments**. 

- In **rare cases**, this effect can also occur for long reads even when the `--circular` option is used. For example, in bacteriophages, single long reads may contain multiple full viral genomes due to concatemers. In such cases, only the portion of the read corresponding to one contig length (or up to twice the contig length when mapping with theBIGbam and the `--circular` option) is mapped, and the remaining sequence is clipped. If many concatemer-derived long reads are present, incompleteness metrics may detect this clipped signal, leading to artificially high estimates of missing sequence.

### Contig topology metrics

The circularity of complete contigs (e.g. full bacterial chromosomes or viral genomes) can only be reliably assessed when the `--circular` option is used.

| Metric              | Definition                                                                                                                                     |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| Circularising reads | Number of reads whose alignment starts near the end of the contig and continues past the contig start                                          |
| Circularity (%)     | Percentage of circularity, calculated as the number of circularising reads divided by the mean coverage at the first and last contig positions |

An approximate estimate of circularity could be derived from mappings performed without the `--circular` option, using supplementary alignments spanning both contig ends and, for paired-end short reads, non-inward read pairs with reasonable insert sizes mapping to both contig ends. However, this approach does not provide an exact estimation of circularity compared to circular mapping with the `--circular` option. For this reason, circularity is not estimated in this case.

### Whole completeness metrics

These metrics assess the overall completeness and contamination of the contig, considering both the contig ends and the internal regions. Only events deemed significant according to the `--coverage_percentage` threshold are taken into account.

| Metric                                                  | Definition                                                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Mismatch frequency (MF)                                 | Per-position mismatch counts (from the MD tag) divided by local coverage                                                                                                                                                                                                                                               |
| Deletion frequency (DF)                                 | Per-position deletion events normalized by local coverage and weighted by the median deletion length at each position (bases that are in reference but not in reads)                                                                                                                                                   |
| Insertion frequency (IF)                                | Per-position insertion events normalized by local coverage and weighted by the median insertion length at each position (bases that are in reads but not in reference)                                                                                                                                                 |
| Read-based clipping frequency (read-based CF)           | Per-position soft/hard clipping events normalized by coverage and weighted by median clip length. Represents sequence present in reads but not aligned to the reference                                                                                                                                                |
| **Completeness (%)**                                    | Where:<br>*Completeness = 100 x (1 - MF + IF + read-based CF)* <br>These represent bases missing from the reference but present in reads (i.e. the reference is incomplete)                                                                                                                                            |
| Reference-based clipping frequency (reference-based CF) | When a right-clipping event is immediately followed by a left-clipping event (sorted by position), the gap between them represents sequence in the reference not covered by reads. The contribution is weighted by the minimum prevalence (count/coverage) of both clips, so the weaker boundary determines confidence |
| **Contamination (%)**                                   | Where: <br>*Contamination = 100 x (1- MF + DF + Reference_based CF)* <br>These represent bases present in the reference but not in reads (i.e. the reference has extra sequence)                                                                                                                                       |

**Warnings:**

- Clipping frequencies are generally less precise than other frequency-based metrics because they capture only the onset of an issue, rather than its full extent, as is the case for mismatches, deletions, or insertions.

- For **Read-based clipping frequency**, the median clipped length provides an estimate of the number of bases missing from the reference. However, this value typically underestimates the true missing sequence, especially for short reads (see the warnings associated with the **Side completeness metrics** for further details).

- For **Reference-based clipping frequency**, the distance between consecutive left and right clipping events can be used to estimate contamination length. This estimate is indicative only and should be interpreted with caution.

## Termini

TODO
