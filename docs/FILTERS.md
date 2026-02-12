# Filtering logic

The Filtering section lets you build complex queries to narrow down which contig/sample pairs are displayed in your plots.

Each filter row consists of:

- A **category**: Contig, Sample, Coverage, Misassembly, Microdiversity, Side misassembly, Topology, or Termini

- A **metric** within that category: examples "Contig length" for Contig, "Coverage mean" for Coverage, "Sequencing type" for Sample

- A **comparison operator**: "=", "!=", ">" and "<" for numeric columns, "=" and "!=" for text columns

- An **input value**: a numeric value for numeric columns, or a text value for other columns. For text inputs, a list of existing column values is provided as autocomplete suggestions

- A **"-" button** to remove the filter row

Multiple conditions can be added within a section (visually grouped by a green bar on the left) using the green "+ Add AND/OR" button.

Multiple sections can be added using the pink "+ Add AND/OR" button.

Clicking "+ Add AND/OR" adds a new **connector** and a new **filter row**. You can choose whether the logical operator is **AND** or **OR**.

<img title="" src="../static/FILTERING.png" alt="image" data-align="center" width="323">

This two-level grouping allows expressive compound queries. For example, the example above is equivalent to the condition:

*(Contig_length > 50kbp AND Sequencing_type = "long") OR (Aligned_fraction > 70% OR Circularising_reads > 10)*

Once your filters are set, the available contigs and samples update accordingly, and only matching pairs will appear in the generated plots. This is particularly useful for focusing your analysis on subsets of interest, such as high-coverage contigs, specific packaging mechanisms, or contigs meeting completeness thresholds.

# Available filters per category

## Contig

Those metrics are available only if a Genbank file was provided (or at least an assembly file for some of them).

| Metric | Definition |
| --- | --- |
| Contig_length | Length of the contig sequence, in base pairs (bp) |
| Duplication_percentage | Proportion of the contig covered by repeats (stored as integer x10; e.g. 125 = 12.5%), calculated from an autoblast of the contig |
| GC_mean | Mean GC content (%) calculated across the contig using non-overlapping 500 bp sliding windows |
| GC_sd | Standard deviation of GC content (stored as integer x100; e.g. 350 = 3.50) |
| GC_skew_amplitude | Amplitude of GC skew, i.e. *max(GC skew) - min(GC skew)* (stored as integer x100), calculated across non-overlapping 1kbp sliding windows |
| Positive_GC_skew_windows_percentage | Percentage of the 1kbp windows with positive skew (stored as integer x10; e.g. 523 = 52.3%) |
| Annotation columns | Columns extracted from the Genbank file: for pharokka for instance, Type, Product, Function and Phrog can be queried |

Additional columns might be available depending on the information you added per contig using **thebigbam add-contig-metadata** command.

## Sample

| Metric | Definition |
| --- | --- |
| Sequencing_type | Type of sequencing technology used for the sample (long, paired-short, single-short) |
| Number_of_reads | Total number of sequencing reads generated for the sample |
| Number_of_mapped_reads | Number of sequencing reads that successfully map to the contigs |

Additional columns might be available depending on the information you added per sample using **thebigbam add-sample-metadata** command.

## Coverage

| Metric | Definition |
| --- | --- |
| Aligned_fraction_percentage | Percentage of the contig length covered by at least one mapped read |
| Above_expected_aligned_fraction | Whether the aligned fraction is above what would be expected given the coverage depth |
| Read_count | Number of reads mapping to this contig for this sample |
| Coverage_mean | Mean number of reads mapping per position |
| Coverage_median | Median number of reads mapping per position |
| Coverage_trimmed_mean | Trimmed mean coverage (excluding extreme values) |
| RPKM | Reads Per Kilobase per Million mapped reads |
| TPM | Transcripts Per Million |
| Coverage_sd | Standard deviation of coverage across positions |
| Coverage_variation | Measures coverage smoothness: *Coverage_variation = sqrt( sum((Coverage[i+1] - Coverage[i])^2) / (n-1) ) / Mean_coverage*. A low value means gradual changes, a high value means sharp jumps |

## Misassembly

Positions with **>=50% prevalence** are counted. Prevalence = event count / local coverage.

| Metric | Definition |
| --- | --- |
| Mismatches_per_100kbp | Number of positions with mismatch prevalence >= 50%, normalized per 100kbp |
| Deletions_per_100kbp | Number of positions with deletion prevalence >= 50%, normalized per 100kbp |
| Insertions_per_100kbp | Number of positions with insertion prevalence >= 50%, normalized per 100kbp |
| Clippings_per_100kbp | Number of positions with clipping prevalence >= 50%, normalized per 100kbp |
| Collapse_bp | Bases present in reads but missing from reference (insertions + clippings + mismatches at >= 50% prevalence, weighted by median event length) |
| Collapse_per_100kbp | Collapse_bp normalized per 100kbp |
| Expansion_bp | Bases present in reference but missing from reads (deletions + mismatches + paired clip distances at >= 50% prevalence, weighted by median event length) |
| Expansion_per_100kbp | Expansion_bp normalized per 100kbp |

## Microdiversity

Same logic as Misassembly but with a **>=10% prevalence** threshold, capturing lower-frequency variants.

| Metric | Definition |
| --- | --- |
| Mismatches_per_100kbp | Number of positions with mismatch prevalence >= 10%, normalized per 100kbp |
| Deletions_per_100kbp | Number of positions with deletion prevalence >= 10%, normalized per 100kbp |
| Insertions_per_100kbp | Number of positions with insertion prevalence >= 10%, normalized per 100kbp |
| Clippings_per_100kbp | Number of positions with clipping prevalence >= 10%, normalized per 100kbp |
| Microdiverse_bp_on_reads | Bases present in reads but not in reference at >= 10% prevalence |
| Microdiverse_bp_per_100kbp_on_reads | Microdiverse_bp_on_reads normalized per 100kbp |
| Microdiverse_bp_on_reference | Bases present in reference but not in reads at >= 10% prevalence |
| Microdiverse_bp_per_100kbp_on_reference | Microdiverse_bp_on_reference normalized per 100kbp |

## Side misassembly

Evaluates completeness at contig extremities based on significant clipping events (>= 50% prevalence) where the median clipped length extends beyond the contig boundary.

| Metric | Definition |
| --- | --- |
| Coverage_first_position | Raw coverage (read depth) at the first position of the contig (before RLE compression) |
| Contig_start_collapse_prevalence | Prevalence of the terminal left-clipping event (clippings / local coverage), i.e. fraction of reads that are clipped at the start |
| Contig_start_collapse_bp | Median clipped length at the terminal left-clipping event (estimated missing sequence at start) |
| Contig_start_expansion_bp | Distance from the terminal left-clipping event to the contig start (contamination length at start) |
| Coverage_last_position | Raw coverage (read depth) at the last position of the contig (before RLE compression) |
| Contig_end_collapse_prevalence | Prevalence of the terminal right-clipping event (clippings / local coverage), i.e. fraction of reads that are clipped at the end |
| Contig_end_collapse_bp | Median clipped length at the terminal right-clipping event (estimated missing sequence at end) |
| Contig_end_expansion_bp | Distance from the terminal right-clipping event to the contig end (contamination length at end) |
| Contig_end_misjoint_mates | Number of read pairs where one mate maps to the contig end and the other maps to a different contig (paired-read only) |
| Normalized_contig_end_misjoint_mates | Contig_end_misjoint_mates normalized by coverage mean (paired-read only) |

**Warnings:**

- Side completeness metrics tend to **underestimate the incompleteness of the contig**, especially when using short reads. Short reads cannot extend far beyond the last properly assembled base, so collapse prevalence and expansion bp are typically slightly lower than the true values. The **collapse bp**, in particular, can be substantially smaller than the actual missing sequence.

- For circular contigs, completeness values can sometimes **overestimate missing sequence**. This typically occurs when mappings are performed without theBIGbam or without the `--circular` option, **because reads that would normally span the junction between the last and first contig positions are clipped: standard mappers do not support circular alignments**.

- In **rare cases**, this effect can also occur for long reads even when the `--circular` option is used. For example, in bacteriophages, single long reads may contain multiple full viral genomes due to concatemers. In such cases, only the portion of the read corresponding to one contig length (or up to twice the contig length when mapping with theBIGbam and the `--circular` option) is mapped, and the remaining sequence is clipped. If many concatemer-derived long reads are present, incompleteness metrics may detect this clipped signal, leading to artificially high estimates of missing sequence.

## Topology

Evaluates the circularity of contigs. Circularity can only be reliably assessed when the `--circular` option is used.

| Metric | Definition |
| --- | --- |
| Circularising_reads | Number of reads whose alignment starts near the end of the contig and continues past the contig start |
| Circularising_reads_prevalence | Circularising reads divided by the mean coverage at the first and last contig positions |
| Circularising_inserts | Number of read pairs spanning the junction between contig end and start (paired-read only) |
| Circularising_insert_size_deviation | Median insert size of circularising pairs minus median insert size of all proper pairs (paired-read only) |
| Normalized_circularising_inserts | Circularising_inserts normalized by coverage mean (paired-read only) |

An approximate estimate of circularity could be derived from mappings performed without the `--circular` option, using supplementary alignments spanning both contig ends and, for paired-end short reads, non-inward read pairs with reasonable insert sizes mapping to both contig ends. However, this approach does not provide an exact estimation of circularity compared to circular mapping with the `--circular` option. For this reason, circularity is not estimated in this case.

## Termini

Phage packaging mechanism detection based on terminus analysis.

| Metric | Definition |
| --- | --- |
| Packaging_mechanism | Detected packaging mechanism type (e.g. headful, cos, DTR, etc.) |
| Left_termini | Comma-separated center positions of kept left (start) terminus peaks |
| Right_termini | Comma-separated center positions of kept right (end) terminus peaks |
| Duplication | Type of terminal repeat: DTR (direct) or ITR (inverted), or NULL if no repeat detected |
| Total_peaks | Total number of kept terminus peaks (both left and right) |
| Repeat_length | Length of the detected terminal repeat in base pairs |
| Terminase_distance | Minimal distance (bp) from any kept terminus center to the nearest terminase gene annotation |
| Terminase_percentage | Terminase_distance expressed as a percentage of contig length |
