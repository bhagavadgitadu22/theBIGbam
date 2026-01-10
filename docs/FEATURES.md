# # What features are inferred from your data?

This document describes all features computed by MGFeatureViewer and stored in the DuckDB database. Features are organized by **module** as they appear in the visualization interface.

6 modules exist at the moment: 5 are defined from the mappings files and the Genome module is defined from the optional annotation file.

---

## Coverage module

Features describing read alignment depth and quality across the genome.

### Subplot: Primary Alignments

#### Primary Reads

- **Description:** The number of reads covering each position, counting only primary alignments
- **How it's computed:** For each position, count reads where this is the best alignment (not secondary or supplementary). This is the standard "coverage" or "read depth"
- **Display:** Shown as a continuous curve

### Subplot: Alignments by Strand

#### Primary Reads (+ only)

- **Description:** Coverage from reads aligned to the forward strand only
- **How it's computed:** Count of primary reads at each position where the read maps to the plus strand
- **Use case:** Strand bias can indicate certain library preparation artifacts or biological features like transcription

#### Primary Reads (- only)

- **Description:** Coverage from reads aligned to the reverse strand only
- **How it's computed:** Count of primary reads at each position where the read maps to the minus strand
- **Use case:** Compare with forward strand to detect strand imbalances

### Subplot: Other Alignments

#### Secondary Reads

- **Description:** Count of secondary alignments at each position
- **How it's computed:** Reads flagged as secondary (SAM flag 0x100) - alternative alignments when a read maps to multiple locations
- **Interpretation:** High secondary alignment counts indicate repetitive or ambiguous regions where reads could map to multiple places

#### Supplementary Reads

- **Description:** Count of supplementary alignments at each position
- **How it's computed:** Reads flagged as supplementary (SAM flag 0x800) - chimeric alignments where different parts of the read map to different locations
- **Interpretation:** High supplementary counts may indicate structural variants, chimeric sequences, or assembly errors

### Subplot: MAPQ

#### Mapping Quality

- **Description:** Average confidence of read alignments at each position
- **How it's computed:** `Average MAPQ = Sum of MAPQ values / Number of primary reads`
- **Value range:** 0-60 (typically), where higher values indicate more confident alignments
- **Interpretation:**
  - **High MAPQ (>30):** Reads align uniquely and confidently
  - **Low MAPQ (<10):** Reads could map to multiple locations; the alignment is ambiguous
  - Note: MAPQ scoring varies between aligners (BWA, Bowtie2, minimap2, etc.)

---

## Long-read Metrics Module

Features specific to long-read sequencing data (PacBio, Nanopore).

### Subplot: Read Lengths

#### Read Lengths

- **Description:** Average length of reads covering each position
- **How it's computed:** `Average = Sum of read lengths at position / Count of reads at position`
- **Use case:** Identify regions covered by shorter or longer reads. Unusually short reads in a region might indicate fragmentation or alignment issues

---

## Paired-read Metrics Module

Features specific to paired-end/mate-pair sequencing data (Illumina).

### Subplot: Insert Sizes

#### Insert Sizes

- **Description:** Average distance between read pairs at each position
- **How it's computed:** For properly paired reads, compute `Average = Sum of insert sizes / Count of proper pairs`
- **Interpretation:**
  - Consistent insert sizes indicate normal library structure
  - Deviations from the expected insert size may indicate structural variants (insertions/deletions)

### Subplot: Non-inward Pairs

#### Non-inward Pairs

- **Description:** Count of reads where the mate is on the same contig but in an unexpected orientation
- **How it's computed:** Count reads where the mate maps to the same contig but the pair orientation is not the expected inward-facing (FR) orientation
- **Interpretation:** High counts suggest:
  - Inversions in the sample relative to the reference
  - Tandem duplications
  - Assembly errors

### Subplot: Mate Not Mapped

#### Missing Mates (unmapped)

- **Description:** Count of reads whose mate failed to align anywhere
- **How it's computed:** Count reads where the mate unmapped flag (0x8) is set
- **Interpretation:** High counts may indicate:
  - Sequence not present in the reference
  - Poor quality mate reads
  - Contamination in the library

### Subplot: Mate on Another Contig

#### Missing Mates (other contig)

- **Description:** Count of reads whose mate aligned to a different contig
- **How it's computed:** Count reads where the mate reference ID differs from the read's reference ID
- **Interpretation:** Can indicate:
  - Chimeric molecules in the library
  - Misassemblies where contigs should be joined
  - Mobile elements or prophages integrated at this position

---

## Mapping Metrics Per Position Module

Features describing alignment anomalies that may indicate assembly issues or biological variation.

### Subplot: Clippings

Soft/hard clipping occurs when part of a read does not align to the reference. The clipped portion represents sequence in the read that has no corresponding match in the reference.

#### Left Clippings

- **Description:** Count and length of clipped bases at the 5' end of reads
- **How it's computed:** At each position, count reads with soft/hard clipping at their left (5') end. Also computes mean, median, and standard deviation of clipping lengths
- **Interpretation:** Indicates sequence present in reads but missing from the left side of the reference at this position. Common at contig ends if the assembly is incomplete

#### Right Clippings

- **Description:** Count and length of clipped bases at the 3' end of reads
- **How it's computed:** At each position, count reads with soft/hard clipping at their right (3') end. Also computes mean, median, and standard deviation of clipping lengths
- **Interpretation:** Indicates sequence present in reads but missing from the right side of the reference at this position

### Subplot: Indels

#### Insertions

- **Description:** Count and length of insertion events at each position
- **How it's computed:** From CIGAR strings, count 'I' operations. Also computes mean, median, and standard deviation of insertion lengths
- **Interpretation:** Sequence present in the reads but absent from the reference. Could indicate:
  - True insertions in the sequenced sample
  - Missing sequence in the reference assembly
  - Sequencing errors (especially in homopolymer regions)

#### Deletions

- **Description:** Count of deletion events at each position
- **How it's computed:** From CIGAR strings, count 'D' operations
- **Interpretation:** Sequence present in the reference but absent from reads. Could indicate:
  - True deletions in the sequenced sample
  - Extra sequence incorrectly included in the reference
  - Alignment artifacts

### Subplot: Mismatches

#### Mismatches

- **Description:** Count of base substitutions at each position
- **How it's computed:** From the MD tag in BAM files, count positions where the read base differs from the reference base
- **Interpretation:** Could indicate:
  - SNPs (true variation between sample and reference)
  - Sequencing errors
  - Alignment errors in repetitive regions

---

## Phage Termini Module

Features designed for detecting phage DNA packaging sites and terminus types. These are particularly useful for analyzing bacteriophage genomes.

### Subplot: Reads Termini

These features count where reads physically start and end, which can reveal DNA packaging cut sites.

#### Read Starts

- **Description:** Count of read 5' ends at each position
- **How it's computed:** Count how many reads have their first aligned base at each position. For paired-end data, uses strand-aware logic to identify true 5' ends
- **Interpretation:** Peaks indicate positions where DNA molecules frequently begin, which may represent:
  - Phage packaging initiation sites
  - DNA cutting/fragmentation sites
  - For random fragmentation (most libraries), should be relatively uniform

#### Read Ends

- **Description:** Count of read 3' ends at each position
- **How it's computed:** Count how many reads have their last aligned base at each position
- **Interpretation:** Peaks indicate positions where DNA molecules frequently end. Combined with Read Starts, helps identify terminus types

### Subplot: Coverage Reduced

#### Coverage Reduced

- **Description:** Coverage counting only "clean" reads without clipping or mismatches at their ends
- **How it's computed:** Count reads that:
  - Start with an exact match (no 5' clipping or mismatch)
  - For long reads: also end with an exact match
- **Use case:** Provides cleaner signal for terminus detection by excluding reads with noisy ends that could obscure true packaging sites
- **Display:** Shown as a continuous curve

### Subplot: Tau

#### Tau (τ)

- **Description:** Normalized measure of read terminus enrichment relative to coverage

- **How it's computed:**
  
  ```
  τ = (Read Starts + Read Ends) / Coverage Reduced
  ```
  
  Averaged over the positions in each compressed data run

- **Value range:** 0 to ~2 (can exceed 2 at very sharp termini)

- **Interpretation:**
  
  - **τ ≈ 0:** Few read termini relative to coverage (typical for internal regions)
  - **τ > 0.5:** Enrichment of read termini, suggesting a potential packaging site
  - **τ approaching 2:** Very strong terminus signal (both starts and ends enriched)

- **Use case:** Tau normalizes terminus counts by coverage, making it easier to identify true packaging sites that would otherwise be obscured by coverage variation

---

## Genome module

Features describing intrinsic genomic properties, independent of sequencing data.

### Subplot: Repeats

#### Direct Repeats

- **Description:** Regions of the genome that are repeated in the same orientation
- **How it's computed:** Detected by self-BLAST of the contig sequence, identifying segments that appear multiple times in the same 5'→3' direction
- **Biological relevance:** Terminal direct repeats (DTR) are characteristic of certain phage packaging mechanisms. Internal direct repeats may indicate mobile elements or gene duplications

#### Inverted Repeats

- **Description:** Regions of the genome that are repeated in opposite orientations
- **How it's computed:** Detected by self-BLAST, identifying segments where one copy is the reverse complement of another
- **Biological relevance:** Inverted terminal repeats (ITR) are found in certain phages and transposons. Internal inverted repeats can form secondary structures (hairpins)

## Database Storage Notes

- All feature values are stored as **INTEGER** in DuckDB for efficient compression
- **Tau** and **MAPQ** values are stored multiplied by 100 (e.g., τ=0.75 stored as 75)
- Features with length statistics (clippings, insertions) include additional **Mean**, **Median**, and **Std** columns
- The **Variable** table contains metadata for each feature including display properties (color, plot type, module assignment)

---

## See Also

- [Assembly Check Metrics](ASSEMBLY_CHECK.md) - Filtering metrics based on these features
- [Phage Packaging](PHAGE_PACKAGING.md) - How packaging mechanisms are detected from terminus patterns

This document describes all features computed by MGFeatureViewer and stored in the DuckDB database. Features are organized by **module** and **subplot** as they appear in the visualization interface.
