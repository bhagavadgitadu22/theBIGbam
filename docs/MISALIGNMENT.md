# Codon Change Annotation for Mismatches

When an annotation file is provided, each mismatch position inside a CDS is annotated with its codon impact: whether the base substitution is **synonymous** (same amino acid) or **non-synonymous** (different amino acid). Positions outside any CDS are labelled **Intergenic**.

## Position-summary codon analysis

Codon changes are computed from **per-position mismatch summaries** during post-processing, rather than per-read during BAM traversal. This is an O(positions_with_mismatches) operation instead of O(reads_with_mismatches), which keeps runtime proportional to genome size rather than sequencing depth.

### Algorithm

After all reads for a contig have been processed:

1. **Identify mismatch positions**: iterate positions where `mismatch_base_counts[pos]` has any non-zero entry (i.e. at least one read carried a mismatch here).
2. **Find dominant base**: at each position, the alternative base (A, C, G, or T) with the highest count is selected. It must exceed the prevalence threshold (relative to primary read depth) to be reported.
3. **Map to CDS**: the position is looked up against a sorted index of CDS intervals (binary search). Positions not in any CDS are recorded as Intergenic.
4. **Build the mutant codon**: starting from the reference codon (extracted from the CDS nucleotide sequence), the dominant mismatch base is substituted at the appropriate codon position. For reverse-strand CDS, the base is complemented before substitution.
5. **Translate and classify**: both reference and mutant codons are translated using the standard genetic code. If the amino acid is unchanged, the change is Synonymous; otherwise Non-synonymous.

### Per-position independence

Each mismatch position is evaluated independently. In the rare case where two mismatches in the same codon originate from the same read (e.g., two adjacent SNPs), the tool evaluates each position separately rather than combining them into a single mutant codon. This affects <0.01% of mismatches and has negligible impact on synonymous/non-synonymous classification, since the dominant mismatch at each position is independent of neighbouring positions.

### Dominant selection

Both `Sequence` and `Codon_change` are derived from the same dominant base at each position:

- **`Sequence` / `Sequence_prevalence`** (tooltip: Sequence / Prevalence): the dominant **nucleotide** at this position, selected by counting how many reads carry each alternative base.

- **`Codon_change` / `AA_change`** (tooltip: Codon / Amino acid): the codon containing the dominant nucleotide substitution at this position, with the rest of the codon taken from the reference sequence.

Because both use the same dominant base, `Sequence` and `Codon_change` are always consistent.

## Database columns

The following columns are added to the `Feature_mismatches` table:

| Column                | Tooltip label | Description                                                                |
| --------------------- | ------------- | -------------------------------------------------------------------------- |
| `Sequence`            | Sequence      | Dominant alternative nucleotide at this position                           |
| `Sequence_prevalence` | Prevalence    | Percentage of reads carrying that nucleotide (relative to coverage)        |
| `Codon_category`      | Category      | `Synonymous`, `Non-synonymous`, or `Intergenic`                            |
| `Codon_change`        | Codon         | The dominant mutant codon (e.g. `ACG`). NULL for Intergenic positions      |
| `AA_change`           | Amino acid    | The resulting amino acid, e.g. `V (Valine)`. NULL for Intergenic positions |

Codon columns are NULL for mismatch positions where no dominant base was identified (i.e. positions stored without a `Sequence` value).

## CDS index

CDS intervals are extracted from the annotation file during parsing. For each contig, a sorted list of CDS intervals (start, end, strand, nucleotide sequence) is built. Lookup uses binary search on CDS start positions, scanning backwards to handle overlapping/nested CDS features correctly.

When multiple isoforms share the same locus tag, only the **longest isoform** has its nucleotide sequence computed and is included in the CDS index for codon analysis. All isoforms are still displayed in the gene map, but codon/amino acid annotations are derived from the longest one only. This simplifies visualization by avoiding conflicting or redundant annotations at overlapping positions.
