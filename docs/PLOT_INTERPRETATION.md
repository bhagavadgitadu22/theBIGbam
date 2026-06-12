# How to interpret the plots

This guide provides practical tips for interpreting the visualizations produced by theBIGbam. For a complete reference of all features, see [FEATURES.md](FEATURES.md).

---

## On misassembly and microdiversity

The misalignment module flags positions where mapped reads diverge substantially from the reference. These positions may reflect subpopula-tions missing a segment of the contig or carrying an extra segment. Dis-crepancies supported by more than 50% of locally aligned reads can also be considered misassemblies when they occur in the sample from which the contig was assembled.
To quantify microdiversity and misassembly per contig per sample, base pairs present in reads but absent from the reference (collapses) and base pairs present in the reference but absent in reads (expansions) are computed at all flagged positions (microdiversity level) and at positions where > 50% of reads support the discrepancy (misassembly level) (Guo et al. 2025). Collapse size at each position is calculated as the average length of insertion and clipping events. The sum of collapses and expan-sions gives contig-level analogs of completeness and contamination.
These metrics should be interpreted as approximations rather than exact measurements, as mapping files only provide contextual information up to the read length (typically ~150 bp for short reads). For instance, a clip-ping event may reflect a small (e.g. ~200 bp) or a large (e.g. >10 kbp) missing segment, yet the inferred length is inherently capped by read length (e.g. ~100 bp for short reads). As a result, those metrics tend to underestimate the size of missing or extra genomic segments.

Complementary information from other modules can further confirm misassemblies or biologically meaningful subpopulations (microdiversity):

- For the **Coverage module**: local peaks and troughs in primary reads depth ; presence of secondary or supplementary alignments suggests repeated regions within a contig or similarity to another contig

- For the **Paired-Reads module**: increase in read insert size indicates that genomic segments may be missing from the reference, whereas a sharp decrease suggests additional sequence in some reads. Non-inward read pairs may indicate inversions, while missing mates point to incomplete assembly. Reads whose mates map to a different contig suggest a relationship between contigs, either due to shared genomic regions or because one represents the continuation of the other in the organism’s genome.



## Coverage vs Coverage Reduced

These two coverage tracks use different counting rules and can yield different values at the same position:

- **Coverage** counts a position as covered only if the read has a true base match at that position (based on CIGAR/MD). Clipped or mismatched bases are not counted.

- **Coverage Reduced** (from the Phage termini module) counts reads that start *and* end with an exact match (no clipping at their termini). For qualifying reads, every position between the first and last aligned base is counted — including positions that may have internal mismatches.

As a result, at any given position either metric can be higher than the other:

| Scenario                              | Which is higher?                       | Likely cause                                                                                      |
| ------------------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------- |
| Coverage > Coverage Reduced           | Coverage Reduced is lower              | Many reads have clipped ends at this position, suggesting the assembly may be missing sequence    |
| Coverage Reduced > Coverage           | Coverage is lower                      | Some reads have internal mismatches at this position, suggesting microdiversity in the population |
| Coverage Reduced << Coverage globally | Coverage Reduced much lower everywhere | Reads were not properly trimmed before mapping — adapter sequences likely remain                  |

---

## Phage packaging signals

When analyzing phage genomes for packaging mechanisms, look at the **Termini module** features:

- **Read starts / Read ends**: Sharp, narrow peaks indicate biological cut sites. Broad or noisy signals suggest random fragmentation.
- **Tau (τ)**: Normalizes terminus counts by coverage. Values above 0.5 suggest enrichment; values approaching 2 indicate very strong terminus signals.
- **Coverage Reduced**: Provides a cleaner baseline for terminus detection by excluding reads with noisy ends.

For detailed information about the classification algorithm and output format, see [PHAGE_PACKAGING.md](PHAGE_PACKAGING.md).

---

## Assembly quality indicators

The **Misalignment module** features help identify problematic regions:

- **Clippings** (left/right): Peaks indicate sequence present in reads but not in the reference. Common at contig ends if the assembly is incomplete, but can also indicate adapter contamination.
- **Insertions / Deletions**: May represent true structural variants, assembly errors, or sequencing artifacts (especially in homopolymer regions).
- **Mismatches**: Can indicate SNPs, sequencing errors, or microdiversity in the population.

Use the **"Plot relative to local coverage"** checkbox to normalize these counts. This reveals whether anomalies are proportional to coverage (likely noise) or represent a consistent signal regardless of depth (likely biological or assembly-related).

---

# 

 
