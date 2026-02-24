# How to interpret the plots

This guide provides practical tips for interpreting the visualizations produced by theBIGbam. For a complete reference of all features, see [FEATURES.md](FEATURES.md).

---

## Coverage vs Coverage Reduced

These two coverage tracks use different counting rules and can yield different values at the same position:

- **Coverage** counts a position as covered only if the read has a true base match at that position (based on CIGAR/MD). Clipped or mismatched bases are not counted.

- **Coverage Reduced** (from the Phage termini module) counts reads that start *and* end with an exact match (no clipping at their termini). For qualifying reads, every position between the first and last aligned base is counted — including positions that may have internal mismatches.

As a result, at any given position either metric can be higher than the other:

| Scenario | Which is higher? | Likely cause |
|---|---|---|
| Coverage > Coverage Reduced | Coverage Reduced is lower | Many reads have clipped ends at this position, suggesting the assembly may be missing sequence |
| Coverage Reduced > Coverage | Coverage is lower | Some reads have internal mismatches at this position, suggesting microdiversity in the population |
| Coverage Reduced << Coverage globally | Coverage Reduced much lower everywhere | Reads were not properly trimmed before mapping — adapter sequences likely remain |

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

## See Also

- [FEATURES.md](FEATURES.md) — Complete reference for all computed features
- [FILTERS.md](FILTERS.md) — Available filtering metrics for narrowing down contigs/samples
- [ASSEMBLY_CHECK.md](ASSEMBLY_CHECK.md) — Detailed explanation of assembly quality metrics
