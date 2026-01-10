# Assembly Check Metrics

This document explains the metrics available in the **Per module:** filtering section of the visualization interface. These metrics help you filter contigs and samples based on coverage patterns and assembly quality indicators.

## Coverage Metrics

These metrics are computed for each contig-sample pair and describe the overall sequencing coverage.

### Coverage Percentage

**What it measures:** The fraction of the contig that has at least one read aligned to it.

**How it's computed:**

```
Coverage Percentage = (Number of positions with coverage >= 1) / (Total contig length) × 100
```

**Interpretation:**

- **100%**: Every base of the contig has at least one read covering it
- **Lower values**: Some regions have no coverage, which may indicate assembly gaps, deletions, or insufficient sequencing depth

### Coverage Mean

**What it measures:** The average sequencing depth across all positions of the contig.

**How it's computed:**

```
Coverage Mean = Sum of coverage at all positions / Contig length
```

**Interpretation:**

- Higher values indicate deeper sequencing
- Useful for filtering out low-depth samples that may have unreliable feature detection
- For example, filtering for Coverage Mean > 50 ensures you're only viewing contigs with robust sequencing support

### Coverage Variation

**What it measures:** How uniform the coverage is along the contig. This metric captures the average change in coverage between adjacent positions.

**How it's computed:**

```
Coverage Variation = (1 / (n-1)) × Σ|coverage(i+1) - coverage(i)|
```

Where n is the contig length and the sum runs over all consecutive position pairs.

**Interpretation:**

- **Low values** (close to 0): Very uniform coverage - the sequencing depth is consistent along the contig
- **High values**: Highly variable coverage - there are significant peaks and valleys in the coverage profile
- Uniform coverage is generally expected for well-assembled contigs; high variation may indicate:
  - Repetitive regions with mapping artifacts
  - GC bias in sequencing
  - Structural variants or assembly errors

---

## Completeness Metrics

These metrics assess assembly quality by analyzing soft-clipped reads (reads that only partially align to the reference). They help identify potential missing sequence at contig ends and internal assembly issues.

### Prevalence Completeness Left / Right

**What it measures:** An estimate of how complete the contig ends are, based on the proportion of reads that are soft-clipped at the terminus.

**How it's computed:**

1. Find the first (outermost) clipping event at each end
2. Check if the median clipped length exceeds the distance to the contig end
3. If yes, compute: `Prevalence = (Number of clipped reads / Total coverage) × 100`
4. The completeness is then: `100 - Prevalence`

**Interpretation:**

- **100%**: No significant clipping at this end - the contig end appears complete
- **Lower values**: A substantial fraction of reads are clipped, suggesting the contig may be truncated and missing sequence beyond this point
- This metric is only computed when the clipping pattern is consistent with a truncated assembly (clipped length > distance to end)

### Percentage Completeness

**What it measures:** A global score estimating what fraction of the "true" sequence is present in the assembly, based on features that suggest missing sequence in reads.

**How it's computed:**

```
Score = Total_mismatches + Total_insertions + Total_reads_clipped
Percentage Completeness = 100 - (Score × 100 / Contig_length)
```

Where each component is computed as the weighted sum across all positions:

- **Mismatches**: Count of mismatched bases normalized by coverage
- **Insertions**: Insertion events weighted by their median length, normalized by coverage
- **Reads clipped**: Soft-clipping events weighted by median clipped length, normalized by coverage

**Interpretation:**

- **100%**: No evidence of missing sequence - reads align perfectly
- **Lower values**: Substantial insertions, mismatches, or clippings suggest the reference may be missing sequence that exists in the reads
- This estimates sequence present in reads but absent from the reference

### Percentage Contamination

**What it measures:** An estimate of how much extra sequence exists in the reference that is not supported by the reads.

**How it's computed:**

```
Score = Total_mismatches + Total_deletions + Total_reference_clipped
Percentage Contamination = (Score × 100 / Contig_length)
```

Where:

- **Mismatches**: Same as above
- **Deletions**: Deletion events weighted by their median length, normalized by coverage
- **Reference clipped**: Paired clipping patterns (right-clip followed by left-clip) that suggest an inserted region in the reference

**Interpretation:**

- **0%**: No evidence of extra sequence - the reference is well-supported by reads
- **Higher values**: The reference may contain sequence not present in the sequenced sample, which could indicate:
  - Contamination in the reference assembly
  - Strain variation (deletions in the sequenced strain)
  - Assembly artifacts

---

## Phage Mechanism Filter

**What it measures:** Filters contigs by their detected DNA packaging mechanism (for phage genomes).

**Available options:** Depends on what was detected during processing. Common mechanisms include:

- **Headful**: Packaging by head-full mechanism (e.g., P22-like phages)
- **Cohesive ends**: Packaging with specific cohesive termini (e.g., lambda-like)
- **DTR** (Direct Terminal Repeats): Terminal redundancy at both ends
- **Host**: No clear phage packaging signal detected

See the [Phage Packaging documentation](PHAGE_PACKAGING.md) for detailed information about how packaging mechanisms are detected.

---

## Using Filters Effectively

### Finding high-quality data

```
Coverage Percentage: 95-100%
Coverage Mean: > 50
Coverage Variation: < 20
Percentage Completeness: > 95%
Percentage Contamination: < 5%
```

### Finding potentially problematic assemblies

```
Percentage Completeness: < 80%
OR
Percentage Contamination: > 10%
```

### Filtering for uniform coverage

```
Coverage Variation: 0-10
```

This is useful when looking for consistent signal across the genome, such as for accurate terminus detection.
