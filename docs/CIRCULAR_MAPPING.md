# On mapping with circular genome support

Circular genome support during mapping allows reads to span contig ends, which is important for assessing contig completeness and topology. This mapping benefits all complete circular elements, particularly smaller contigs, where terminal regions constitute a larger fraction, and long-read datasets, where reads are more likely to span contig boundaries.

## Principle

Reads spanning contig ends are typically splitted into a primary and a chimeric supplementary alignment with standard mappers. The supplementary alignment is then often discarded by default.

With `thebigbam mapping-per-sample --circular` mode, each contig is duplicated prior to alignment, enabling seamless mapping across the junction. Artificial secondary and supplementary alignments arising from the duplication are removed. Finally, a modulo operation reassigns reads to their correct positions on the original contig before generating the output BAM files.

This approach preserves consistent coverage at contig ends of circular genomes, whereas standard mappers produce a coverage drop. For linear genomes, circular mapping only marginally affects terminal bases. The probability that a read would randomly extend ≥10 bp beyond a true contig end is extremely low. Under a random nucleotide model, a 10 bp match is expected with a frequency of (1/4)¹⁰ (~10⁻⁷), making boundary effects negligible in practice.

## Impact of circular mapping on the HK97 genome

We mapped both short and long reads to the 39.7 kbp genome of the model bacteriophage HK97 (NC_002167.1) using `thebigbam mapping-per-sample` command in default and circular modes to assess the difference. The reference genome was rotated by 1 kbp to prevent the under-covered COS 3′ termini from appearing at the contig ends.

Circular mapping approximately doubled coverage at both contig ends. At the 3′ end, circular mapping increased coverage from 100 to 197 for short reads and from 2,482 to 6,337 for long reads. The region of improved coverage spanned approximately half the average read length (e.g. 75 bp for short reads) from each contig end. Circular mapping using long reads increased the global mean coverage from 5,000 to 5,700.

## SAM specification for circular genomes

The output BAM follows the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) (section 2.7) for circular references:

> *"As usual POS should be between 1 and the @SQ header’s LN value, but POS plus the sum of the lengths of M/=/X/D/N CIGAR operations may exceed LN. Coordinates greater than LN are interpreted by subtracting LN so that bases at LN+1, LN+2, LN+3, . . . are considered to be mapped at positions 1, 2, 3, . . .; thus each (1-based) position p is interpreted as ((p − 1) mod LN) + 1."*

This means reads that cross the origin of a circular genome will have their alignment extend past the contig length in the BAM file, which is valid per the specification.

**Compatibility warning:** although this encoding is SAM-spec-compliant, most downstream tools (IGV, samtools depth, variant callers, etc.) do not yet implement circular coordinate handling. Programs that are unaware of this part of the specification will likely produce incorrect results for reads wrapping around the origin. Verify that your downstream tool explicitly supports circular genomes before relying on the output outside of theBIGbam.

## MAPQ limitation

in circular mode the mapper sees the doubled reference as two distinct loci, so every read produces at least one ghost alignment that competes with the real one. The mapper therefore assigns lower MAPQ scores than the same read would receive against a single-copy reference, because MAPQ quantifies the confidence that the *reported* alignment is the only good one — and in a doubled reference there are always at least two.

The circular-BAM conversion step removes these ghost secondary and supplementary alignments, restoring the correct read count. For reads whose *only* competing alignments were ghosts (i.e. after removal they have no remaining secondary or supplementary records), the conversion sets MAPQ to 60, which is the conventional value for a uniquely mapping read. Reads that retain legitimate secondary or supplementary alignments after ghost removal keep their original MAPQ unchanged.

Full MAPQ recalculation is not possible because both minimap2 and bwa-mem2 derive MAPQ from internal values that are never written to BAM output:

- **minimap2** uses chaining scores (*f*₁, *f*₂) and anchor counts to compute MAPQ; these are discarded after alignment.
- **bwa-mem2** uses chaining scores, seed coverage, and the repetitive fraction (*frac_rep*); these are likewise internal only.

The common practice of using MAPQ values to filter alignments makes this limitation frustrating. However, the impact is mitigated by the fact that MAPQ scores are inherently software-dependent and should always be interpreted with caution.
