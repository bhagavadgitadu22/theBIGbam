# theBIGbam visualisation

## Hidden rules

When multiple isoforms share the same locus tag, only the **longest isoform** has its nucleotide sequence computed and is included in the CDS index for codon analysis. All isoforms are still displayed in the gene map, but codon/amino acid annotations are derived from the longest one only. This simplifies visualization by avoiding conflicting or redundant annotations at overlapping positions.

For MAG dot colors and gene map features, when multiple coloring rules apply to the same element, the first matching rule is applied. If a feature has multiple annotations (e.g., CRISPR and restriction–modification activity), all values are shown in labels (e.g., `CRISPR^RM` when using `activity`). However, color assignment (e.g., “Use random colors”) uses a single category per feature, even if multiple annotations are present.
