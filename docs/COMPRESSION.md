# Compression and information loss

A few additional information losses occur during `thebigbam calculate`, beyond those described on the [main page](../README.md#database-compression).

## GC content and GC skew windowing

GC content is computed in 500 bp non-overlapping windows and GC skew in 1 kbp non-overlapping windows. Each stored value represents the window average, not a per-base value. Sub-window variation is not recoverable from the database. For most use cases — broad GC composition, horizontal gene transfer detection — this resolution is sufficient.

## Dominant variant sequences are truncated

For clippings, insertions, and mismatches, theBIGbam tracks which sequence variant occurs most frequently at each position. At most 10 unique sequence variants are tracked per position (`MAX_SEQS_PER_POS = 10`). Once 10 variants have been seen, new unique sequences are dropped, but counts for already-tracked variants continue to accumulate. The dominant sequence is determined at the end of the BAM pass as the variant with the highest count. In practice, the true biological dominant variant appears early and frequently, so the 10-variant cap rarely affects which variant is identified as dominant.

For clippings (left and right), sequences are truncated to the first 20 bp before tracking (`MAX_CLIP_SEQ_LEN = 20`). For insertions, sequences up to 40 bp are stored in full; longer insertions are truncated to the first 20 bp and the last 20 bp. Clipping and insertion lengths are recorded separately and are not truncated.

## Codon change approximation

For mismatch annotation, theBIGbam classifies variants as synonymous, non-synonymous, or intergenic using per-position mismatch summaries rather than individual read sequences. Each genomic position is evaluated independently: the dominant mismatch base observed at that position is substituted into the reference codon, and the resulting codon change is classified.

This approach does not preserve the linkage between mutations occurring on the same read. Consequently, if two or more mismatches occur within the same codon and are carried by the same read, they are not combined into a single mutant codon. As a result, the inferred amino acid change may occasionally differ from the true amino acid change present in the sequenced molecule.
