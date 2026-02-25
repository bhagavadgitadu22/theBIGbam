# Compression & Memory Optimizations

## Dominant Sequence Tracking

For each assembly-check feature (insertions, left/right clips, start/end clips), theBIGbam tracks the **dominant sequence** — the most frequent sequence variant at each genomic position. This is used to populate the `sequence` column in the database.

### Memory bounds

Two caps keep memory usage bounded during the BAM pass:

1. **Max 10 unique sequence variants per position** (`MAX_SEQS_PER_POS = 10`).
   At any position, only the first 10 distinct sequences are stored. Once 10 variants are tracked, new unique sequences are dropped. However, **counts for already-tracked sequences continue to increment without limit**. This means the dominant sequence count remains accurate as long as the dominant variant appears within the first 10 unique sequences seen — which is virtually guaranteed for biological data, where the true dominant sequence appears early and frequently.

2. **Insertion sequences capped to 43 bytes** (first 20bp + `...` + last 20bp).
   Long-read insertions can be thousands of base pairs. Sequences longer than 40bp are truncated to the first 20bp and last 20bp, separated by `...` to indicate truncation. The full insertion *length* is still recorded separately in `insertion_lengths` and is unaffected by this cap.

### Prevalence calculation

The prevalence (stored as `percentage_x10` in the database) is computed as:

```
prevalence = dominant_count / total_primary_reads_at_position
```

The denominator is the total number of primary reads covering that position (from the coverage array), **not** the sum of tracked variant counts. This means the 10-variant cap does not inflate prevalence — if the dominant sequence appeared in 800 out of 1000 reads, prevalence is 0.80 regardless of how many other variants were tracked or dropped.

### Allocation avoidance

The `track_sequence` helper accepts a byte slice (`&[u8]`) rather than an owned `Vec<u8>`. A heap allocation (`to_vec()`) only occurs when a genuinely new variant needs to be inserted into the map. In the common case — the sequence is already tracked (just increment) or the map is full (just drop) — no allocation happens. This avoids millions of unnecessary per-read allocations during the BAM pass.
