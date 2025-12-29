# Phage Packaging Classification Algorithm

The packaging mechanism of a contig is determined through a multi-step peak detection and classification process:

1. Signal Detection: For each position in the contig, we count read starts (5' ends) and read ends (3' ends) from the BAM alignment. These signals indicate where DNA molecules begin and end, which corresponds to packaging termini in phage genomes.

2. Bar Compression: Positions are filtered to retain only significant signals where count > coverage_reduced × bar_ratio. This removes noise at high-coverage regions while preserving signals that stand out relative to local coverage.

3. Clipping Filter: Peaks near significant soft-clipping events (≥ min_events) are discarded. Clippings indicate reads that don't fully align, often due to assembly artifacts rather than true termini. Left clippings filter read_starts; right clippings filter read_ends.

4. Two-Pass Peak Merging:
- First pass (circular-aware): Nearby peaks within max_distance_peaks (default 20bp) are merged, keeping the strongest. For circular genomes, peaks at position 1 and near the genome end are recognized as adjacent and merged. This happens on ORIGINAL data before any DTR region zeroing.
- Second pass (DTR-aware): For genomes with detected terminal repeats (from autoblast, ≥90% identity), equivalent positions in the first and second repeat regions are recognized as the same biological signal and merged.

5. Classification: Based on the final merged peaks:
- 0 peaks: No_packaging
- 1 peak (start or end only): PAC (headful packaging with single initiation site)
- 2 peaks (1 start + 1 end): Classified by distance between them:
    - Distance < 2bp: COS (cohesive ends, exact cut)
    - Distance ≤ 20bp: COS_3' or COS_5' (cohesive ends with overhang, direction indicates which strand protrudes)
    - Distance ≤ 1000bp: DTR_short_3' or DTR_short_5' (short direct terminal repeat)
    - Distance ≤ 10% genome: DTR_long_3' or DTR_long_5' (long direct terminal repeat)
    - Distance > 10% genome: DTR_outlier_3' or DTR_outlier_5'
    - For peaks in ITR (inverted repeat) regions: ITR_short, ITR_long, or ITR_outlier variants
- >2 peaks: Unknown_packaging

The _3' or _5' suffix indicates directionality: which terminus comes first when traversing the genome in the forward direction.