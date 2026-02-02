# Database command control

### Updating the database

Example command:

```sh
thebigbam add-variable examples/outputs/HK97/HK97.db test bars "#f60820" "Test Bar" examples/inputs/HK97/variable_test.csv
thebigbam add-variable examples/outputs/HK97/HK97.db test2 curve "#aef1c2" "Test Curve" examples/inputs/HK97/variable_test2.csv
```

Database compression

For each sample, Rust computations are performed for all contigs detected as present in the BAM file (i.e., those exceeding the --min_coverage threshold). After computation, the final database is further reduced through a compression step. Coverage metrics are compressed using run-length encoding (RLE): consecutive positions with similar values (within the --compress_ratio threshold) are stored as a single entry, preserving the overall signal while drastically reducing storage size. For phagetermini and assemblycheck metrics, values are compared to the local coverage and discarded if they fall below the --compress_ratio threshold, ensuring that only meaningful peaks are retained. When consecutive positions remain after filtering, RLE is applied to group them.

The output is an SQLite database containing all computed values. This database is typically ~1000× smaller than the original BAM files while retaining the essential characteristics of the mapping data.

TODO: explain possibility to add, remove, list variables, contigs-(metadata), samples(-metadata)
