# Database command control

### Updating the database

Example command:

```sh
thebigbam add-variable examples/outputs/HK97/HK97.db test bars "#f60820" "Test Bar" examples/inputs/HK97/variable_test.csv
thebigbam add-variable examples/outputs/HK97/HK97.db test2 curve "#aef1c2" "Test Curve" examples/inputs/HK97/variable_test2.csv
```

Database compression

For each sample, Rust computations are performed for all contigs detected as present in the BAM file (i.e., those exceeding the --min_aligned_fraction and --min_coverage_depth thresholds). After computation, the final database is further reduced through a compression step. Coverage metrics are compressed using run-length encoding (RLE): consecutive positions with similar values (within the --compress_ratio threshold) are stored as a single entry, preserving the overall signal while drastically reducing storage size. For phagetermini and assemblycheck metrics, values are compared to the local coverage and discarded if they fall below the --compress_ratio threshold, ensuring that only meaningful peaks are retained. When consecutive positions remain after filtering, RLE is applied to group them.

The output is a DuckDB database containing all computed values. This database is typically ~100× smaller than the original BAM files while retaining the essential characteristics of the mapping data.

### Listing database contents

```sh
thebigbam list-samples -d my_database.db
thebigbam list-contigs -d my_database.db
thebigbam list-variables -d my_database.db          # add --detailed for full metadata
thebigbam list-sample-metadata -d my_database.db     # user-added columns on Sample table
thebigbam list-contig-metadata -d my_database.db     # user-added columns on Contig table
```

### Adding and removing metadata

Sample and contig tables can be extended with user-defined metadata columns from CSV files:

```sh
thebigbam add-sample-metadata -d my_database.db --csv metadata.csv
thebigbam add-contig-metadata -d my_database.db --csv metadata.csv
```

To remove a user-added metadata column:

```sh
thebigbam remove-sample-metadata -d my_database.db --colname my_column
thebigbam remove-contig-metadata -d my_database.db --colname my_column
```

Custom variables (added via `add-variable`) can be removed with:

```sh
thebigbam remove-variable -d my_database.db --name my_variable
```

### Removing samples and contigs

Entire samples or contigs can be removed from the database along with all their associated data:

```sh
thebigbam remove-sample -d my_database.db --name sample_to_remove
thebigbam remove-contig -d my_database.db --name contig_to_remove
```

These commands delete the sample/contig row and all rows referencing it across every table in the database (coverage, features, packaging, etc.).

**Effect on metrics:** All derived metrics such as RPKM, TPM, and normalized counts (per 100kbp) are computed on-the-fly through SQL views — they are not stored as fixed values. This means that removing a sample or contig does **not** make any remaining metric values incorrect. The views automatically recalculate from whatever data remains in the database.

### Extending an existing database (`--extend`)

For large datasets, it may be practical to build the database incrementally rather than processing all BAM files in a single run. The `--extend` flag allows adding new samples (and optionally new contigs) to an existing database:

```sh
# Initial database creation
thebigbam calculate -b first_batch/ -g contigs.gbk -o my_database.db

# Add more samples later
thebigbam calculate --extend -b second_batch/ -o my_database.db

# Add samples with new contigs
thebigbam calculate --extend -b third_batch/ -g new_contigs.gbk -o my_database.db
```

**Important considerations when using `--extend`:**

- **Contig name consistency:** Contig names must be identical between former and new BAM mappings. If the same biological contig has different names across your BAM files, results will be split across separate entries.

- **Module inheritance:** The modules computed for the existing database (Coverage, Misalignment, etc.) are automatically reused. Any `--modules` flag provided during extend is ignored — you cannot change which modules are computed after database creation.

- **Sample name collisions:** If a new BAM file has the same stem name as an existing sample, `--extend` will refuse to proceed. If you want to replace a sample with updated mappings, first remove it with `thebigbam remove-sample` and then re-run `--extend`.

- **Custom variables and metadata:** If you previously added custom variables (via `add-variable`) or metadata columns (via `add-sample-metadata` / `add-contig-metadata`) for existing samples or contigs, you will need to add those for the newly added data as well. They are not automatically propagated.

- **New contigs and existing samples:** When new contigs are introduced during an extend operation, existing samples will **not** have mapping data for those contigs. If you need coverage data for all samples on all contigs, you should remove the affected samples with `thebigbam remove-sample` and re-add them with BAM files that were mapped against the complete set of contigs (old + new).

- **New contig handling depends on database origin:**
  
  - If the database was created with annotation files (`-g`/`-a`), new contigs are only added if you also provide `-g`/`-a` during extend.
  - If the database was created in BAM-only mode (no annotation files), new contigs found in BAM headers are automatically added.
