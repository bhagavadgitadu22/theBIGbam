# Database command control

Several commands are available to manipulate the database.

### Extending an existing database

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

- **Module inheritance:** The modules computed for the existing database (Coverage, Misalignment, etc.) are automatically reused. Any `--modules` flag provided during extend is ignored, you cannot change which modules are computed after database creation.

- **Sample name collisions:** If a new BAM file has the same stem name as an existing sample, `--extend` will refuse to proceed. If you want to replace a sample with updated mappings, first remove it with `thebigbam remove-sample` and then re-run `--extend`.

- **Custom variables and metadata:** If you previously added custom variables (via `add-variable`) or metadata columns (via `add-sample-metadata` / `add-contig-metadata`) for existing samples or contigs, you will need to add those for the newly added data as well. They are not automatically propagated.

- **New contigs and existing samples:** When new contigs are introduced during an extend operation, existing samples will **not** have mapping data for those contigs. If you need coverage data for all samples on all contigs, you should remove the affected samples with `thebigbam remove-sample` and re-add them with BAM files that were mapped against the complete set of contigs (old + new).

- **New contig handling depends on database origin:**
  
  - If the database was created with annotation files (`-g`/`-a`), new contigs are only added if you also provide `-g`/`-a` during extend.
  - If the database was created in BAM-only mode (no annotation files), new contigs found in BAM headers are automatically added.

## Adding content to the database

Sample, contig and MAG tables can be extended with user-defined metadata columns from CSV files:

```sh
thebigbam add-sample-metadata -d my_database.db --csv metadata.csv
thebigbam add-contig-metadata -d my_database.db --csv metadata.csv
thebigbam add-mag-metadata -d my_database.db --csv metadata.csv
```

To add sample metadata, metadata.csv has a header like Sample,Var1,Var2,... with Var1, Var2, etc. the new metadata columns. The same for contig and MAG metadata with Contig,Var1,Var2,... or MAG,Var1,Var2,...

New custom variables can also be added:

```sh
thebigbam add-variable --db my_database.db --name test --type bars --color "#f60820" --title "Test Bar" --csv variable_test.csv
thebigbam add-variable --db my_database.db --name test2 --type curve --color "#aef1c2" --title "Test Curve" --csv variable_test2.csv
```

The first example adds a new variable test that will be plotted using bars. The CSV containing the data must look like Contig,Sample,First_position,Last_position,Value. A row without Sample value is interpreted as a contig feature (like the GC content of a contig). The second example does the same but adding a variable that will be plotted as a curve.

## Removing content from the database

Entire samples or contigs can be removed from the database along with all their associated data:

```sh
thebigbam remove-sample -d my_database.db --name sample_to_remove
thebigbam remove-contig -d my_database.db --name contig_to_remove
thebigbam remove-mag -d my_database.db --name mag_to_remove
```

These commands delete the sample/contig/MAG row and all rows referencing it across every table in the database (coverage, features, packaging, etc.). Some metrics must be recomputed:

- `Number of samples` where a contig/MAG is found is recomputed after erasing a sample

- `Number of mapped reads` per sample is recomputed after removing a contig/MAG

To remove a user-added metadata column (added via `add-sample-metadata`, `add-contig-metadata` or `add-mag-metadata`):

```sh
thebigbam remove-sample-metadata -d my_database.db --colname my_column
thebigbam remove-contig-metadata -d my_database.db --colname my_column
thebigbam remove-mag-metadata -d my_database.db --colname my_column
```

Custom variables (added via `add-variable`) can be removed with:

```sh
thebigbam remove-variable -d my_database.db --name my_variable
```

### Listing database contents

```sh
thebigbam list-samples -d my_database.db
thebigbam list-contigs -d my_database.db
thebigbam list-mags -d my_database.db
thebigbam list-sample-metadata -d my_database.db
thebigbam list-contig-metadata -d my_database.db
thebigbam list-mag-metadata -d my_database.db
thebigbam list-variables -d my_database.db
```

`list-variables --detailed`` returns a description of each variable in addition to the list of variables.
