<div align="center">
  <img src="static/LOGO.png" alt="image" width="300" />
</div>

Interactive visualization of mapping-derived features for genomes. Designed to work efficiently even with large metagenomic datasets.

theBIGbam processes BAM alignment files to extract and visualize genomic features including:

- **Coverage tracks** (primary, secondary, supplementary reads)
- **Assembly quality features** (clippings, indels, mismatches, read lengths, paired-end properties)
- **Termini detection** (positions of reads starts/ends, useful for identifying phage packaging sites)

Built with **Rust** for fast BAM processing and **Python/Bokeh** for interactive visualization.

---

# Quick Start

You need conda installed. If you don't have it yet, you can install Miniconda following instructions at: https://docs.conda.io/en/latest/miniconda.html

Then in command-line:

```bash
conda env create -f thebigbam_env.yaml
conda activate thebigbam
pip install git+https://github.com/bhagavadgitadu22/theBIGbam
```

Be aware that it may take quite some time to compile the Rust code during installation.

To reduce the installation size, you can at the end remove the rust cache:

```bash
cargo clean
```

### Errors

To avoid potential compilation errors:

- On Linux: Install development headers: `sudo apt-get install libbz2-dev liblzma-dev zlib1g-dev clang libclang-dev`
- On macOS: `brew install xz bzip2 zlib`
- On HPC clusters: you may need to load the LLVM module first: `module load llvm`

If you cannot use any of these methods and still get errors related to `libclang` during installation such as :

```bash
thread 'main' (3889591) panicked at /home/gmichoud/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/bindgen-0.69.5/lib.rs:622:31:
        Unable to find libclang: "couldn't find any valid shared libraries matching: ['libclang.so', 'libclang-*.so', 'libclang.so.*', 'libclang-*.so.*'], set the `LIBCLANG_PATH` environment variable to a path where one of these files can be found (invalid: [])"
        
```

Please set the `LIBCLANG_PATH` environment variable to the path where your `libclang` library is located. For example:

```bash
export LIBCLANG_PATH=$(python -c "import os; print(os.environ['CONDA_PREFIX'] + '/lib')")
```

Then try installing again.

```bash
pip install git+https://github.com/bhagavadgitadu22/theBIGbam
```

## Check main command works

```bash
thebigbam -h
```

## Run quick test with example data

```bash
thebigbam calculate \
 -g tests/HK97/HK97_GCF_000848825.1_pharokka.gbk \
 -b tests/HK97/ \
 -o tests/HK97/test.db \
 -t 2
```

## Visualize interactively the test data

```bash
thebigbam serve --db tests/HK97/test.db --port 5006
```

If you're using wsl, you may need to install the package `wslu` to allow opening the browser from the wsl terminal:

```bash
sudo apt install wslu
```

Open browser to http://localhost:5006

See [the installation guide](docs/INSTALL.md) for more detailed instructions.

---

# Main usage

theBIGbam consists of two main steps: 

- Generation of a DuckDB database summarizing your genomic and mapping features with Rust, using your input files
- Interactive visualization of the DuckDB database content using Python and Bokeh

## Database computation

This is the core of theBIGbam. It takes a list of BAM files containing read mappings against contigs of interest and extracts relevant features to store them in a DuckDB database. Individual read information is discarded in favor of lightweight per-position averages.

A genbank file containing annotations of your contigs of interest can also be provided to be added in the database.

### Quick usage with the HK97 test data

Calculate command takes at least a directory of mapping files (-b) and an output path for the database (-o):

```sh
thebigbam calculate \
  -b examples/inputs/HK97/bams \
  -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk \
  -m "Coverage,Misalignment" \
  -o examples/outputs/HK97/HK97.db \
  -t 4
```

Here several optional parameters are added:

- -g option to provide an annotation file (GenBank `.gbk`/`.gbff` or GFF3 `.gff` format)

- -m option to select which modules to compute (comma-separated). If omitted, all modules are computed

- -t option to specify the number of threads available to speed up computation

- -s / --sequencing_type to force the sequencing type (`long`, `paired-short`, or `single-short`). If omitted, it is auto-detected per sample from BAM flags

**Note:** If your BAM files were generated using `thebigbam mapping-per-sample --circular`, the circular mapping information is automatically detected from the BAM headers during `calculate` — no additional flag is needed.

### What input files do I need?

You need at least one of the following: BAM mapping files (`-b`) or an annotation file (`-g`). If only an annotation file is provided, contig-level data (annotations, GC content, repeats) is populated without any sample-level mapping features.

#### Mapping files:

**Parameter:** --bam_files DIRECTORY, short-version -b

Your mapping files need to be **sorted BAM files with MD tags**. If you have mapping files but not in the right format, Samtools is your friend! 

```bash
# to convert a SAM/BAM file in a sorted BAM file
samtools view -bS example.sam | samtools sort -o example.sorted.bam

# to add an index file to your BAM file
samtools index example.sorted.bam

# to add MD tags to your BAM file 
# you also need the fasta file used during the mapping step
samtools calmd -b example.sorted.bam ref.fasta > example.sorted.md.bam
```

Alternatively, if you do not have sorted mapping files with MD tags (.bam extension), you can generate them using the scripts provided in [the preprocessing section](docs/PREPROCESSING.md). The mapping scripts use a modified version of the standard mapper (minimap2) to allow for seamless circular mapping. 

**Note:** If you use those mapping scripts with the `--circular` option, the circular mapping information is embedded in the BAM file headers and will be automatically detected during `calculate`.

#### Annotation file:

**Parameter (optional, strongly recommended):** --genbank FILE, short-version -g

Annotation file should be in GenBank (`.gbk`, `.gbff`, `.gb`) or GFF3 (`.gff`, `.gff3`) format, made with the tool of your choice: bakta for bacteria, pharokka or phold for phages, eggnog-mapper, etc.

### Which features can I calculate?

**Parameter (optional):** --modules COMMA-SEPARATED LIST, short version -m

theBIGbam performs fast Rust-based computations on your BAM files to extract relevant values and discard irrelevant information. All modules are computed and stored in the database unless you provide a specific subset of modules.

5 modules exist at the moment:

- **Coverage**: computes per-position coverage for primary, secondary, and supplementary reads, as well as the mapping quality (MAPQ)
- **Misalignment:** computes per-position number of clippings, insertions, deletions and mismatches
- **Long-reads:** computes per-position average length of reads
- **Paired-reads:** computes per-position average insert size of reads along with the number of incorrect pair orientations (non-inward pairs, mate unmapped or mapping or another contig)
- **Phage termini:** compute per-position coverage for primary-reads starting with an exact match. Among those reads, the number of mapped reads starting and ending is computed along with the tau ratio calculating the proportion of reads terminating at each position relative to the coverage

If an annotation file is provided, the **Genome** module is also computed. It keeps track of the contig annotations (positions of the coding sequences and their functions) and calculates the repeats contained within each contig using an autoblast.

A more detailed explanation of the modules and the features it contains is available in [the features section](docs/FEATURES.md).

### Database compression

**Parameters (optional):** --min_aligned_fraction, --min_coverage_depth, --variation_percentage, --contig_variation_percentage, --coverage_percentage,

Discarding the reads to only keep the main features of the mappings (like the coverage per position) already allows the DuckDB database to be way lighter than the original BAM file. The database itself is also structured to be as light as possible. 

First, the database is organised per contig per sample (qualified as a contig/sample pair thereafter). Only pairs relative to a contig present in a sample are stored in the database. The definition of a presence can be tweaked via two parameters: **--min_aligned_fraction** controls the minimum percentage of positions that received reads (default 50%, meaning a contig is considered present only if more than half of it received reads), and **--min_coverage_depth** sets the minimum mean coverage depth required for contig inclusion (default 0, i.e. disabled — set to e.g. 5 to filter out contigs with very low depth that produce noisy signals).

To further reduce the size of the database, values per feature are compressed rather than saving all positions. The type of compression depends on the type of plots:

- A **Run-Length Encoding approach (RLE)** is applied to the Curve plots (features from Coverage, Paired-reads and Long-reads module, "Coverage reduced" feature in Phage termini module). RLE stores consecutive genomic positions with similar values as a single entry, preserving the overall signal while substantially reducing storage size. The algorithm tracks the minimum and maximum values within each run, and a new entry is created when the range of values in the run exceeds a threshold relative to the smallest absolute value, defined as:
  
  $\text{max}(\text{run}) - \text{min}(\text{run}) > r \times \min(|\text{min}(\text{run})|, |\text{max}(\text{run})|)$
  
  where $r$ is the allowed variation ratio. This range-based criterion is symmetric (independent of position order) and prevents drift on gradual monotonic changes. The stored value for each run is the average of all values in the run. The allowed percentage of variation can be adjusted using the **--variation_percentage** parameter (default 50% ie 0.5) for mapping-related features, and the **--contig_variation_percentage** parameter (default: 10% ie 0.1) for contig-related features.
  
  Contig features use a separate parameter with a lower default value because only one value needs to be computed per contig and per position (O(n²)), whereas mapping features require computing one value per contig, per position, and per sample in which the contig is present (O(n³)).

- Only positions with values above a defined percentage of the local coverage are retained for Bar plots (Misalignment and Phage termini module except for "Coverage reduced" feature). For each position, values are compared to the local coverage and discarded if they fall below the **--coverage_percentage** threshold (default 10%), ensuring that only meaningful peaks are preserved.

The output is a DuckDB database that is typically ~100× smaller than the original BAM files, while retaining the essential characteristics of the mapping data.

## Visualization

Once the database has been computed, it can be visualized interactively using `thebigbam serve` command. This starts a local web server that hosts the interactive plots.

Example command:

```bash
thebigbam serve --db examples/outputs/HK97/HK97.db --port 5006
```

When accessing the web server (http://localhost:5006), you will be presented with a web interface:

<div align="center">
  <img src="static/VISUALIZATION.png" alt="image" width="800" />
</div>

### Selection panel

#### One Sample mode

You are initially in the **One Sample** mode, which allows exploration of all computed features for a single sample. Several sections on the left panel control what is plotted:

- **Filtering**: Only pairs of contig/samples matching the selected filters are available in the **Contigs** and **Samples** sections. For instance, if the contig length filter is set to >10 kbp, only contigs longer than this threshold will appear in the **Contigs** section, and only samples containing at least one such contig will appear in the **Samples** section. To consult the list of filters available have a look at [the filtering page](docs/FILTERS.md)

- **Contigs**: Select the contig you want to explore. If an annotation file was provided when creating the database, genomic features (gene maps, repeats, GC content, GC skew) can be selected for plotting by clicking on the contig features

- **Samples**: Select the sample you want to explore

- **Variables**: Select the features to plot. You can either use the checkboxes to select all features from a module or click individual features within a module

Finally, click **Apply** to visualize the requested features for the selected contig and sample. Alternatively, click **Peruse Data** to display tables containing the metrics and feature values.

#### All Samples mode

**All Samples** mode enables comparison of a specific feature across multiple samples. Compared to the **One Sample** mode, the **Samples** section is omitted, and only a single feature can be selected in the **Variables** section.

### Plotting

All plots leverage the full capabilities of Bokeh: you can pan, zoom, and hover over specific points to inspect local values (genomic position in base pairs and corresponding y-values).

Buttons in the top-right section allow you to disable pan, zoom, or hover interactions, reset the plots to their original state, and **export the current view as a PNG image**.

#### Adaptive resolution rendering

theBIGbam adapts the level of detail based on the viewing window size to ensure smooth, responsive plotting:

- **Full resolution (≤ 10 kb windows)**: All compressed RLE data points are expanded and plotted. Tooltips are enabled, showing exact position and value when hovering over data points. This provides maximum detail for focused exploration of specific genomic regions.

- **Downsampled view (> 10 kb windows)**: SQL-side binning reduces the number of points sent to the browser:
  - The visible window is divided into **1000 fixed-width bins**
  - Each RLE segment is assigned to a bin based on its **midpoint position**
  - Values within each bin are **averaged** to produce a single representative point per bin
  - Tooltips are disabled to improve performance
  - The x-coordinate for each bin is its **center position**

This two-tier approach balances detail and performance: small windows provide single-base resolution for detailed inspection, while large windows show trends across megabase-scale regions without overwhelming the browser.

#### Coverage normalization

For features that depend on local coverage depth (clippings, indels, mismatches, reads starts/ends, paired-read anomalies), an optional **"Plot relative to local coverage"** checkbox normalizes values to facilitate comparison across regions with varying coverage:

- **Enabled**: Y-axis shows the **ratio** of events to local coverage (value ÷ coverage), scaled between 0 and 1. For example, a value of 0.10 means 10% of reads at that position have the feature (e.g., 10 clippings per 100× coverage, or 50 clippings per 500× coverage).

- **Disabled** (default): Y-axis shows **absolute counts** (e.g., number of clippings, insertions, mismatches).

This normalization reveals whether anomalies are proportional to coverage (expected sequencing noise) or represent true biological signal or assembly errors that persist regardless of depth.

---

## Additional utilities

### Preprocessing

Consult the [PREPROCESSING.md](docs/PREPROCESSING.md) for additional scripts to help with assembly annotation and read mapping.

### Database maintenance

Consult [DATABASE.md](docs/DATABASE.md) for instructions on how to read and modify your database after its initial creation.

### Exporting data

Export any metric as a contig x sample TSV matrix:

```bash
thebigbam export -d my_database.db --metric Coverage_mean -o coverage.tsv
```

Run `thebigbam export -h` to see the full list of available metrics.
