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

You need Rust and Python installed. Install Rust following instructions at: https://rust-lang.org/tools/install/

Then in command-line:

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
 --circular \
 -t 2
```

## Visualize interactively the test data

```bash
thebigbam serve --db tests/HK97/test.db --port 5006
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
thebigbam calculate -b examples/inputs/HK97/bams --circular -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk --annotation_tool pharokka -m "Coverage","Misalignment"-o examples/outputs/HK97/HK97.db -t 4 
```

Here several optional parameters are added:

- --circular is necessary if the mapping files were generated using the --circular option of the tool

- -g option to provide an annotation file (GENBANK .gbk extension)

- --annotation_tool to specify the bioinformatic tool used to generate the annotation file: it is used to color the genes on the genome track during the visualization

- -t option to specify the number of threads available to speed up computation

### What input files do I need?

TODO: allow only genbank file without mapping files

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

**Warning:** If you use those mapping scripts with the --circular option, do not forget to specify the --circular flag as well when computing the database.

#### Annotation file:

**Parameter (optional, strongly recommended):** --genbank FILE, short-version -g

Annotation assembly file should be a GENBANK file (.gbk extension) made with the tool of your choice: bakta for bacteria, pharokka or phold for phages, eggnog-mapper, etc.

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

First, the database is organised per contig per sample (qualified as a contig/sample pair thereafter). Only pairs relative to a contig present in a sample are stored in the database. The definition of a presence can be tweaked via two parameters: **--min_aligned_fraction** controls the minimum percentage of positions that received reads (default 50%, meaning a contig is considered present only if more than half of it received reads), and **--min_coverage_depth** sets the minimum mean coverage depth required for contig inclusion (default 5×), filtering out contigs with very low depth that produce noisy signals.

To further reduce the size of the database, values per feature are compressed rather than saving all positions. The type of compression depends on the type of plots:

- A **Run-Length Encoding approach (RLE)** is applied to the Curve plots (features from Coverage, Paired-reads and Long-reads module, "Coverage reduced" feature in Phage termini module). RLE consists in storing consecutive positions with similar values as a single entry, preserving the overall signal while drastically reducing storage size. RLE stores consecutive genomic positions with similar values as a single entry, preserving the overall signal while substantially reducing storage size. By default, a new entry is created only when the variation between two consecutive values exceeds a threshold, defined as: 
  
  $∣x[i]−x[i−1]∣>0.5×min(x[i],x[i−1])$
  
  The allowed percentage of variation can be adjusted using the **--variation_percentage** parameter (default 50% ie 0.5) for mapping-related features, and the **--contig_variation_percentage** parameter (default: 10% ie 0.1) for contig-related features.
  
  Contig features use a separate parameter with a lower default value because only one value needs to be computed per contig and per position (O(n²)), whereas mapping features require computing one value per contig, per position, and per sample in which the contig is present (O(n³)).

- Only positions with values above a defined percentage of the local coverage are retained for Bar plots (Misalignment and Phage termini module except for "Coverage reduced" feature). For each position, values are compared to the local coverage and discarded if they fall below the **--compress_ratio** threshold (default 10% ie 0.1), ensuring that only meaningful peaks are preserved.

The output is a DuckDB database that is typically ~100× smaller than the original BAM files, while retaining the essential characteristics of the mapping data.

## Visualization

Once the database has been computed, it can be visualized interactively using `thebigbam serve` command. This starts a local web server that hosts the interactive plots.

Example command:

```bash
# the copy-paste of the file is necessary in my windows setup to avoid issues
cp examples/outputs/HK97/HK97.db ~/HK97.db
thebigbam serve --db ~/HK97.db --port 5006
```

When accessing the web server (http://localhost:5006), you will be presented with a web interface:

<div align="center">
  <img src="static/VISUALIZATION.png" alt="image" width="800" />
</div>

### Selection panel

#### One Sample mode

You are initially in the **One Sample** mode, which allows exploration of all computed features for a single sample. Several sections on the left panel control what is plotted:

- **Filtering**: Only pairs of contig/samples matching the selected filters are available in the **Contigs** and **Samples** sections. For instance, if the contig length filter is set to >10 kbp, only contigs longer than this threshold will appear in the **Contigs** section, and only samples containing at least one such contig will appear in the **Samples** section. To consult the list of filters available have a look at [the filtering page](docs/FILTERS.md)

- **Contigs**: Select the contig you want to explore. If a GenBank file was provided when creating the database, genomic features can be selected for plotting by clicking on the contig features. Currently, gene maps and repeats are supported

TODO: update list of features available and explain contigs part with figure

- **Samples**: Select the sample you want to explore

- **Variables**: Select the features to plot. You can either use the checkboxes to select all features from a module or click individual features within a module

Finally, click **Apply** to visualize the requested features for the selected contig and sample. Alternatively, click **Peruse Data** to display tables containing the metrics and feature values.

#### All Samples mode

**All Samples** mode enables comparison of a specific feature across multiple samples. Compared to the **One Sample** mode, the **Samples** section is omitted, and only a single feature can be selected in the **Variables** section.

TODO: add Contigs only mode as a simpler mode

### Plotting

All plots leverage the full capabilities of Bokeh: you can pan, zoom, and hover over specific points to inspect local values (genomic position in base pairs and corresponding y-values).

Buttons in the top-right section allow you to disable pan, zoom, or hover interactions, reset the plots to their original state, and **export the current view as a PNG image**.

TODO: add image with zoom on the tools, explain the download possibilities as well

---

## Additional utilities

### Preprocessing

Consult the [PREPROCESSING.md](docs/PREPROCESSING.md) for additional scripts to help with assembly annotation and read mapping.

### Database maintenance

Consult [DATABASE.md](docs/DATABASE.md) for instructions on how to read and modify your database after its initial creation. 