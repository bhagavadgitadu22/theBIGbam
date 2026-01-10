# MGFeatureViewer

Interactive visualization of mapping-derived features for genomes. Designed to work efficiently even with large metagenomic datasets.

MGFeatureViewer processes BAM alignment files to extract and visualize genomic features including:

- **Coverage tracks** (primary, secondary, supplementary reads)
- **Assembly quality metrics** (clippings, indels, mismatches, orientations, insert sizes)
- **Termini detection** (positions of reads starts/ends, useful for identifying phage packaging sites)

Built with **Rust** for fast BAM processing and **Python/Bokeh** for interactive visualization.

---

## Quick Start

You need Rust and Python installed. Install Rust following instructions at: https://rust-lang.org/tools/install/

Then in command-line:

```bash
pip install git+https://github.com/bhagavadgitadu22/MGFeatureViewer

# Check main command works
mgfeatureviewer -h

# Run quick test with example data
mgfeatureviewer calculate \
  -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk \
  -b examples/inputs/HK97/ \
  -m coverage,phagetermini,assemblycheck \
  -o examples/outputs/HK97/test.db \
  --circular \
  -t 2

# Visualize interactively the test data
mgfeatureviewer serve --db examples/outputs/HK97/test.db --port 5006
# Open browser to http://localhost:5006
```

See [the installation guide](docs/INSTALL.md) for more detailed instructions.

---

## Main usage

The MGFeatureViewer pipeline consists of two main steps: 

- Generation of a DuckDB database summarizing your genomic features with Rust, using your input files
- Interactive visualization of the DuckDB database content using Python and Bokeh

### Database computation

This is the core of the MGFeatureViewer pipeline. It takes a list of BAM files containing read mappings against contigs of interest and extracts relevant features to store them in a lightweight DuckDB database. 

A genbank file containing annotations of your contigs of interest can also be provided to be added in the database.

#### Quick usage with the HK97 test data

Calculate command takes at least a directory of mapping files (-b) and an output path for the database (-o):

```sh
mgfeatureviewer calculate -b examples/inputs/HK97/bams --circular -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk --annotation_tool pharokka -o examples/outputs/HK97/HK97.db -t 4 
```

Here several optional parameters are added:

- --circular is necessary if the mapping files were generated using the --circular option of the tool

- -g option to provide an annotation file (here with genbank format .gbk)

- --annotation_tool to specify the bioinformatic tool used to generate the annotation file: it is used to color the genes on the genome track during the visualization

- -t option to specify the number of threads available to speed up computation

#### What input files do I need?

##### Mapping files:

**Parameter:** --bam_files DIRECTORY, short-version -b

Your mapping files need to be **sorted BAM files with MD tags**. If you have mapping files but not in the right format, Samtools is your friend! 

```shell
# to convert a SAM/BAM file in a sorted BAM file
samtools view -bS example.sam | samtools sort -o example.sorted.bam

# to add an index file to your BAM file
samtools index example.sorted.bam

# to add MD tags to your BAM file 
# you also need the fasta file used during the mapping step
samtools calmd -b example.sorted.bam ref.fasta > example.sorted.md.bam
```

Alternatively, if you do not have sorted mapping files with MD tags (.bam), you can generate them using the scripts provided in [the preprocessing section](docs/PREPROCESSING.md). The mapping scripts use a modified version of the standard mapper (minimap2) to allow for seamless circular mapping. 

**Warning:** If you use those mapping scripts with the --circular option, do not forget to specify the --circular flag as well when computing the database.

##### Annotation file:

**Parameter:** --genbank FILE, short-version -g

Annotation assembly file should be a .gbk or .gff file made with the tool of your choice: bakta for bacteria, pharokka or phold for phages, eggnog-mapper, etc.)

TODO: what files can be inputted exactly?

#### Which features can I calculate?

**Parameter:** --modules LIST, short version -m

MGFeatureViewer performs fast Rust-based computations on your BAM files to extract values corresponding to the modules you request. All modules are stored in the database unless you provide a specific list of modules.

TODO: how to specify modules

5 modules exist at the moment:

- **Coverage**: computes per-position coverage for primary, secondary, and supplementary reads
- **Mapping metrics per position:**
- **Long-read metrics:**
- **Paired-read metrics:**
- **Phage termini:**
- **assemblycheck**: computes various metrics to assess contig completeness and contamination. This includes counts of clippings, indels, and mismatches. For long reads, read lengths are computed; for paired-end short reads, insert sizes and incorrect pair orientations are reported
- **phagetermini**: analyzes the start and end positions of mapped reads, particularly useful for identifying phage termini (see the [PhageTerm publication](https://www.nature.com/articles/s41598-017-07910-5)). A stringent coverage metric is calculated using only reads that begin (and for long reads, end) with a match. For each position, the number of read starts, read ends, and the ratio $τ = (reads\_starts + reads\_ends) / stringent\_coverage$ are reported

The --modules (-m) option allows you to specify which of these features to compute. You can select one or more features by providing a comma-separated list (e.g., --modules coverage,phagetermini,assemblycheck).

If an annotation file is provided, the Genome module is additionnally computed. It keeps track of the contig annotations (positions of the coding sequences and their functions) and calculates the repeats contained within each contig using an autoblast.

A more detailed explanation of the modules and the features it contains is available in the preprocessing section](docs/FEATURES.md).

#### Database compression

For each sample, Rust computations are performed for all contigs detected as present in the BAM file (i.e., those exceeding the --min_coverage threshold). After computation, the final database is further reduced through a compression step. Coverage metrics are compressed using run-length encoding (RLE): consecutive positions with similar values (within the --compress_ratio threshold) are stored as a single entry, preserving the overall signal while drastically reducing storage size. For phagetermini and assemblycheck metrics, values are compared to the local coverage and discarded if they fall below the --compress_ratio threshold, ensuring that only meaningful peaks are retained. When consecutive positions remain after filtering, RLE is applied to group them.

The output is an SQLite database containing all computed values. This database is typically ~1000× smaller than the original BAM files while retaining the essential characteristics of the mapping data.

### Visualisation

Once you have computed your database, you can visualize it interactively using the `mgfeatureviewer serve` command. This launches a local web server that hosts the interactive plots.

Example command:

```sh
# the copy-paste of the file is necessary in my windows setup to avoid issues
cp examples/outputs/HK97/HK97.db ~/HK97.db
mgfeatureviewer serve --db ~/HK97.db --port 5006
```

The visualisation has 2 modes: per-sample and all-samples. "One sample" mode allows you to explore all computed features for a single sample, while all-samples mode enables comparison of a specific feature across multiple samples.

TODO: explain what we can do in this visualisation (filtering, zooming, exporting...)
TODO: also possibility to export from the server directly

Instead of exploring the plots interactively in your browser, you can also generate standalone HTML files containing the plots. You need to provide the database path, the contig of interest, the sample name or the feature to plot (for per-sample and all-samples plots, respectively).

Example commands:

```sh
mgfeatureviewer plot-per-sample -d examples/outputs/HK97/HK97.db -v "Coverage,Phage termini,Assembly check,test" --contig NC_002167.1 --sample HK97_R1_illumina_mapped_on_HK97_GCF_000848825.1 --html examples/outputs/HK97/HK97_illumina_per_sample.html

mgfeatureviewer plot-all-samples -d examples/outputs/HK97/HK97.db -v "Primary alignments" --contig NC_002167.1 --html examples/outputs/HK97/HK97_illumina_all_samples.html
```

TODO: check filtering.md that includes metrics explanation for what you could look for, interpretation.md

---

## Additional utilities

### Preprocessing

Consult the [PREPROCESSING.md](docs/PREPROCESSING.md) for additional scripts to help with assembly annotation and read mapping.

### Database maintenance

Consult [DATABASE.md](docs/DATABASE.md) for instructions on how to read and modify your database after its initial creation. 