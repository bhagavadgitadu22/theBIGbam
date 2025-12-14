# MGFeatureViewer

Interactive visualization of mapping-derived features for genomes. Designed to work efficiently even with large metagenomic datasets.

MGFeatureViewer processes BAM alignment files to extract and visualize genomic features including:
- **Coverage tracks** (primary, secondary, supplementary reads)
- **Phage termini detection** (read starts/ends for identifying packaging sites)
- **Assembly quality metrics** (clippings, indels, mismatches, orientations, insert sizes)

Built with **Rust** for fast BAM processing and **Python/Bokeh** for interactive visualization.

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

See **[the installation guide](docs/INSTALL.md)** for more detailed instructions.

## Main usage

The MGFeatureViewer pipeline consists of two main steps: 
- Generation of an SQLite database summarizing your genomic features with Rust, using your mapping files (and optionally an annotation file and homemade csv files)
- Interactive visualization of the SQLite database content using Python and Bokeh

If you do not have an assembly annotation file (.genbank) or sorted mapping files with MD tags (.bam), you can generate them using the scripts provided in [the preprocessing section](docs/PREPROCESSING.md).

TODO: you can also skip assembly altogether by not providing a genbank file

### Database computation

This is the core of the MGFeatureViewer pipeline. It takes a genbank file containing the contigs of interest and a list of BAM files containing mappings against all—or a subset of—those contigs. 

If your mappings were performed in circular mode, do not forget to specify the --circular flag as described in the [the preprocessing section](docs/PREPROCESSING.md).

#### Features computed

MGFeatureViewer performs fast Rust-based computations on your BAM files to extract values corresponding to the features you request. Three types of features can currently be computed:
- **coverage**: computes per-position coverage for primary, secondary, and supplementary reads
- **phagetermini**: lanalyzes the start and end positions of mapped reads, particularly useful for identifying phage termini (see the [PhageTerm publication](https://www.nature.com/articles/s41598-017-07910-5)). A stringent coverage metric is calculated using only reads that begin (and for long reads, end) with a match. For each position, the number of read starts, read ends, and the ratio $τ = (reads\_starts + reads\_ends) / stringent\_coverage$ are reported
- **assemblycheck**: computes various metrics to assess contig completeness and contamination. This includes counts of clippings, indels, and mismatches. For long reads, read lengths are computed; for paired-end short reads, insert sizes and incorrect pair orientations are reported

The --modules (-m) option allows you to specify which of these features to compute. You can select one or more features by providing a comma-separated list (e.g., --modules coverage,phagetermini,assemblycheck).

#### Database compression

For each sample, Rust computations are performed for all contigs detected as present in the BAM file (i.e., those exceeding the --min_coverage threshold). After computation, the final database is further reduced through a compression step. Coverage metrics are compressed using run-length encoding (RLE): consecutive positions with similar values (within the --compress_ratio threshold) are stored as a single entry, preserving the overall signal while drastically reducing storage size. For phagetermini and assemblycheck metrics, values are compared to the local coverage and discarded if they fall below the --compress_ratio threshold, ensuring that only meaningful peaks are retained. When consecutive positions remain after filtering, RLE is applied to group them.

The output is an SQLite database containing all computed values. This database is typically ~1000× smaller than the original BAM files while retaining the essential characteristics of the mapping data.

#### Example with HK97 phage test data
```sh
mgfeatureviewer calculate -t 4 -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk --annotation_tool pharokka -b examples/inputs/HK97/bams -m coverage,phagetermini,assemblycheck -o examples/outputs/HK97/HK97.db --circular --compress-ratio=10 --min-coverage=50
```

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

## Additional utilities

### Preprocessing

Consult the [PREPROCESSING.md](docs/PREPROCESSING.md) for additional scripts to help with assembly annotation and read mapping.

### Database maintenance

Consult [DATABASE.md](docs/DATABASE.md) for instructions on how to read and modify your database after its initial creation.



## How to interpret the plots?

## On the coverage:
coverage and coverage_reduced use different kind of filtering.

coverage:
- only counts one cover if a read truly matches the given basepair
coverage_reduced:
- considers less reads than coverage because it only considers reads starting and finishing with a match (not clipping authorised at the ends of the read)
- counts one cover for all the positions between the first basepair and the last basepair of the read

Thus for a given position both coverage can be higher than the other one depending on the scenarii:
- coverage > coverage_reduced when part of the genome is missing in the assembly resulting in a lot of clippings at the position considered
- coverage_reduced > coverage when the basepair is missing in some of the organism' population

A particularly low coverage_reduced compared to coverage suggests your reads were not trimmed properly before mappings: you likely still have adapter sequences.

# TODO:
explain database structure
explain how to interrogate database with sql

mention how secondary is an approximation in case of doubled mapping because some primary at the ends of the contig would not have associated secondary. as i keep a minimum of 0 secondary, i will not see if some true secondary appear there. but this is minor and allows to go way faster than checking each read individually to remove secondary associated exactly to a primary read

lention run-pipeline command