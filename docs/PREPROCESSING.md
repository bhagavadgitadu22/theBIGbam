# Preprocessing your data

## Mapping

You can start the pipeline with your reference contigs and read files. You may map a single sample or multiple samples at once. For each input file, the pipeline generates a sorted .bam file along with its corresponding index (.bam.bai).

To map multiple samples, provide a CSV file with four comma-separated columns:
$read1,read2,sequencing\_type,assembly\_file$

Only read1 is required.

- read2 is used only for paired-end reads
- sequencing_type and assembly_file may instead be supplied through the global options -s and -a if the same values apply to all samples

Example of csv file (available in examples/inputs/HK97/mapping_rows.csv):

```
examples/inputs/HK97/HK97_R1_illumina.fastq.gz,examples/inputs/HK97/HK97_R2_illumina.fastq.gz,short,examples/inputs/HK97/HK97_GCF_000848825.1.fasta
examples/inputs/HK97/HK97_nanopore.fastq.gz,,long,examples/inputs/HK97/HK97_GCF_000848825.1.fasta
```

The mapping strategy used here supports circularized mapping via the --circular flag. In this mode, each contig is concatenated to itself, effectively doubling the assembly size. This enables reads spanning the contig boundaries to map seamlessly without being split at the ends. If you use this option, ensure the --circular flag is also specified during database computation to correctly interpret the mappings.

Another feature of this mapping process is the computation of MD tags using samtools. MD tags are useful to quickly identify mismatches among the mapped reads.

Examples of commands to map a single sample or several:

```sh
DIR="examples/inputs/HK97"
mgfeatureviewer mapping-per-sample -s "short" -r1 "${DIR}/HK97_R1_illumina.fastq.gz" -r2 "${DIR}/HK97_R2_illumina.fastq.gz" -a "${DIR}/HK97_GCF_000848825.1.fasta" -o "${DIR}/HK97_GCF_000848825.1_with_MD.bam" --circular

mgfeatureviewer mapping-all-samples --csv "${DIR}/mapping_rows.csv" -o "${DIR}/bams" --circular
```

## Assembly annotation

MGFeatureViewer requires a genbank file containing all your contigs of interest. 

If you performed mappings for your different samples against several assemblies (typically one assembly per sample), you can concatenate all those assemblies into a single genbank file (example `cat assemblies/*.gbk > all_assemblies.gbk`). That can be useful if you want to annotate your contigs with several different software (typically bakta for bacteria and pharokka for phages).

To avoid potential pitfalls down the read because of the use of different annotation software or conflicting metadata, MGFeatureViewer provides a utility to re-annotate your assemblies in a consistent manner. You can either provide a fasta file containing all your contigs of interest with the --assembly option or the same csv file used for mapping with the --csv option (only the 4th column assembly will be used then to identify all the contig files).

For the moment, you can pick between 2 annotation softares using the --annotation_tool option: bakta for bacteria and pharokka for phages. If your assemblies contain 

Examples of commands to annotate your contigs of interest:

```sh
DIR="examples/inputs/HK97"
DB_DIR="/mnt/c/Users/boutroux/Documents/databases/pharokka_db"
mgfeatureviewer annotate-assemblies --threads 4 -a "${DIR}/HK97_GCF_000848825.1.fasta" --annotation_tool pharokka --annotation_db ${DB_DIR} -g "${DIR}/HK97.gbk"

# alternatively using a csv file
mgfeatureviewer annotate-assemblies  --csv "${DIR}/mapping_rows.csv" --annotation_tool pharokka --annotation_db ${DB_DIR} -g "${DIR}/HK97.gbk"
```

## Pipeline

You can start with a csv file listing your samples to map and annotate your assemblies in one go using the `mgfeatureviewer run-pipeline` command. This command will first map your reads, then annotate your assemblies, and finally compute the database and serve the result to your local browser.

Example commands:

```sh
DIR="examples/inputs/HK97"
DB_DIR="/mnt/c/Users/boutroux/Documents/databases/pharokka_db"

# Single-sample pipeline
mgfeatureviewer run-pipeline \
  --read1 "${DIR}/HK97_R1_illumina.fastq.gz" \
  --read2 "${DIR}/HK97_R2_illumina.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -s short \
  --annotation_tool pharokka \
  --annotation_db ${DB_DIR} \
  -m coverage \
  -o examples/outputs/HK97/pipeline_single_sample

# Multi-sample pipeline
mgfeatureviewer run-pipeline \
  --csv "${DIR}/mapping_rows.csv" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -s short \
  --annotation_tool pharokka \
  --annotation_db ${DB_DIR} \
  -m coverage,phagetermini,assemblycheck \
  --circular \
  -o examples/outputs/HK97/pipeline_multi_sample
```