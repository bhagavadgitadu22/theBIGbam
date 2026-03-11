# How to run theBIGbam on big projects?

## Use case 1: medium dataset, 1 contig, 252 samples

The temporal dataset comes from a microbial evolution experiment in which *Escherichia coli* K-12 was grown on an agar megaplate with increasing concentrations of trimethoprim [(Baym et al. 2016)](https://www.science.org/doi/10.1126/science.aag0822).

252 samples were collected and sequenced with paired-end short-reads, allowing reconstruction of evolutionary trajectories with mutations and selections happening along the expanding microbial front.

So we have:

- 1 contig (*E. coli* K-12): we download the associated FASTA sequence and GenBank annotation file from NCBI (U00096.2)

- 252 samples: we can download them using the SRA ([prefetch and fasterq-dump commands](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump))

### Mapping

Then we first need to perform an alignment. The bacterial chromosome is complete (one circular contig) so we can use the circular mode of theBIGbam mapping:

```shell
cat samples.csv | while read -r read1 read2; do
    thebigbam mapping-per-sample -t 4 --circular -r1 ${read1} -r2 ${read2} -a Ecoli_K12_U00096.2.fasta -o ${read1/_1.fastq/.bam}
done
```

samples.csv contains one row per sample, with two columns specifying the paths to the paired FASTQ files

### Calculate

First we annotate the bacterial chromosome using [Bakta](https://github.com/oschwengers/bakta):

```shell
bakta --db <path_to_bakta_db> --output bakta_Ecoli_K12_U00096.2 --prefix Ecoli_K12_U00096.2 --complete --threads 10 Ecoli_K12_U00096.2.fasta
```

Then we provide the directory containing the alignment files we produced (bams_circular) and the bakta .gbff file to the calculate command:

```shell
thebigbam calculate -t 8 -b bams_circular -g bakta_Ecoli_K12_U00096.2/Ecoli_K12_U00096.2.gbff -o ecoli_megaplate.db
```

The command takes 12 minutes and 19 GB of RAM to compute with 8 threads. The final DuckDB database takes 872 MB to compare with the initial 35 GB of alignment files and 6 MB of annotations (40x compression factor).

### Exploration

We complement the database with extra metadata for each sample to know the ancestors and descendants of each sample on the megaplate:

```shell
thebigbam add-sample-metadata --db ecoli_megaplate.db --csv metadata_samples.csv
```

The metadata file looks like:

```csv
Sample,Instrument,Bases,Experiment,Isolate,Id,Parent,Root,Root_distance,Parent_nodes,Children_nodes,Terminal_nodes,MIC_color,MIC_level,MIC_approximative
Manifold222,Illumina HiSeq 2500,152393400,Multistep,Manifold222,222,217,187,3,"217;193;187","223;224;226","223;224;226",turquoise,5,100
Manifold223,Illumina HiSeq 2500,137961800,Multistep,Manifold223,223,222,187,4,"222;217;193;187",,223,blue,4,30
```

When serving the database to a localhost, that allows to pick any evolutionary path and see how the *E. coli* evolved throughout this path:

```shell
thebigbam serve --db ecoli_megaplate.db
```

On the figure below, we selected the evolutionary path leading to node 111 in the **Filtering** section of the **All Samples** mode. We picked the Mismatches variable as the plotted variable to see the mismatches appearing and disappearing along this path. To make it easily readable, we order samples by Root_distance in the **Plotting parameters**: 

![image](https://raw.githubusercontent.com/bhagavadgitadu22/theBIGbam/master/docs/images/VISUALIZATION.png)

On the figure, we zoomed in around position 1,237,440 to see the appearance of a particular mismatch. When zooming, we need to re-click APPLY to refresh the plots. If zooming close enough, as here, we see the sequence and codons appearing below the gene map on the top tracks.

---

## Use case 2: large dataset, 50k contigs, 192 samples

We tested the tool on an in-house large dataset of 192 glacier-fed stream metagenomes collected worldwide. We first assembled all metagenomes using MEGAHIT and identified all viral contigs > 10 kbp using geNomad. We combined the viral contigs found in all samples and dereplicated them using scripts from CheckV (at 95% identity and 75% coverage). We ended up with 50,043 viral contigs > 10 kbp (mostly viruses of bacteria a.k.a. bacteriophages).

So we have

- 50,043 contigs

- 192 samples

### Mapping

First we perform the mapping of each sample against all viral contigs. Complete phage genomes are circular during an infection cycle so it makes sense to use theBIGbam circular mode: 

```shell
cat samples.csv | while read -r read1 read2; do 
    sbatch -q serial -N 1 -n 1 -c 16 -t 10:00:00 --wrap="thebigbam mapping-per-sample -t 16 --circular -r1 ${read1} -r2 ${read2} -a viral_contigs_bigger_than_10kbp.fasta -o ${read1/_1.fastq/.bam}""
done
```

Considering the size of each sample, it is better practice to perform the mapping on a HPC cluster using slurm or an equivalent as done here.

### Calculate

First we annotate all viral contigs using [pharokka](https://github.com/gbouras13/pharokka):

```shell
sbatch -q serial -N 1 -n 1 -c 20 -t 20:00:00 --wrap="pharokka.py --meta -d <path_to_pharokka_db> -i viral_contigs_bigger_than_10kbp.fasta -o pharokka_viral_contigs_bigger_than_10kbp -t 20"
```

Then we provide the directory containing the alignment files we produced (bams_circular) and the pharokka .gbk file to the calculate command:

```shell
sbatch -q serial -N 1 -n 1 -c 32 -t 60:00:00 --wrap="thebigbam calculate -t 32 -b bams_circular -g pharokka_viral_contigs_bigger_than_10kbp/pharokka.gbk --annotation_tool pharokka -o NOMIS_VIRUS_corrected_26_2_2026_32threads_per_sample.db"
```

The command takes 13 hours and 61 GB of RAM to compute with 32 threads. The final DuckDB database takes 30 GB to compare with the initial 350 GB of alignment files and 3 GB of annotations (12x compression factor).

---

## Use case 3: biobank collection extended regularly

Not all projects compute the database once for a fixed set of contigs and samples. In long-term projects, new contigs and samples may be added regularly. The `--extend` flag of the `calculate` command allows an existing theBIGbam database to be updated with these additional data.

We use this function for an in-house phage collection called ALP. Whenever we sequence a new phage with Nanopore long-reads:

- We assemble it using Flye

- We perform the mapping using theBIGbam circular mode

- We annotate the contig using pharokka

And finally we update the ALP database with the new sample:

```shell
thebigbam calculate -t 8 -b new_phage.bam -g pharokka_new_phage/pharokka.gbk -o ALP.db --extend
```

When extending a database, previously processed samples will not contain mapping information for newly added contigs. If this information is required, those samples must be removed and reprocessed using updated BAM files that include the new contigs as references.

This is not useful for the collection presented here, but it may be relevant in **Use Case 2**. For example, if a new glacier-fed stream metagenome is sequenced, additional viral contigs may be assembled that were present but not assembled in previously sequenced metagenomes. Detecting these would require remapping all samples to the updated reference set. However, remapping all samples can be computationally expensive and should only be performed if it is expected to substantially affect the final results.
