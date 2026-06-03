# Applying theBIGbam to real projects?

## Use case 1: medium dataset, 1 contig, 252 samples

The temporal dataset comes from a microbial evolution experiment in which *Escherichia coli* K-12 was grown on an agar megaplate with increasing concentrations of trimethoprim [(Baym et al. 2016)](https://www.science.org/doi/10.1126/science.aag0822).

252 samples were collected and sequenced with paired-end short-reads, allowing reconstruction of evolutionary trajectories with mutations and selections happening along the expanding microbial front.

So we have:

- **1 contig** (*E. coli* K-12): we download the associated FASTA sequence and GenBank annotation file from NCBI (U00096.2)

- **252 samples**: we can download them using the SRA ([prefetch and fasterq-dump commands](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump))

### Mapping

We first need to perform an alignment. The bacterial chromosome is complete (one circular contig) so we can use the circular mode of theBIGbam mapping:

```shell
mkdir bams_circular
cat samples.csv | while read -r read1 read2; do
    thebigbam mapping-per-sample -t 4 --circular -r1 ${read1} -r2 ${read2} -a Ecoli_K12_U00096.2.fasta -o bams_circular/${read1/_1.fastq/.bam}
done
```

samples.csv contains one row per sample, with two columns specifying the paths to the paired FASTQ files

### Calculate

To provide context to our read mapping, we annotate the bacterial chromosome using [Bakta](https://github.com/oschwengers/bakta):

```shell
bakta --db <path_to_bakta_db> --output bakta_Ecoli_K12_U00096.2 --prefix Ecoli_K12_U00096.2 --complete --threads 10 Ecoli_K12_U00096.2.fasta
```

Then we provide the directory containing the alignment files we produced (bams_circular) and the Bakta .gbff file to the calculate command:

```shell
thebigbam calculate -t 8 -b bams_circular -g bakta_Ecoli_K12_U00096.2/Ecoli_K12_U00096.2.gbff -o ecoli_megaplate.db
```

The command takes 10 minutes and 19 GB of RAM to compute with 8 threads. The final DuckDB database takes 2.2 GB to compare with the initial 25 GB of alignment files and 6 MB of annotations (11x compression factor).

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

On the figure, we zoomed in around position 1,237,440 to see the prevalence of a mismatch across generations. **When zooming, we need to re-click APPLY to refresh the plots**. If zooming close enough, as here, we see the sequence and codons appearing below the gene map on the top tracks.

---

## Use case 2: large dataset, 50k contigs, 192 samples

We tested the tool on an in-house large dataset of **192 glacier-fed stream metagenomes** collected worldwide. We first assembled all metagenomes using MEGAHIT and identified all viral contigs > 10 kbp using geNomad. We combined the viral contigs found in all samples and dereplicated them using scripts from CheckV (at 95% identity and 75% coverage). We ended up with **50,043 viral contigs > 10 kbp** (mostly viruses of bacteria a.k.a. bacteriophages).

### Mapping

We perform the mapping of each sample against all viral contigs. Complete phage genomes are circular during an infection cycle so it makes sense to use theBIGbam circular mode: 

```shell
mkdir bams_circular
cat samples.csv | while read -r read1 read2; do 
    sbatch -q serial -N 1 -n 1 -c 8 -t 10:00:00 --wrap="thebigbam mapping-per-sample -t 16 --circular -r1 ${read1} -r2 ${read2} -a viral_contigs_bigger_than_10kbp.fasta -o bams_circular/${read1/_1.fastq/.bam}"
done
```

Considering the size of each sample, it is better practice to perform the mapping on a HPC cluster using slurm or an equivalent as done here.

### Calculate

We annotate all viral contigs using [pharokka](https://github.com/gbouras13/pharokka):

```shell
sbatch -q serial -N 1 -n 1 -c 20 -t 20:00:00 --wrap="pharokka.py --meta -d <path_to_pharokka_db> -i viral_contigs_bigger_than_10kbp.fasta -o pharokka_viral_contigs_bigger_than_10kbp -t 20"
```

Then we provide the directory containing the alignment files we produced (bams_circular) and the pharokka .gbk file to the calculate command:

```shell
sbatch -q serial -N 1 -n 1 -c 16 -t 60:00:00 --wrap="thebigbam calculate -t 16 -b bams_circular -g pharokka_viral_contigs_bigger_than_10kbp/pharokka.gbk --annotation_tool pharokka -o thebigbam_virus.db"
```

The command takes 1 hour 28 min and 27 GB of RAM to compute with 16 threads. The final DuckDB database takes 36 GB to compare with the initial 350 GB of alignment files and 3 GB of annotations (10x compression factor).

---

## Use case 3: large dataset, 3000 MAGs, 192 samples

We used the tool on the same **192 samples** as in [Use case 2](#use-case-2-large-dataset-50k-contigs-192-samples), but focusing on the bacterial MAGs. Binning of the MEGAHIT assemblies was done with Binny, MetaBAT2, MaxBin2. Results from the 3 tools were integrated with DAS Tool and qualitative bins were filtered using CheckM2 (see [Michoud et al. 2025](https://www.nature.com/articles/s41564-024-01874-9) for more details). This resulted in a final dataset of **2,837 bacterial MAGs** for 1,582,129 contigs.

### Mapping

Mapping of each sample against all bacterial contigs was made with default minimap2 without circular mode as MAGs are fragmented and incomplete:

```shell
mkdir bams_linear
cat samples.tsv | while read -r sample r1 r2; do 
    sbatch -q serial -N 1 -n 1 -c 8 -t 10:00:00 --wrap="thebigbam mapping-per-sample -t 8 --circular -r1 $r1 -r2 $r2 -a finalBins -o bams_linear/${sample}.bam";
done
```

### Annotate

To make the most of thebigbam database, it is very useful to provide context by including genomic annotations among the input files. Those annotations can be a mix of different tools assembled together before running calculate.

Here we start by annotations the MAGs with Bakta:

```bash
mkdir bakta
for f in $(find finalBins -name "*.fa"); do
    name=$(basename $f | sed 's/.fa//');
    sbatch -q serial -N 1 -n 1 -c 4 -t 10:00:00 --wrap="bakta --db <path_to_bakta_db> --verbose --output bakta/${name}_bakta --prefix ${name} --meta --threads 4 $f";
done
```

We copy the 3,287 Bakta files in a new directory annotations_on_mags:

```bash
mkdir annotations_on_mags
find bakta -name "*.gff3" -exec cp {} annotation_on_mags \;
```

We also run defence-finder to identify defense and antidefense genes present in our MAGs. We first extract the list of proteins found by Bakta (ie by prodigal) and run [defense-finder](https://github.com/mdmparis/defense-finder) on this list:

```bash
find bakta -name "*.faa" -not -name "*.hypotheticals.faa" > mags_proteins.txt
sbatch --array=1-100 -q serial -N 1 -n 1 -c 1 -t 30:00:00 --wrap="sh defensefinder_array.sh mags_proteins.txt"
```

The script defensefinder_array.sh allows to use a job array to distribute jobs efficiently on the cluster:

```bash
#!/bin/bash
INPUT_LIST="$1"
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")
BASENAME=$(basename "$FILE" .faa)
OUTDIR="defenseFinder_${BASENAME}"
echo "Running on: $FILE"
echo "Output: $OUTDIR"

defense-finder run \
    -w 1 \
    --antidefensefinder \
    -o "$OUTDIR" \
    "$FILE"
```

We can then update our Bakta annotation files with the newly identified defense proteins. To do so, we use the `add-contig-annotations` utility from `thebigbam`. This tool takes a directory containing initial annotation files (one per MAG, here in .gff3 format) and a CSV file listing the annotations to be added. The `--match-by` option specifies the columns used to identify features to extend in the annotation files. Here, we match features using the `feature_type` column (`CDS`) and the `locus_tag` column, which uniquely identifies each CDS in Bakta.

First, we need to reformat the DefenseFinder output to match the input format expected by `add-contig-annotations`:

```bash
echo "replicon,locus_tag,gene_name,hit_pos,model_fqn,sys_id,sys_loci,locus_num,sys_wholeness,sys_score,sys_occ,hit_gene_ref,hit_status,hit_seq_len,hit_i_eval,hit_score,hit_profile_cov,hit_seq_cov,hit_begin_match,hit_end_match,counterpart,used_in,type,subtype,activity" > MAGs_defense_finder_genes_reformatted.csv
# replace "," by "|" and replace "\t" by ","
sed 's/,/|/g' MAGs_defense_finder_genes.tsv | sed $'s/\t/,/g' >> MAGs_defense_finder_genes_reformatted.csv
# adding column feature_type with value "CDS" for all rows
awk 'BEGIN{FS=OFS=","} NR==1{print $0,"feature_type"; next} {print $0,"CDS"}' MAGs_defense_finder_genes_reformatted.csv > MAGs_defense_finder_genes_reformatted_CDS_only.csv
```

Then we use `add-contig-annotations` to get a new directory of annotation files comprising both the Bakta and the DefenseFinder qualifiers for each coding sequence:

```bash
sbatch -q serial -N 1 -n 1 -c 1 -t 10:00:00 --wrap="thebigbam add-contig-annotations -g annotations_on_mags --csv MAGs_defense_finder_genes_reformatted_CDS_only.csv --match-by feature_type,locus_tag -o annotations_on_mags_enriched"
```

### Calculate

We provide the directory containing the alignment files we produced (bams_linear) and the enriched annotation files file to the calculate command:

```shell
sbatch -q serial -N 1 -n 1 -c 16 -t 30:00:00 --wrap="thebigbam calculate -t 16 -g annotations_on_mags_enriched -b mappings_on_mags --view mag -o thebigbam_bacteria_mag_view.db --min_aligned_fraction 10 --min_occurrences 5"
```

We used specific options for this use case: 

- `--view mag` allows to calculate metrics per MAG in addition to per contig. Information like coverage per MAG, etc. are recorded in the database

- `--min_aligned_fraction 10` to consider a MAG present in a sample when more than 10% of its nucleotides received at least one mapping. MAGs are way bigger than the viral genomes from [Use case 2](#use-case-2-large-dataset-50k-contigs-192-samples), so a lower aligned fraction typically is enough to consider something present

- `--min_occurrences 5` to only keep events occuring more than 5 times. This means a mismatch will be recorded at a position if more than 5 reads harbor this mismatch there

Despite mapping files 10 times bigger (~ 3 Tb) than in [Use case 2](#use-case-2-large-dataset-50k-contigs-192-samples), the command still runs quite fast: it takes 23 hours and 37 GB of RAM to compute with 16 threads. The final DuckDB database takes 473 GB to compare with the initial 2.9 TB of alignment files and 23 GB of annotations (6.2x compression factor).

### Exploration

The database obtained is quite massive so we keep it stored remotely on the cluster. To visualize it we serve it on the cluster:

```bash
sbatch -q serial -N 1 -n 1 -c 8 -t 10:00:00 --wrap="thebigbam serve --db /work/river/NOMIS_VIRUS/05_MAGs/thebigbam_bacteria_mag_view.db"
```

Such a massive database needs about 5 minutes for preloading before the webpage is ready, after which the log should print: "Server ready".

We need to find the node where the job is running (here jst234) and relay it locally via ssh. For us, `ssh -N -L 5006:jst234:5006 user@remote.server.com` fails because the cluster firewall blocks direct access to login nodes. If you meet a similar problem you can try to jump via the login node with:

```bash
ssh -J user@remote.server.com -N -L 5006:localhost:5006 user@jst234
```

We can now open `localhost:5006` in a web browser. For MAG datasets, theBIGbam can visualize complete MAGs instead of individual contigs. Here, the coverage track of a MAG is displayed, with the MAG track above it indicating the boundaries of the constituent contigs:

![image](https://raw.githubusercontent.com/bhagavadgitadu22/theBIGbam/master/docs/images/VISUALIZATION_MAG.png)

---

## Use case 4: biobank collection extended regularly

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

This is not useful for the collection presented here, but it may be relevant in [Use case 2](#use-case-2-large-dataset-50k-contigs-192-samples). For example, if a new glacier-fed stream metagenome is sequenced, additional viral contigs may be assembled that were present but not assembled in previously sequenced metagenomes. Detecting these would require remapping all samples to the updated reference set. However, remapping all samples can be computationally expensive and should only be performed if the final results are expected to substantially differ.
