sh /home/boutroux/scripts/pipeline_for_plot_start_reads_and_misassembly.sh 4 sr /home/boutroux/PhageTermKmer/assemblies_NCBI/HK97_GCF_000848825.1.fasta /home/boutroux/PhageTermKmer/pharokka/pharokka_HK97_GCF_000848825.1_no_dnaapler/pharokka.gbk /home/boutroux/PhageTermKmer/reads_illumina/HK97_R1_illumina.fastq.gz /home/boutroux/PhageTermKmer/reads_illumina/HK97_R2_illumina.fastq.gz
sh /home/boutroux/scripts/pipeline_for_plot_start_reads_and_misassembly.sh 4 ont /home/boutroux/PhageTermKmer/assemblies_NCBI/HK97_GCF_000848825.1.fasta /home/boutroux/PhageTermKmer/pharokka/pharokka_HK97_GCF_000848825.1_no_dnaapler/pharokka.gbk /home/boutroux/PhageTermKmer/reads_nanopore/HK97_nanopore.fastq

sh /home/boutroux/scripts/pipeline_for_plot_start_reads_and_misassembly.sh 4 ont /work/river/AKIRA/bacterial_genomes/assemblies/Sphingomonas_sp_5_plasmid.fasta /work/river/AKIRA/bacterial_genomes/bakta/Sphingomonas_sp_5_plasmid_bakta/Sphingomonas_sp_5_plasmid.gbff /work/river/AKIRA/bacterial_genomes/reads/Sphingomonas_sp_5.fastq.gz
sh /home/boutroux/scripts/pipeline_for_plot_start_reads_and_misassembly.sh 4 ont /work/river/AKIRA/bacterial_genomes/assemblies/Sphingomonas_sp_5_chromosome.fasta /work/river/AKIRA/bacterial_genomes/bakta/Sphingomonas_sp_5_bakta/Sphingomonas_sp_5.gbff /work/river/AKIRA/bacterial_genomes/reads/Sphingomonas_sp_5.fastq.gz

# To generate metagenomic test data
cd /scratch/boutroux/PhageTermKmer/test_with_virome_phageterm
# First needed to update the names to remove spaces in fasta headers
cat MOCK_VIROME_REF_PHAGE_FOR_MARTIN/Individual_Phages_Genomes_REF_AND_ASSEMBLY/*_assembly.fasta > virome_assembly.fasta
conda activate bioinfo
sbatch -N 1 -n 1 -c 10 -t 10:00:00 --wrap="sh /home/boutroux/scripts/circular_mapping_on_doubled_assembly.sh 10 short virome_assembly.fasta MOCK_VIROME_REF_PHAGE_FOR_MARTIN/VIROME_DATASET_AND_ASSEMBLIES_1_MILLIONS_READS_EACH_PHAGE/R1_1_MILLION_READS_FROM_EACH_PHAGE.fastq MOCK_VIROME_REF_PHAGE_FOR_MARTIN/VIROME_DATASET_AND_ASSEMBLIES_1_MILLIONS_READS_EACH_PHAGE/R2_1_MILLION_READS_FROM_EACH_PHAGE.fastq"
conda activate pharokka
sbatch -q serial -N 1 -n 1 -c 10 -t 10:00:00 --wrap="pharokka.py -t 10 --meta -d /work/river/Databases/pharokka_db/ -i virome_assembly.fasta -o pharokka_virome_assembly"

# number of reads in SAM file
samtools view -c examples/inputs/HK97_GCF_000848825.1_with_MD.bam
samtools view -c examples/inputs/Sphingomonas_sp_5_plasmid_with_MD.bam
# reads starting with soft-clipping
samtools view examples/inputs/HK97_GCF_000848825.1_with_MD.bam | awk '$6 ~ /^([0-9]+)S/' | wc -l
samtools view examples/inputs/Sphingomonas_sp_5_plasmid_with_MD.bam | awk '$6 ~ /^([0-9]+)S/' | wc -l
# reads ending with soft-clipping
samtools view examples/inputs/HK97_GCF_000848825.1_with_MD.bam | awk '$6 ~ /S$/' | wc -l
samtools view examples/inputs/Sphingomonas_sp_5_plasmid_with_MD.bam | awk '$6 ~ /S$/' | wc -l

# Plotting tests
# For HK97
python calculating_data.py -t 4 -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk -a pharokka -b examples/inputs/HK97 -m coverage,phagetermini,assemblycheck -d examples/outputs/HK97/HK97.db
python add_variable.py examples/outputs/HK97/HK97.db test bars "#f60820" "Test Bar" examples/inputs/HK97/variable_test.csv
python add_variable.py examples/outputs/HK97/HK97.db test2 curve "#aef1c2" "Test Curve" examples/inputs/HK97/variable_test2.csv

python plotting_data_per_sample.py -d examples/outputs/HK97/HK97.db -v "Coverage,Phage termini,Assembly check,test" --contig NC_002167.1 --sample HK97_GCF_000848825.1_with_MD --html examples/outputs/HK97/HK97_illumina_per_sample.html
python plotting_data_all_samples.py -d examples/outputs/HK97/HK97.db -v "Coverage" --contig NC_002167.1 --html examples/outputs/HK97/HK97_illumina_all_samples.html
python start_bokeh_server.py --db examples/outputs/HK97/HK97.db --port 5006

# For HK97 with mapping
# mgfeatureviewer mapping-per-sample --read1 examples/inputs/HK97/HK97_R1_illumina.fastq.gz --read2 examples/inputs/HK97/HK97_R2_illumina.fastq.gz --assembly examples/inputs/HK97/HK97_GCF_000848825.1.fasta --sequencing-type short --threads 4 --circular --output examples/outputs/HK97/test_mapping.bam
INPUT_DIR="examples/inputs/HK97"
OUTPUT_DIR="examples/outputs/HK97/test_pipeline"
mgfeatureviewer run-pipeline --threads 4 --csv ${INPUT_DIR}/mapping_rows.csv --circular --annotation-tool pharokka --annotation-db /mnt/c/Users/boutroux/Documents/databases/pharokka_db --modules coverage,phagetermini,assemblycheck --db ${OUTPUT_DIR}/HK97.db

mgfeatureviewer mapping-all-samples --threads 4 --csv ${INPUT_DIR}/mapping_rows.csv --circular --output-dir ${OUTPUT_DIR}
mgfeatureviewer annotate-assemblies --threads 4 --csv ${INPUT_DIR}/mapping_rows.csv --tool pharokka --db /mnt/c/Users/boutroux/Documents/databases/pharokka_db --output-genbank ${OUTPUT_DIR}/combined_annotations.gbk

mgfeatureviewer calculate -t 4 -g ${INPUT_DIR}/HK97_GCF_000848825.1_pharokka.gbk -a pharokka -b ${OUTPUT_DIR} -m coverage,phagetermini,assemblycheck -d ${OUTPUT_DIR}/HK97.db
mgfeatureviewer serve --db ${OUTPUT_DIR}/HK97.db --port 5006

# For Sphingomonas sp 5
python calculating_data.py -t 4 -g examples/inputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5.gbff -b examples/inputs/Sphingomonas_sp_5_chromosome -m coverage,phagetermini,assemblycheck -d examples/outputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5_chromosome.db
python start_bokeh_server.py --db examples/outputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5_chromosome.db --port 5006

# For virome
python calculating_data.py -t 6 -g examples/inputs/phageterm_virome/virome_assembly_pharokka.gbk -b examples/inputs/phageterm_virome -m coverage,phagetermini,assemblycheck -d examples/outputs/phageterm_virome/phageterm_virome.db
python start_bokeh_server.py --db examples/outputs/phageterm_virome/phageterm_virome.db --port 5006

# For phageterm references (long reads vs short reads and different library preps)
python calculating_data.py -t 4 -g examples/inputs/phageterm_long_short_and_nextera/pharokka.gbk -a pharokka -b examples/inputs/phageterm_long_short_and_nextera/bams_with_MD -m coverage,phagetermini,assemblycheck -d examples/outputs/phageterm_long_short_and_nextera/phageterm_long_short_and_nextera.db
python start_bokeh_server.py --db examples/outputs/phageterm_long_short_and_nextera/phageterm_long_short_and_nextera.db --port 5006

# For Sphingomonas sp 5 with mapping
# mgfeatureviewer mapping-per-sample --read1 examples/inputs/HK97/HK97_R1_illumina.fastq.gz --read2 examples/inputs/HK97/HK97_R2_illumina.fastq.gz --assembly examples/inputs/HK97/HK97_GCF_000848825.1.fasta --sequencing-type short --threads 4 --circular --output examples/outputs/HK97/test_mapping.bam
INPUT_DIR="examples/inputs/Sphingomonas_sp_5_chromosome"
OUTPUT_DIR="examples/outputs/Sphingomonas_sp_5_chromosome/test_mapping"
mgfeatureviewer mapping-all-samples --threads 4 --csv ${INPUT_DIR}/mapping_rows.csv --circular --output-dir ${OUTPUT_DIR}
mgfeatureviewer annotate-assemblies --threads 4 --csv ${INPUT_DIR}/mapping_rows.csv --tool bakta --db /mnt/c/Users/boutroux/Documents/databases/bakta_db/db-light --output-genbank ${OUTPUT_DIR}/combined_annotations.gbk

# For AKIRA database
mgfeatureviewer calculate -t 32 -g examples/inputs/AKIRA/pharokka.gbk -b examples/inputs/AKIRA/bams_filtered -m coverage,phagetermini,assemblycheck -o examples/outputs/AKIRA/db  # This is now a path to a folder instead of a SQLite filename
mgfeatureviewer serve --db examples/outputs/AKIRA/AKIRA/db --port 5006  # Use the new folder path
