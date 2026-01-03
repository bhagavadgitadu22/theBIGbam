To test compressing ratios I checked number of points I get with different values for HK97:
DIR="examples/inputs/HK97"
mgfeatureviewer calculate -t 4 -b ${DIR}/HK97_GCF_000848825.1_with_MD.bam -m coverage,phagetermini,assemblycheck -o examples/outputs/HK97/pipeline_single_sample/features_no_gbk_curve10_bar100.db --circular --curve_ratio 0 --bar_ratio 0
mgfeatureviewer calculate -t 4 -b ${DIR}/HK97_GCF_000848825.1_with_MD.bam -m coverage,phagetermini,assemblycheck -o examples/outputs/HK97/pipeline_single_sample/features_no_gbk_curve10_bar100.db --circular --curve_ratio 10 --bar_ratio 10
mgfeatureviewer calculate -t 4 -b ${DIR}/HK97_GCF_000848825.1_with_MD.bam -m coverage,phagetermini,assemblycheck -o examples/outputs/HK97/pipeline_single_sample/features_no_gbk_curve10_bar100.db --circular --curve_ratio 10 --bar_ratio 100


Test with various inputs:
mgfeatureviewer calculate -t 4 -b examples/inputs/HK97 -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka -o examples/outputs/HK97/HK97.db
mgfeatureviewer serve --db examples/outputs/HK97/HK97.db --port 5006

mgfeatureviewer calculate -t 4 -b examples/inputs/Sphingomonas_sp_5_plasmid -g examples/inputs/Sphingomonas_sp_5_plasmid/Sphingomonas_sp_5_plasmid.gbff \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka -o examples/outputs/Sphingomonas_sp_5_plasmid/Sphingomonas_sp_5_plasmid.db
mgfeatureviewer serve --db examples/outputs/Sphingomonas_sp_5_plasmid/Sphingomonas_sp_5_plasmid.db --port 5006

mgfeatureviewer calculate -t 4 -b examples/inputs/Sphingomonas_sp_5_chromosome -g examples/inputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5.gbff \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka -o examples/outputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5_chromosome.db
mgfeatureviewer serve --db examples/outputs/Sphingomonas_sp_5_chromosome/Sphingomonas_sp_5_chromosome.db --port 5006

Test with several samples/contigss:
mgfeatureviewer calculate -t 4 -b examples/inputs/AKIRA/bams_filtered -g examples/inputs/AKIRA/pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka -o examples/outputs/AKIRA/akira.db
mgfeatureviewer serve --db examples/outputs/AKIRA/akira.db --port 5006

mgfeatureviewer calculate -t 4 -b examples/inputs/phageterm_virome -g examples/inputs/phageterm_virome/virome_assembly_pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka -o examples/outputs/phageterm_virome/phageterm_virome.db
mgfeatureviewer serve --db examples/outputs/phageterm_virome/phageterm_virome.db --port 5007

mgfeatureviewer calculate -t 4 -b examples/inputs/phageterm_long_short_and_nextera/bams_with_MD \
    -g examples/inputs/phageterm_long_short_and_nextera/pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --circular --annotation_tool pharokka \
    -o examples/outputs/phageterm_long_short_and_nextera/phageterm_long_short_and_nextera.db

mgfeatureviewer serve --db examples/outputs/phageterm_virome/phageterm_virome.db --port 5007


# Test on the cluster:
mgfeatureviewer installed in base conda environment
sbatch -q serial -N 1 -n 1 -c 32 -t 48:00:00 --wrap="mgfeatureviewer calculate \
    -t 32 -b /scratch/boutroux/test_MGFeatureViewer/bams_NOMIS \
    -g /work/river/NOMIS_VIRUS/15_final_viral_dataset/pharokka_viral_contigs_biggerthan10000_or_HQ_NOMIS/pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --annotation_tool pharokka \
    -o /scratch/boutroux/test_MGFeatureViewer/NOMIS_VIRUS.db"

Test with NOMIS first I made mapping files against doubled assembly using a sbatch array:
using list of samples and read files contained in indexes_used_for_KAUST_2025.xlsx
python double_assembly.py -a /work/river/NOMIS_VIRUS/15_final_viral_dataset/all_viruses_dereplicated_biggerthan10000_present_in_sediments.fna -o doubled_viruses_NOMIS.fasta
sbatch --array=1-192 -q serial -N 1 -n 1 -c 1 -t 48:00:00 --wrap="sh mgfeatureviewer_array.sbatch"

sbatch -q serial -N 1 -n 1 -c 32 -t 48:00:00 --wrap="mgfeatureviewer calculate \
    -t 32 -b /scratch/boutroux/test_MGFeatureViewer/NOMIS_mapping/ \
    -g /work/river/NOMIS_VIRUS/15_final_viral_dataset/pharokka_viral_contigs_biggerthan10000_or_HQ_NOMIS/pharokka.gbk \
    -m coverage,phagetermini,assemblycheck --annotation_tool pharokka --circular \
    -o /scratch/boutroux/test_MGFeatureViewer/NOMIS_VIRUS_duckdb.db"

# BEWARE: for LEA Genome I need eukaryotic annotations
# Use list_contig_filter_5kb_50X
python /scratch/boutroux/test_MGFeatureViewer/NOMIS_mapping/double_assembly.py -a contigs_hydrurus.fasta -o doubld_contigs_hydrurus.fasta

sbatch -q serial -N 1 -n 1 -c 20 -t 20:00:00 --wrap="minimap2 -ax ont -t 20 doubled_contigs_hydrurus.fasta /work/river/ICEBIO/Hannes/Hydrurus.fastq.gz | samtools view -@ 20 -bS -F 4 | samtools sort -@ 20 -o contigs_hydrurus.bam"
sbatch -q serial -N 1 -n 1 -c 12 -t 20:00:00 --wrap="emapper.py --cpu 12 -i contigs_hydrurus.fasta --itype genome --data_dir /work/river/Databases/eggnog-mapper-2.1.12 -o  contigs_hydrurus"

sbatch -q serial -N 1 -n 1 -c 8 -t 48:00:00 --wrap="mgfeatureviewer run-pipeline \
  -t 8 --read1 /work/river/ICEBIO/Hannes/Hydrurus.fastq.gz \
  -a /work/river/MGFeatureViewer/polished.assembly.fasta.gz \
  -m coverage,assemblycheck \
  -s short \
  --annotation_tool pharokka \
  --annotation_db ${DB_DIR} \
  -o examples/outputs/HK97/pipeline_single_sample
