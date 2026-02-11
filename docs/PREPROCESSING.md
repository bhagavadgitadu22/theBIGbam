# Preprocessing your data

## Mapping

You can map reads to a reference assembly using `mapping-per-sample`. The command generates a sorted .bam file along with its corresponding index (.bam.bai).

Available mapper choices: `minimap2-sr`, `bwa-mem2`, `minimap2-ont`, `minimap2-pb`, `minimap2-hifi`, `minimap2-no-preset`, `minimap2-sr-secondary`

The mapping strategy used here supports circularized mapping via the --circular flag. In this mode, each contig is concatenated to itself, effectively doubling the assembly size. This enables reads spanning the contig boundaries to map seamlessly without being split at the ends. If you use this option, ensure the --circular flag is also specified during database computation to correctly interpret the mappings.

Another feature of this mapping process is the computation of MD tags using samtools. MD tags are useful to quickly identify mismatches among the mapped reads.

You can pass extra parameters to the underlying mapper with `--minimap2-params` or `--bwa-params`. Extra parameters are appended after the preset flags, so if you specify a flag that already appears in the preset, your value overrides it (last-flag-wins behaviour in both minimap2 and bwa-mem2).

Examples:

```sh
DIR="examples/inputs/HK97"
thebigbam mapping-per-sample --mapper minimap2-sr \
  -r1 "${DIR}/HK97_R1_illumina.fastq.gz" \
  -r2 "${DIR}/HK97_R2_illumina.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_GCF_000848825.1_with_MD.bam" --circular

# Using interleaved paired-end FASTQ
thebigbam mapping-per-sample --mapper bwa-mem2 \
  --interleaved "${DIR}/HK97_interleaved.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_bwa.bam"

# Passing extra mapper parameters
thebigbam mapping-per-sample --mapper minimap2-sr \
  -r1 "${DIR}/HK97_R1_illumina.fastq.gz" \
  -r2 "${DIR}/HK97_R2_illumina.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_custom.bam" \
  --minimap2-params "--secondary=no -N5"

```
