First make sure reference fasta is UNIX compatible and without line breaks between sequences.

# Paired short reads 1 (PSR1):
Split fasta in 50 pieces and make pairs of reads of 100bp with an insert size of 200bp (no gap between mates):
```bash
awk '
BEGIN {
  READLEN = 100
  INSERT  = 200
  NPAIRS  = 50
  R1FILE = "50_read_pairs_for_test_10kbp_R1.fastq"
  R2FILE = "50_read_pairs_for_test_10kbp_R2.fastq"
}
# read single-contig FASTA - strip all whitespace
!/^>/ { 
  gsub(/[[:space:]]/, "")
  seq = seq $0 
}
# reverse complement function
function revcomp(s,    i, rc, c, comp) {
  comp["A"] = "T"; comp["T"] = "A"
  comp["G"] = "C"; comp["C"] = "G"
  comp["N"] = "N"
  rc = ""
  for (i = length(s); i >= 1; i--) {
    c = substr(s, i, 1)
    rc = rc comp[c]
  }
  return rc
}
END {
  L = length(seq)
  need = NPAIRS * INSERT
  if (L < need) {
    print "ERROR: contig too short (" L "), need " need > "/dev/stderr"
    exit 1
  }
  # pre-generate quality string (all "I" = Q40)
  qual = ""
  for (i = 1; i <= READLEN; i++) qual = qual "I"
  
  for (i = 0; i < NPAIRS; i++) {
    start = i * INSERT + 1
    frag  = substr(seq, start, INSERT)
    r1 = substr(frag, 1, READLEN)
    r2_fwd = substr(frag, READLEN + 1, READLEN)
    r2 = revcomp(r2_fwd)  # R2 is reverse complement
    id = "@pair_" (i+1)
    
    # R1 file
    print id "/1" > R1FILE
    print r1 > R1FILE
    print "+" > R1FILE
    print qual > R1FILE
    
    # R2 file
    print id "/2" > R2FILE
    print r2 > R2FILE
    print "+" > R2FILE
    print qual > R2FILE
  }
}' test_10000bp.fasta
```

# Paired short reads 2 (PSR2):
Same but inverting strand + and -
```bash
awk '
BEGIN {
  READLEN = 100
  INSERT  = 200
  NPAIRS  = 50
  R1FILE = "50_read_pairs_for_test_10kbp_inverted_R1.fastq"
  R2FILE = "50_read_pairs_for_test_10kbp_inverted_R2.fastq"
}
# read single-contig FASTA - strip all whitespace
!/^>/ {
  gsub(/[[:space:]]/, "")
  seq = seq $0
}
# reverse complement function
function revcomp(s,    i, rc, c, comp) {
  comp["A"] = "T"; comp["T"] = "A"
  comp["G"] = "C"; comp["C"] = "G"
  comp["N"] = "N"
  rc = ""
  for (i = length(s); i >= 1; i--) {
    c = substr(s, i, 1)
    rc = rc comp[c]
  }
  return rc
}
END {
  L = length(seq)
  need = NPAIRS * INSERT
  if (L < need) {
    print "ERROR: contig too short (" L "), need " need > "/dev/stderr"
    exit 1
  }

  # pre-generate quality string (Q40)
  qual = ""
  for (i = 1; i <= READLEN; i++) qual = qual "I"

  for (i = 0; i < NPAIRS; i++) {
    start = i * INSERT + 1
    frag  = substr(seq, start, INSERT)

    # original forward reads
    r1_fwd = substr(frag, 1, READLEN)
    r2_fwd = substr(frag, READLEN + 1, READLEN)

    # INVERT STRANDS
    r1 = revcomp(r1_fwd)   # R1 now reverse-complement
    r2 = r2_fwd            # R2 now forward

    id = "@pair_" (i+1)

    # R1 file
    print id "/1" > R1FILE
    print r1 > R1FILE
    print "+" > R1FILE
    print qual > R1FILE

    # R2 file
    print id "/2" > R2FILE
    print r2 > R2FILE
    print "+" > R2FILE
    print qual > R2FILE
  }
}' test_10000bp.fasta
```

# Paired short reads 3 (PSR3):
Take 100 times the contig with start circularly permuted and split it in 50 pairs of reads of 100bp:
```bash
awk '
BEGIN {
  srand(12345)      # reproducible
  READLEN = 100
  INSERT  = 200
  NITER   = 100
  NPAIRS  = 50
  R1FILE = "5000_read_pairs_for_test_10kbp_concatenated_100_times_R1.fastq"
  R2FILE = "5000_read_pairs_for_test_10kbp_concatenated_100_times_R2.fastq"
}
# read FASTA (single contig) - strip all whitespace
!/^>/ { 
  gsub(/[[:space:]]/, "")
  seq = seq $0 
}
# reverse complement function
function revcomp(s,    i, rc, c, comp) {
  comp["A"] = "T"; comp["T"] = "A"
  comp["G"] = "C"; comp["C"] = "G"
  comp["N"] = "N"
  rc = ""
  for (i = length(s); i >= 1; i--) {
    c = substr(s, i, 1)
    rc = rc comp[c]
  }
  return rc
}
END {
  L = length(seq)
  if (L < INSERT) {
    print "ERROR: contig shorter than insert size" > "/dev/stderr"
    exit 1
  }
  # pre-generate quality string (all "I" = Q40)
  qual = ""
  for (i = 1; i <= READLEN; i++) qual = qual "I"
  
  for (it = 1; it <= NITER; it++) {
    # 1) rotate contig by random amount (circular permutation)
    offset = int(rand() * L)
    perm = substr(seq, offset + 1) substr(seq, 1, offset)
    
    # 2) draw 50 paired reads SYSTEMATICALLY (tiling) from rotated contig
    for (p = 0; p < NPAIRS; p++) {
      start = (p * INSERT) + 1  # systematic tiling, not random
      frag  = substr(perm, start, INSERT)
      r1 = substr(frag, 1, READLEN)
      r2_fwd = substr(frag, READLEN + 1, READLEN)
      r2 = revcomp(r2_fwd)  # R2 is reverse complement
      id = "@iter" it "_pair" (p+1)
      
      # R1 file
      print id "/1" > R1FILE
      print r1 > R1FILE
      print "+" > R1FILE
      print qual > R1FILE
      
      # R2 file
      print id "/2" > R2FILE
      print r2 > R2FILE
      print "+" > R2FILE
      print qual > R2FILE
    }
  }
}' test_10000bp.fasta
```

# Long reads 1 (LR1):
Randomly break contig 100 times in 10 fragments to get fake fastq of 1000 long reads:
```bash
awk '
BEGIN{srand()}
!/^>/ {seq=seq$0}
END{
  for (r=1; r<=100; r++) {
    n=length(seq)
    for (i=1;i<=9;i++) cuts[i]=int(rand()*n)
    asort(cuts)
    cuts[0]=0; cuts[10]=n
    for (i=1;i<=10;i++) {
      frag=substr(seq,cuts[i-1]+1,cuts[i]-cuts[i-1])
      print "@read_"r"_"i
      print frag
      print "+"
      print gensub(/./,"I","g",frag)
    }
  }
}' test_10000bp.fasta > 1000_long_reads_for_test_10kbp.fastq
```

# Long reads 2 (LR2):
This time first contenate contig 100 times to itself and then break it in 1000 fragments:
```bash
awk '
BEGIN{srand()}

# read fasta, ignore header
!/^>/ {seq=seq$0}

END{
  # 1) concatenate contig 100 times
  big=""
  for(i=1;i<=100;i++) big=big seq
  L=length(big)

  # 2) generate 99 random cut points
  for(i=1;i<=99;i++) cuts[i]=int(rand()*L)+1

  # 3) sort cut points
  n=asort(cuts)

  # add boundaries
  cuts[0]=1
  cuts[100]=L+1

  # 4) emit 100 FASTQ reads
  for(i=1;i<=100;i++){
    len=cuts[i]-cuts[i-1]
    if(len<=0) continue
    frag=substr(big,cuts[i-1],len)

    print "@read_"i
    print frag
    print "+"
    gsub(/./,"I",frag)
    print frag
  }
}' test_10000bp.fasta > 100_long_reads_for_test_10kbp_concatenated_100_times.fastq
```

# Utils:
To check max length of generated reads:
```bash
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 50_read_pairs_for_test_10kbp_R1.fastq
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 50_read_pairs_for_test_10kbp_R2.fastq
# 100 OK
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 5000_read_pairs_for_test_10kbp_concatenated_100_times_R1.fastq
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 5000_read_pairs_for_test_10kbp_concatenated_100_times_R2.fastq
# 100 OK
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 1000_long_reads_for_test_10kbp.fastq
# 7420
awk 'NR%4==2 { if (length($0) > max) max = length($0) } END { print max }' 100_long_reads_for_test_10kbp_concatenated_100_times.fastq
# 44217
```

To see cigars:
```bash
samtools view tests/linear_bams/1000_long_reads.bam | awk '{print $1, $6}' | head -100
```

To get coverages:
```bash
samtools depth -a linear_bams/50_read_pairs_inverted.bam | cut -f 3 > coverages/coverage_50_read_pairs_inverted.csv
samtools depth -a linear_bams/50_read_pairs.bam | cut -f 3 > coverages/coverage_50_read_pairs.csv
samtools depth -a linear_bams/5000_read_pairs_concat.bam | cut -f 3 > coverages/coverage_5000_read_pairs_concat.csv
samtools depth -a linear_bams/100_long_reads_concat.bam | cut -f 3 > coverages/coverage_100_long_reads_concat.csv
samtools depth -a linear_bams/1000_long_reads.bam | cut -f 3 > coverages/coverage_1000_long_reads.csv

samtools depth -a circular_bams/1000_long_reads_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_1000_long_reads_circular.csv
samtools depth -a circular_bams/100_long_reads_concat_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_100_long_reads_concat_circular.csv
samtools depth -a circular_bams/5000_read_pairs_concat_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_5000_read_pairs_concat_circular.csv
samtools depth -a circular_bams/50_read_pairs_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_50_read_pairs_circular.csv
samtools depth -a circular_bams/50_read_pairs_inverted_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_50_read_pairs_inverted_circular.csv

samtools depth -a circular_bams/50_read_pairs_circular.bam | cut -f 3 | awk -v N=10000 '{ sum[(NR-1)%N] += $1 } END { for (i=0; i<N; i++) print sum[i]+0 }' > coverages/coverage_50_read_pairs_circular.csv
```