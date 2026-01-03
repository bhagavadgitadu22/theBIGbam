https://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
MAPQ scores vary between different programs and you should not directly compare results from, say, Bowtie 2 and BWA.
You should look at your MAPQ scores and potentially filter out the really bad alignments.
Bioinformatics software documentation can often omit some really important details (see also my last blog post on this subject)

if linear mappings and no clippings on the sides -> linear contig
if circular mappings and no clippings on the sides and reads passing from one end to the other -> circular contig^^