# TODO List

# Calculate:

## On computation:
- Threads param from calculate do not seem to work it is overriden by number threads i give with cluster

## On compressing:
- Check RLE acts symmetrically if I turn a contig around: there should be no difference in compression

## On tests:
- Add plenty of additional tests during computation
- add set of long reads where no read passes more than twice end of contig -> check then than coverage_reduced = coverage and that no clippings

## On coverage features:
- Tools like megadepth, mosdepth, samtools depth, etc., do not treat a read as one continuous start–end interval, they break the read into aligned blocks, based on the CIGAR string -> check whether i do the same?
- Issues warnings if coverage is too low for reliable analysis for assembly check and phage termini?
- Identifies extended regions (>100 bp) where coverage drops below 2/3 of average. These regions may indicate deletions, gaps, or low-quality areas. Have metric representing percentage of whole contig that is low coverage?
- Problem because of duplication contig: supplementary, secondary and MAPQ do not represent reality anymore. How to fix that? Maybe a different-mapper would do a better job?
- Is my coverage_variation index good?

## On assembly check features:
- Better index is number of bp unmatched divide by number of basepairs matched? -> average completeness
- Add mapping metrics per position conclusions:
For example: contig likely complete (no major clippings and insertions)
contig likely misassembled (clippings, insertions, deletions, mismatches that represent more than 50%)
-> status but also minimal estimation of bases missing and bases added and bases changed
contig likely microdiverse (clippings, insertions, deletions, mismatches that do not represent full coverage)
-> mention number of bases that are microdiverse, percentage of total and number of genes?
- Insertions should be saved ah half positions

## On paired-end features:
- Missing mates and non-inward reads can take a stupid amount of data
- MAPQ has weird shapes maybe I should check it as well
- What kind of useful metric can I get from paired-reads features?
- Have a module table with order between modules for visualization, otherwise long read metrics misplaces
- Allow inner autocomplete, not only for works starting with words typed?

## On phage termini:
- To avoid losing so many reads at the termini spot, find a way to keep reads if one end starts with a match (rather than both ends). Maybe a way is to split long reads in 2? or keep only first 1000 bp on both sides?
- Why so many clippings for long reads at the termini spot? would help to keep only first 1000bp of long reads?
- Maybe I should accept tiny clippings for long reads, instead of having exact_read_match -> relax a bit the threshold
- I should understand where all reads start, sometimes weird that starts+clippings do not bring all answers
- Problem when start_reads one bp away from other start_reads but on the other side and in a duplication -> becomes one duplication away and thus not merged
- Phage packaging determination affected by coverage_percentage maybe should not
- I can use division by cov+ and cov- to get stronger signal for tau? Would imply different mathematical rules
- Normally reads covering a position RP = reads starting with a match RM + reads starting with a clipping RC -> I could a check like that telling how reads start of average
- For long reads i lose SO MANY of the reads by forcing the read to start with a match: about 80% of the reads lost because of that for nanopore! Need to be more tolerant with reads i consider and allow some clippings?
- Plutôt que de regarder zones de départ pourquoi ne pas regarder reads qui crossent tel/tel endroit -> I should be able to know how each read is starting and based on that i can have view to calculate percentage or reads passing each spot (number starting there divided by number continuing-ish)

# Utils:

## On custom variables:
- I already added blast parser for duplications, I can use to allow user to load blast results
- Also add inStrain results (NEW microdiversity module!)
- Add possibility for user to supply table of additional info for contigs or for samples. Can be used to expand the filtering options if user did it correctly. Should be csv file that is provided with one row per Contig or Sample
- Add kmers: useful to debug PhageTermKmer?

## On other utils:
- Fix minimap2 problem, mapping-per-sample do not work at the moment
- Possibility to provide minimap2 parameters with tool? at leats to say if map-ont pacbio...
- Removed ax sr from mapping to get secondary reads check that is ok?
- --sequencing-type {long,paired-short,single-short} but message says Sequencing type: use "long" or "short" -> discrepancy
- During mapping discard reads that pass twice the termini ending, cause they will create fake clippings in circular mode just because they are too long to map on the reference. Or discard those reads directly in calculate?
- Add possibility for user to annotate eukaryotic genome with eggnogg-mapper?
- Think about proper way to calculate sequencing_type per sample?

## On database management:
- add and remove functions for variable, would take similar csv compared to mapping
- also for add_variable have optional --csv_types to give all the parameters for the subplots 
- also a add function for metadata to put sample and contig data in the database with csv like sample var1 var2 ... 
- have add_data as well for new contigs and samples: takes bams and gbk and adds all in db if not exists, produces a lot of warnings though to say custom_variables will be missing and should be added manually...
- also functions remove_contig and remove_samples

# On visualization:
- Small bug per-variable -> if i click directly on Contig field list of contigs is not actualised
- Sometimes autocomplete for Contig field (and others?) do not work at first use but work after clicking away and coming back
- Non_inward_pairs, mate_not_mapped, mate_on_another_contig are Curves in plots, can we see it when small range from from far away?
- If no genbank provided you should not show gene map option
- Filter to only show part of the contig possible (from bp x to bp y) -> could be very useful for eukaryotic contigs
- Filter contigs by genes present to look for a specific function for example
- Make 2 buttons instead of just Apply: "Peruse" and "Plot". Peruse lets you see one tables with as many columns as plots requested. Plus show general statistics associated to this contig in this sample at the top. Similarly adapted for All-samples view with table showing differents statistics per samples as well
- Possibility to extract tables with data as well

# On publication: 

## On cluster:
- Sending one plot via ssh all right it is the second one which fails

## On release:
- Conda provides precompiled binaries, avoiding source builds entirely
- People should be able to deposit their dataset on the website made by evan?
- Way to restart a run when it stopped halfway -> need of checkpoints

## On documentation:
- Beware what is what is written in sam: is mapping from read or contig point of view? refer to SAM documentation
- Specify that you need MD tag in mappings for phagetermini, add warming if MD not found
- Docu on metrics like coverage variation (how calculated?), completeness, contamination...
