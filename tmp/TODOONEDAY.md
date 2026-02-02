# TODO List

# Calculate:

## On computation:

* Threads param from calculate do not seem to work it is overriden by number threads i give with cluster

## On tests:

* Check RLE acts symmetrically if I turn a contig around: there should be no difference in compression
* Add plenty of additional tests during computation
* Add set of long reads where no read passes more than twice end of contig -> check then than coverage\_reduced = coverage and that no clippings

## On coverage features:

* Tools like megadepth, mosdepth, samtools depth, etc., do not treat a read as one continuous start–end interval, they break the read into aligned blocks, based on the CIGAR string -> check whether i do the same?
* Issues warnings if coverage is too low for reliable analysis for assembly check and phage termini?
* Identifies extended regions (>100 bp) where coverage drops below 2/3 of average. These regions may indicate deletions, gaps, or low-quality areas. Have metric representing percentage of whole contig that is low coverage?
* Problem because of duplication contig: supplementary, secondary and MAPQ do not represent reality anymore. How to fix that? Maybe a different-mapper would do a better job?

## On assembly check features:

* !important! Better index is number of bp unmatched divide by number of basepairs matched? -> average completeness
* !important! Add misalignment conclusions:
  For example: contig likely complete (no major clippings and insertions)
  contig likely misassembled (clippings, insertions, deletions, mismatches that represent more than 50%)
  -> status but also minimal estimation of bases missing and bases added and bases changed
  contig likely microdiverse (clippings, insertions, deletions, mismatches that do not represent full coverage)
  -> mention number of bases that are microdiverse, percentage of total and number of genes?
* Insertions should be saved at half positions

## On paired-end features:

* Missing mates and non-inward reads can take a stupid amount of data
* MAPQ has weird shapes maybe I should check it as well
* What kind of useful metric can I get from paired-reads features?
* Allow inner autocomplete, not only for words starting with words typed?

## On phage termini:

* I should understand where all reads start, sometimes weird that starts+clippings do not bring all answers
* !important! Problem when start\_reads one bp away from other start\_reads but on the other side and in a duplication -> becomes one duplication away and thus not merged
* I can use division by cov+ and cov- to get stronger signal for tau? Would imply different mathematical rules
* Normally reads covering a position RP = reads starting with a match RM + reads starting with a clipping RC -> I could a check like that telling how reads start of average
* Plutôt que de regarder zones de départ pourquoi ne pas regarder reads qui crossent tel/tel endroit -> I should be able to know how each read is starting and based on that i can have view to calculate percentage or reads passing each spot (number starting there divided by number continuing-ish)
* !important! Say that if too much coverage variation I cannot conclude safely on phage packaging: raise warning or stop?

# Utils:

## On custom variables:

* I already added blast parser for duplications, I can use to allow user to load blast results
* Also add inStrain results (NEW microdiversity module!)
* Add kmers: useful to debug PhageTermKmer?

## On other utils:

* Possibility to provide minimap2 parameters with tool? at least to say if map-ont pacbio...
* During mapping discard reads that pass twice the termini ending, cause they will create fake clippings in circular mode just because they are too long to map on the reference. Or discard those reads directly in calculate?
* Add possibility for user to annotate eukaryotic genome with eggnogg-mapper?
* Think about proper way to calculate sequencing\_type per sample?

## On database management:

* Also for add\_variable have optional --csv\_types to give all the parameters for the subplots
* Have add\_data as well for new contigs and samples: takes bams and gbk and adds all in db if not exists, produces a lot of warnings though to say custom\_variables will be missing and should be added manually...
* Also functions remove\_contig and remove\_samples

# On visualization:

* Non\_inward\_pairs, mate\_not\_mapped, mate\_on\_another\_contig are Curves in plots, can we see it when small range from from far away?
* !important! Be able to slide down in autocomplete lists -> pagination easy to implement?
* Flavobacterium147\_T9\_12\_phage2 -> no sites kept because many small site -> i could plot sites post termini classification via view? idea highlight on plot of phagetermini via point or similar?
* Show loading logo when waiting for plot to be generated
* Change title Bokeh Application in the web page

# On publication:

## On cluster:

* Sending one plot via ssh all right it is the second one which fails

## On release:

* Conda provides precompiled binaries, avoiding source builds entirely
* People should be able to deposit their dataset on the website made by evan?
* Way to restart a run when it stopped halfway -> need of checkpoints

## !important! On documentation:

* Beware what is what is written in sam: is mapping from read or contig point of view? refer to SAM documentation
* Specify that you need MD tag in mappings for phagetermini, add warning if MD not found
* Beware to keep circular option for both mapping and calculate!
* Say circular necessary for phagetermini because?

⦁	explain database structure

⦁	mention how secondary is an approximation in case of doubled mapping because some primary at the ends of the contig would not have associated secondary. as i keep a minimum of 0 secondary, i will not ⦁	see if some true secondary appear there. but this is minor and allows to go way faster than checking each read individually to remove secondary associated exactly to a primary read

⦁	mention run-pipeline command, test this command as well

