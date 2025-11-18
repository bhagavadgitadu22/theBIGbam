# Test:
git clone https://github.com/bhagavadgitadu22/MGFeatureViewer
cd MGFeatureViewer
python -m pip install -e .

mgfeatureviewer -h
# or alternatively: python -m mgfeatureviewer.cli -h

mgfeatureviewer calculate -h
mgfeatureviewer calculate -t 4 -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk -a pharokka -b examples/inputs/HK97 -m coverage,phagetermini,assemblycheck -d examples/outputs/HK97/HK97.db

mgfeatureviewer add-variable examples/outputs/HK97/HK97.db test bars "#f60820" "Test Bar" examples/inputs/HK97/variable_test.csv
mgfeatureviewer add-variable examples/outputs/HK97/HK97.db test2 curve "#aef1c2" "Test Curve" examples/inputs/HK97/variable_test2.csv

mgfeatureviewer plot-per-sample -d examples/outputs/HK97/HK97.db -v "Coverage,Phage termini,Assembly check,test" --contig NC_002167.1 --sample HK97_GCF_000848825.1_with_MD --html examples/outputs/HK97/HK97_illumina_per_sample.html
mgfeatureviewer plot-all-samples -d examples/outputs/HK97/HK97.db -v "Coverage" --contig NC_002167.1 --html examples/outputs/HK97/HK97_illumina_all_samples.html

mgfeatureviewer serve --db examples/outputs/HK97/HK97.db --port 5006


# How to interpret the plots?

## On the coverage:
coverage and coverage_reduced use different kind of filtering.

coverage:
- only counts one cover if a read truly matches the given basepair
coverage_reduced:
- considers less reads than coverage because it only considers reads starting and finishing with a match (not clipping authorised at the ends of the read)
- counts one cover for all the positions between the first basepair and the last basepair of the read

Thus for a given position both coverage can be higher than the other one depending on the scenarii:
- coverage > coverage_reduced when part of the genome is missing in the assembly resulting in a lot of clippings at the position considered
- coverage_reduced > coverage when the basepair is missing in some of the organism' population

A particularly low coverage_reduced compared to coverage suggests your reads were not trimmed properly before mappings: you likely still have adapter sequences.
