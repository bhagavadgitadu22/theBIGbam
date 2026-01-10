## How to interpret the plots?

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
