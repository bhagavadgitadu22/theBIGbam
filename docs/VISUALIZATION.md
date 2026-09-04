# theBIGbam visualisation

## History

A right-hand panel, hidden by default, can be expanded from the web interface. It records each **APPLY FILTERS** and **APPLY** action performed during the session. You can restore the state corresponding to any recorded action at any time.

To save the complete history of actions performed during a session, click **SAVE SESSION**. This creates a JSON file in the directory from which the database is being served. You can later restore the session with: `thebigbam serve --db <db> --json <session>`. 

Alternatively, click **SAVE SETTINGS** to save the current settings to a JSON file in the directory from which the database is being served. Unlike **SAVE SESSION**, **SAVE SETTINGS** also records settings that have been modified but not yet applied. You can later restart with the same history using: `thebigbam serve --db <db> --json <settings>`.

## Hidden rules

When multiple isoforms share the same locus tag, only the **longest isoform** has its nucleotide sequence computed and is included in the CDS index for codon analysis. All isoforms are still displayed in the gene map, but codon/amino acid annotations are derived from the longest one only. This simplifies visualization by avoiding conflicting or redundant annotations at overlapping positions.

For MAG dot colors and gene map features, when multiple coloring rules apply to the same element, the first matching rule is applied. If a feature has multiple annotations (e.g., CRISPR and restriction–modification activity), all values are shown in labels (e.g., `CRISPR^RM` when using `activity`). However, color assignment (e.g., “Use random colors”) uses a single category per feature, even if multiple annotations are present.
