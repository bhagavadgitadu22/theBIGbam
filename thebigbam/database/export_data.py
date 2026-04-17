"""Export a metric from the database as a contig × sample (or MAG × sample) TSV matrix."""

import csv
import sys

import duckdb

# Mapping: user-facing metric name → (contig_view, sql_column, mag_view_or_None)
# mag_view is None when the metric has no MAG-level equivalent.
METRIC_TO_VIEW = {
    # Explicit_coverage / Explicit_coverage_per_MAG
    "Aligned_fraction_percentage": ("Explicit_coverage", "Aligned_fraction_percentage", "Explicit_coverage_per_MAG"),
    "Above_expected_aligned_fraction": ("Explicit_coverage", "Above_expected_aligned_fraction", "Explicit_coverage_per_MAG"),
    "Read_count": ("Explicit_coverage", "Read_count", "Explicit_coverage_per_MAG"),
    "Coverage_mean": ("Explicit_coverage", "Coverage_mean", "Explicit_coverage_per_MAG"),
    "Coverage_median": ("Explicit_coverage", "Coverage_median", "Explicit_coverage_per_MAG"),
    "Coverage_trimmed_mean": ("Explicit_coverage", "Coverage_trimmed_mean", "Explicit_coverage_per_MAG"),
    "RPKM": ("Explicit_coverage", "RPKM", "Explicit_coverage_per_MAG"),
    "TPM": ("Explicit_coverage", "TPM", "Explicit_coverage_per_MAG"),
    "Coverage_coefficient_of_variation": ("Explicit_coverage", "Coverage_coefficient_of_variation", "Explicit_coverage_per_MAG"),
    "Relative_coverage_roughness": ("Explicit_coverage", "Relative_coverage_roughness", "Explicit_coverage_per_MAG"),
    # Explicit_misassembly / Explicit_misassembly_per_MAG
    "Misassembly_mismatches_per_100kbp": ("Explicit_misassembly", "Mismatches_per_100kbp", "Explicit_misassembly_per_MAG"),
    "Misassembly_deletions_per_100kbp": ("Explicit_misassembly", "Deletions_per_100kbp", "Explicit_misassembly_per_MAG"),
    "Misassembly_insertions_per_100kbp": ("Explicit_misassembly", "Insertions_per_100kbp", "Explicit_misassembly_per_MAG"),
    "Misassembly_clippings_per_100kbp": ("Explicit_misassembly", "Clippings_per_100kbp", "Explicit_misassembly_per_MAG"),
    "Collapse_bp": ("Explicit_misassembly", "Collapse_bp", "Explicit_misassembly_per_MAG"),
    "Collapse_per_100kbp": ("Explicit_misassembly", "Collapse_per_100kbp", "Explicit_misassembly_per_MAG"),
    "Expansion_bp": ("Explicit_misassembly", "Expansion_bp", "Explicit_misassembly_per_MAG"),
    "Expansion_per_100kbp": ("Explicit_misassembly", "Expansion_per_100kbp", "Explicit_misassembly_per_MAG"),
    # Explicit_microdiversity / Explicit_microdiversity_per_MAG
    "Microdiversity_mismatches_per_100kbp": ("Explicit_microdiversity", "Mismatches_per_100kbp", "Explicit_microdiversity_per_MAG"),
    "Microdiversity_deletions_per_100kbp": ("Explicit_microdiversity", "Deletions_per_100kbp", "Explicit_microdiversity_per_MAG"),
    "Microdiversity_insertions_per_100kbp": ("Explicit_microdiversity", "Insertions_per_100kbp", "Explicit_microdiversity_per_MAG"),
    "Microdiversity_clippings_per_100kbp": ("Explicit_microdiversity", "Clippings_per_100kbp", "Explicit_microdiversity_per_MAG"),
    "Microdiverse_bp_on_reads": ("Explicit_microdiversity", "Microdiverse_bp_on_reads", "Explicit_microdiversity_per_MAG"),
    "Microdiverse_bp_per_100kbp_on_reads": ("Explicit_microdiversity", "Microdiverse_bp_per_100kbp_on_reads", "Explicit_microdiversity_per_MAG"),
    "Microdiverse_bp_on_reference": ("Explicit_microdiversity", "Microdiverse_bp_on_reference", "Explicit_microdiversity_per_MAG"),
    "Microdiverse_bp_per_100kbp_on_reference": ("Explicit_microdiversity", "Microdiverse_bp_per_100kbp_on_reference", "Explicit_microdiversity_per_MAG"),
    # Explicit_side_misassembly (contig-only)
    "Coverage_first_position": ("Explicit_side_misassembly", "Coverage_first_position", None),
    "Contig_start_collapse_prevalence": ("Explicit_side_misassembly", "Contig_start_collapse_prevalence", None),
    "Contig_start_collapse_bp": ("Explicit_side_misassembly", "Contig_start_collapse_bp", None),
    "Contig_start_expansion_bp": ("Explicit_side_misassembly", "Contig_start_expansion_bp", None),
    "Coverage_last_position": ("Explicit_side_misassembly", "Coverage_last_position", None),
    "Contig_end_collapse_prevalence": ("Explicit_side_misassembly", "Contig_end_collapse_prevalence", None),
    "Contig_end_collapse_bp": ("Explicit_side_misassembly", "Contig_end_collapse_bp", None),
    "Contig_end_expansion_bp": ("Explicit_side_misassembly", "Contig_end_expansion_bp", None),
    "Contig_end_misjoint_mates": ("Explicit_side_misassembly", "Contig_end_misjoint_mates", None),
    "Normalized_contig_end_misjoint_mates": ("Explicit_side_misassembly", "Normalized_contig_end_misjoint_mates", None),
    # Explicit_topology (contig-only)
    "Circularising_reads": ("Explicit_topology", "Circularising_reads", None),
    "Circularising_reads_prevalence": ("Explicit_topology", "Circularising_reads_prevalence", None),
    "Circularising_inserts": ("Explicit_topology", "Circularising_inserts", None),
    "Circularising_insert_size_deviation": ("Explicit_topology", "Circularising_insert_size_deviation", None),
    "Normalized_circularising_inserts": ("Explicit_topology", "Normalized_circularising_inserts", None),
    # Explicit_phage_mechanisms (contig-only)
    "Packaging_mechanism": ("Explicit_phage_mechanisms", "Packaging_mechanism", None),
    "Repeat_length": ("Explicit_phage_mechanisms", "Repeat_length", None),
    "Terminase_distance": ("Explicit_phage_mechanisms", "Terminase_distance", None),
    "Terminase_percentage": ("Explicit_phage_mechanisms", "Terminase_percentage", None),
}

AVAILABLE_METRICS = sorted(METRIC_TO_VIEW.keys())

_METRIC_EPILOG = """\
Available metrics:

  coverage (MAG-mode only):
    Aligned_fraction_percentage, Above_expected_aligned_fraction,
    Read_count, Coverage_mean, Coverage_median, Coverage_trimmed_mean,
    RPKM, TPM, Coverage_coefficient_of_variation, Relative_coverage_roughness

  misassembly (MAG-mode only):
    Misassembly_mismatches_per_100kbp, Misassembly_deletions_per_100kbp,
    Misassembly_insertions_per_100kbp, Misassembly_clippings_per_100kbp,
    Collapse_bp, Collapse_per_100kbp, Expansion_bp, Expansion_per_100kbp

  microdiversity (MAG-mode only):
    Microdiversity_mismatches_per_100kbp, Microdiversity_deletions_per_100kbp,
    Microdiversity_insertions_per_100kbp, Microdiversity_clippings_per_100kbp,
    Microdiverse_bp_on_reads, Microdiverse_bp_per_100kbp_on_reads,
    Microdiverse_bp_on_reference, Microdiverse_bp_per_100kbp_on_reference

  side misassembly:
    Coverage_first_position, Contig_start_collapse_prevalence,
    Contig_start_collapse_bp, Contig_start_expansion_bp,
    Coverage_last_position, Contig_end_collapse_prevalence,
    Contig_end_collapse_bp, Contig_end_expansion_bp,
    Contig_end_misjoint_mates, Normalized_contig_end_misjoint_mates

  topology:
    Circularising_reads, Circularising_reads_prevalence,
    Circularising_inserts, Circularising_insert_size_deviation,
    Normalized_circularising_inserts

  phage mechanisms:
    Packaging_mechanism, Repeat_length,
    Terminase_distance, Terminase_percentage

  NOTE: not all metrics may be available depending on which modules
  were calculated. export reports an error for metrics that are not available.

examples:
  # Export per-contig coverage (default)
  thebigbam export -d my.db --metric Coverage_mean -o coverage.tsv

  # Export per-MAG coverage (requires a MAG-mode database)
  thebigbam export -d my.db --metric Coverage_mean --view mag -o mag_coverage.tsv
"""


def add_export_args(parser):
    """Register CLI arguments for the export subcommand."""
    import argparse

    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.epilog = _METRIC_EPILOG
    parser.add_argument('-d', '--db', required=True, help='Path to DuckDB database')
    parser.add_argument(
        '--metric', required=True, choices=AVAILABLE_METRICS,
        metavar='METRIC',
        help='Metric to export (see list below)',
    )
    parser.add_argument('-o', '--output', required=True, help='Output TSV file path')
    parser.add_argument(
        '--view', choices=['contig', 'mag'], default='contig',
        help="Aggregation level: 'contig' (default) exports per-contig values; "
             "'mag' exports per-MAG values (requires a MAG-mode database)",
    )


def run_export(args):
    """Export a metric as a contig × sample (or MAG × sample) TSV matrix. Returns 0 on success, 2 on error."""
    metric = args.metric
    contig_view, sql_column, mag_view = METRIC_TO_VIEW[metric]
    view_mode = getattr(args, 'view', 'contig')

    try:
        conn = duckdb.connect(args.db, read_only=True)
    except Exception as e:
        print(f"Error opening database: {e}", file=sys.stderr)
        return 2

    try:
        if view_mode == 'mag':
            if mag_view is None:
                print(
                    f"Error: metric '{metric}' is not available at the MAG level.",
                    file=sys.stderr,
                )
                return 2

            from thebigbam.database.database_getters import is_mag_mode
            if not is_mag_mode(conn):
                print(
                    "Error: --view mag requires a MAG-mode database. "
                    "This database was built in contig mode.",
                    file=sys.stderr,
                )
                return 2

            view_name = mag_view
            entity_col = "MAG_name"
            entity_label = "MAG"
        else:
            view_name = contig_view
            entity_col = "Contig_name"
            entity_label = "Contig"

        try:
            conn.execute(f"SELECT 1 FROM {view_name} LIMIT 0")
        except duckdb.Error:
            print(
                f"Error: view '{view_name}' not found in database. "
                f"Was the corresponding module calculated?",
                file=sys.stderr,
            )
            return 2

        rows = conn.execute(
            f"SELECT {entity_col}, Sample_name, {sql_column} "
            f"FROM {view_name} ORDER BY {entity_col}, Sample_name"
        ).fetchall()

        if not rows:
            print(f"No data found for metric '{metric}'.", file=sys.stderr)
            return 2

        entities_seen = {}
        samples_seen = {}
        matrix = {}
        for entity, sample, value in rows:
            entities_seen.setdefault(entity, None)
            samples_seen.setdefault(sample, None)
            matrix[(entity, sample)] = value

        entities = list(entities_seen.keys())
        samples = list(samples_seen.keys())

        with open(args.output, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow([entity_label] + samples)
            for entity in entities:
                row = [entity] + [
                    matrix.get((entity, s), '') for s in samples
                ]
                writer.writerow(row)

        print(f"Exported '{metric}' ({len(entities)} {entity_label.lower()}s x {len(samples)} samples) to {args.output}")
        return 0

    except Exception as e:
        print(f"Error exporting metric: {e}", file=sys.stderr)
        return 2
    finally:
        conn.close()
