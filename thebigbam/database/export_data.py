"""Export a metric from the database as a contig × sample TSV matrix."""

import csv
import sys

import duckdb

# Mapping: user-facing metric name → (view_name, sql_column_name)
METRIC_TO_VIEW = {
    # Explicit_coverage
    "Aligned_fraction_percentage": ("Explicit_coverage", "Aligned_fraction_percentage"),
    "Above_expected_aligned_fraction": ("Explicit_coverage", "Above_expected_aligned_fraction"),
    "Read_count": ("Explicit_coverage", "Read_count"),
    "Coverage_mean": ("Explicit_coverage", "Coverage_mean"),
    "Coverage_median": ("Explicit_coverage", "Coverage_median"),
    "Coverage_trimmed_mean": ("Explicit_coverage", "Coverage_trimmed_mean"),
    "RPKM": ("Explicit_coverage", "RPKM"),
    "TPM": ("Explicit_coverage", "TPM"),
    "Coverage_coefficient_of_variation": ("Explicit_coverage", "Coverage_coefficient_of_variation"),
    "Coverage_relative_coverage_roughness": ("Explicit_coverage", "Coverage_relative_coverage_roughness"),
    # Explicit_misassembly
    "Misassembly_mismatches_per_100kbp": ("Explicit_misassembly", "Mismatches_per_100kbp"),
    "Misassembly_deletions_per_100kbp": ("Explicit_misassembly", "Deletions_per_100kbp"),
    "Misassembly_insertions_per_100kbp": ("Explicit_misassembly", "Insertions_per_100kbp"),
    "Misassembly_clippings_per_100kbp": ("Explicit_misassembly", "Clippings_per_100kbp"),
    "Collapse_bp": ("Explicit_misassembly", "Collapse_bp"),
    "Collapse_per_100kbp": ("Explicit_misassembly", "Collapse_per_100kbp"),
    "Expansion_bp": ("Explicit_misassembly", "Expansion_bp"),
    "Expansion_per_100kbp": ("Explicit_misassembly", "Expansion_per_100kbp"),
    # Explicit_microdiversity
    "Microdiversity_mismatches_per_100kbp": ("Explicit_microdiversity", "Mismatches_per_100kbp"),
    "Microdiversity_deletions_per_100kbp": ("Explicit_microdiversity", "Deletions_per_100kbp"),
    "Microdiversity_insertions_per_100kbp": ("Explicit_microdiversity", "Insertions_per_100kbp"),
    "Microdiversity_clippings_per_100kbp": ("Explicit_microdiversity", "Clippings_per_100kbp"),
    "Microdiverse_bp_on_reads": ("Explicit_microdiversity", "Microdiverse_bp_on_reads"),
    "Microdiverse_bp_per_100kbp_on_reads": ("Explicit_microdiversity", "Microdiverse_bp_per_100kbp_on_reads"),
    "Microdiverse_bp_on_reference": ("Explicit_microdiversity", "Microdiverse_bp_on_reference"),
    "Microdiverse_bp_per_100kbp_on_reference": ("Explicit_microdiversity", "Microdiverse_bp_per_100kbp_on_reference"),
    # Explicit_side_misassembly
    "Coverage_first_position": ("Explicit_side_misassembly", "Coverage_first_position"),
    "Contig_start_collapse_prevalence": ("Explicit_side_misassembly", "Contig_start_collapse_prevalence"),
    "Contig_start_collapse_bp": ("Explicit_side_misassembly", "Contig_start_collapse_bp"),
    "Contig_start_expansion_bp": ("Explicit_side_misassembly", "Contig_start_expansion_bp"),
    "Coverage_last_position": ("Explicit_side_misassembly", "Coverage_last_position"),
    "Contig_end_collapse_prevalence": ("Explicit_side_misassembly", "Contig_end_collapse_prevalence"),
    "Contig_end_collapse_bp": ("Explicit_side_misassembly", "Contig_end_collapse_bp"),
    "Contig_end_expansion_bp": ("Explicit_side_misassembly", "Contig_end_expansion_bp"),
    "Contig_end_misjoint_mates": ("Explicit_side_misassembly", "Contig_end_misjoint_mates"),
    "Normalized_contig_end_misjoint_mates": ("Explicit_side_misassembly", "Normalized_contig_end_misjoint_mates"),
    # Explicit_topology
    "Circularising_reads": ("Explicit_topology", "Circularising_reads"),
    "Circularising_reads_prevalence": ("Explicit_topology", "Circularising_reads_prevalence"),
    "Circularising_inserts": ("Explicit_topology", "Circularising_inserts"),
    "Circularising_insert_size_deviation": ("Explicit_topology", "Circularising_insert_size_deviation"),
    "Normalized_circularising_inserts": ("Explicit_topology", "Normalized_circularising_inserts"),
    # Explicit_phage_mechanisms
    "Packaging_mechanism": ("Explicit_phage_mechanisms", "Packaging_mechanism"),
    "Repeat_length": ("Explicit_phage_mechanisms", "Repeat_length"),
    "Terminase_distance": ("Explicit_phage_mechanisms", "Terminase_distance"),
    "Terminase_percentage": ("Explicit_phage_mechanisms", "Terminase_percentage"),
}

AVAILABLE_METRICS = sorted(METRIC_TO_VIEW.keys())

_METRIC_EPILOG = """\
available metrics (grouped by module):

  coverage:
    Aligned_fraction_percentage, Above_expected_aligned_fraction,
    Read_count, Coverage_mean, Coverage_median, Coverage_trimmed_mean,
    RPKM, TPM, Coverage_coefficient_of_variation, Coverage_relative_coverage_roughness

  misassembly:
    Misassembly_mismatches_per_100kbp, Misassembly_deletions_per_100kbp,
    Misassembly_insertions_per_100kbp, Misassembly_clippings_per_100kbp,
    Collapse_bp, Collapse_per_100kbp, Expansion_bp, Expansion_per_100kbp

  microdiversity:
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


def run_export(args):
    """Export a metric as a contig x sample TSV matrix. Returns 0 on success, 2 on error."""
    metric = args.metric
    view_name, sql_column = METRIC_TO_VIEW[metric]

    try:
        conn = duckdb.connect(args.db, read_only=True)
    except Exception as e:
        print(f"Error opening database: {e}", file=sys.stderr)
        return 2

    try:
        # Check the view exists
        try:
            conn.execute(f"SELECT 1 FROM {view_name} LIMIT 0")
        except Exception:
            print(
                f"Error: view '{view_name}' not found in database. "
                f"Was the corresponding module calculated?",
                file=sys.stderr,
            )
            return 2

        # Get all contigs and samples present in the view
        rows = conn.execute(
            f"SELECT Contig_name, Sample_name, {sql_column} "
            f"FROM {view_name} ORDER BY Contig_name, Sample_name"
        ).fetchall()

        if not rows:
            print(f"No data found for metric '{metric}'.", file=sys.stderr)
            return 2

        # Collect ordered unique contigs and samples
        contigs_seen = {}
        samples_seen = {}
        matrix = {}  # (contig, sample) → value
        for contig, sample, value in rows:
            contigs_seen.setdefault(contig, None)
            samples_seen.setdefault(sample, None)
            matrix[(contig, sample)] = value

        contigs = list(contigs_seen.keys())
        samples = list(samples_seen.keys())

        # Write TSV
        with open(args.output, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['Contig'] + samples)
            for contig in contigs:
                row = [contig] + [
                    matrix.get((contig, s), '') for s in samples
                ]
                writer.writerow(row)

        print(f"Exported '{metric}' ({len(contigs)} contigs x {len(samples)} samples) to {args.output}")
        return 0

    except Exception as e:
        print(f"Error exporting metric: {e}", file=sys.stderr)
        return 2
    finally:
        conn.close()
