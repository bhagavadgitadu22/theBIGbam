"""
Autoblast: BLAST each reference sequence against itself to identify duplicated regions.

This module performs self-BLAST for each reference sequence independently,
which is useful for detecting terminal repeats and other duplicated regions
relevant for phage termini analysis.
"""

from pathlib import Path
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _blast_single_contig(record: SeqRecord, evalue: float) -> List[str]:
    """Perform self-BLAST for a single contig.

    Args:
        record: BioPython SeqRecord to BLAST against itself
        evalue: E-value threshold for BLAST

    Returns:
        List of BLAST output lines (filtered to remove exact self-hits)
    """
    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write single sequence to temp file
        seq_fasta = Path(tmpdir) / "seq.fasta"
        SeqIO.write([record], str(seq_fasta), "fasta")

        # Run BLASTn of this sequence against itself
        blast_output = Path(tmpdir) / "blast_out.tsv"
        blastn_cmd = [
            "blastn",
            "-query", str(seq_fasta),
            "-subject", str(seq_fasta),
            "-evalue", str(evalue),
            "-outfmt", "6",
            "-out", str(blast_output)
        ]
        subprocess.run(blastn_cmd, check=True, capture_output=True)

        # Read results and filter out exact self-hits
        if blast_output.exists():
            with open(blast_output, 'r') as blast_handle:
                for line in blast_handle:
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        qstart, qend = int(fields[6]), int(fields[7])
                        sstart, send = int(fields[8]), int(fields[9])
                        # Skip if positions are identical (exact self-hit)
                        if qstart == sstart and qend == send:
                            continue
                    results.append(line)

    return results


def extract_fasta_from_genbank(genbank_path: Path, output_fasta: Path) -> bool:
    """Extract sequences from a GenBank file and write to FASTA format.

    Args:
        genbank_path: Path to GenBank file (.gbk, .gbff, .gb)
        output_fasta: Path to output FASTA file

    Returns:
        True if sequences were extracted, False if no sequences found
    """
    records = list(SeqIO.parse(str(genbank_path), "genbank"))
    if not records:
        return False

    # Check if sequences have actual content (not just annotations)
    records_with_seq = [r for r in records if len(r.seq) > 0 and str(r.seq) != "N" * len(r.seq)]
    if not records_with_seq:
        return False

    SeqIO.write(records_with_seq, str(output_fasta), "fasta")
    return True


def perform_autoblast(threads: int, assembly_path: Path, output_path: Path, evalue: float = 1e-10):
    """Perform self-BLAST for each reference sequence independently.

    Each sequence is BLASTed against itself only (not against other sequences),
    which identifies internal duplications like terminal repeats.

    Args:
        threads: Number of threads for BLAST (parallelizes across contigs)
        assembly_path: Path to FASTA file with reference sequences
        output_path: Path to output file with BLAST results (outfmt 6)
        evalue: E-value threshold for BLAST (default 1e-10)
    """
    # Parse all sequences from the assembly
    sequences = list(SeqIO.parse(str(assembly_path), "fasta"))

    if not sequences:
        raise ValueError(f"No sequences found in {assembly_path}")

    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run BLAST for each contig in parallel
    all_results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all BLAST jobs
        future_to_record = {
            executor.submit(_blast_single_contig, record, evalue): record.id
            for record in sequences
        }

        # Collect results as they complete
        for future in as_completed(future_to_record):
            contig_id = future_to_record[future]
            try:
                results = future.result()
                all_results.extend(results)
            except Exception as e:
                print(f"Warning: BLAST failed for {contig_id}: {e}")

    # Write all results to output file
    with open(output_path, 'w') as out_handle:
        for line in all_results:
            out_handle.write(line)

    return output_path


def add_autoblast_args(parser):
    """Add autoblast-specific arguments to an argument parser."""
    parser.add_argument('-a', '--assembly', required=True, help='Assembly file (FASTA format)')
    parser.add_argument('-o', '--output', required=True, help='Name of output file with BLASTn results')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Threads to pass to BLAST (default: 1)')
    parser.add_argument('-e', '--evalue', type=float, default=1e-10, help='E-value threshold (default: 1e-10)')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="BLASTn references against themselves to identify duplicated regions (e.g., terminal repeats)."
    )
    add_autoblast_args(parser)
    args = parser.parse_args()

    result = perform_autoblast(
        args.threads,
        Path(args.assembly),
        Path(args.output),
        evalue=args.evalue
    )
    print(f"Autoblast results written to: {result}")
