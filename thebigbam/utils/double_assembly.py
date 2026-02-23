from pathlib import Path

from Bio import SeqIO

def _double_fasta(in_path: Path, out_path: Path) -> None:
    records = list(SeqIO.parse(str(in_path), "fasta"))
    with out_path.open("w") as fh:
        for rec in records:
            rec.seq = rec.seq + rec.seq
            SeqIO.write(rec, fh, "fasta")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Double fasta sequences for circular mapping")
    parser.add_argument('-a', '--assembly', required=True, help='Reference assembly to map against (fasta file)')
    parser.add_argument('-o', '--output', required=True, help='Output BAM path (will be written)')
    args = parser.parse_args()
    _double_fasta(Path(args.assembly), Path(args.output))
