from pathlib import Path
import subprocess
import shutil
import tempfile
import argparse
import csv
from typing import Optional

from Bio import SeqIO

def _double_fasta(in_path: Path, out_path: Path) -> None:
    records = list(SeqIO.parse(str(in_path), "fasta"))
    with out_path.open("w") as fh:
        for rec in records:
            rec.seq = rec.seq + rec.seq
            SeqIO.write(rec, fh, "fasta")

def add_mapping_per_sample_args(parser):
    parser.add_argument('-r1', '--read1', required=True, help='Read1 (fastq/fastq.gz)')
    parser.add_argument('-r2', '--read2', help='Read2 (fastq/fastq.gz; optional)')
    parser.add_argument('-a', '--assembly', required=True, help='Reference assembly to map against (fasta file)')
    parser.add_argument('--sequencing-type', required=True, choices=['long', 'paired-short', 'single-short'], help='Sequencing type: use "paired-short", "single-short" or "long"')
    parser.add_argument('--circular', action='store_true', help='Concatenate each contig to itself during the mapping to circularize it')
    parser.add_argument('-o', '--output', required=True, help='Output BAM path (will be written)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads to pass to minimap2 and samtools (default: 4)')

def run_mapping_per_sample(args):
    read2 = Path(args.read2) if getattr(args, 'read2', None) else None
    return map_with_mapper(
        args.threads,
        Path(args.assembly),
        args.sequencing_type,
        Path(args.read1),
        read2,
        Path(args.output),
        circular=bool(getattr(args, 'circular', False)),
    ) or 0

def map_with_mapper(threads: int, assembly_file: Path, sequencing_type: str, read1: Path,
                      read2: Optional[Path], output_file: Path, circular: bool = False) -> None:
    """Run mapper + samtools pipeline and produce final indexed BAM at `output_file`.

    This function expects `minimap2` and `samtools` to be on PATH.
    """
    for exe in ("minimap2", "samtools"):
        if shutil.which(exe) is None:
            raise FileNotFoundError(f"Required executable not found on PATH: {exe}")

    assembly = Path(assembly_file)
    if not assembly.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly}")

    work_assembly = assembly
    temp_files = []
    try:
        if circular:
            doubled = Path(tempfile.mkstemp(prefix=assembly.stem + "_doubled_", suffix=".fasta")[1])
            _double_fasta(assembly, doubled)
            work_assembly = doubled
            temp_files.append(doubled)

        # Build minimap2 command
        # from https://github.com/lh3/minimap2/blob/6d49eb690f3c32ae2b95a796951397bf598396f0/minimap2.1#L721
        # minimap2 options for sr 
        # -k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -r100 -p.5 -N20 -f1000,5000 -n2 -m25 -s40 -g100 -2K50m --heap-sort=yes --secondary=no
        # minimap2 options for map-ont
        # default mode -> no particular options needed
        if sequencing_type == "long":
            mapper_cmd = [
                "minimap2", "-ax", "map-ont", "-t", str(threads), str(work_assembly), str(read1)
            ]
        else:
            # just removed secondary=no to keep secondary alignments for short reads
            # -a flag is required to output SAM format (otherwise minimap2 outputs PAF)
            mapper_cmd = [
                "minimap2", "-a", "-k21", "-w11", "--sr", "--frag=yes", "-A2", "-B8", "-O12,32", "-E2,1", "-r100", "-p.5", "-N20", "-f1000,5000", "-n2", "-m25", "-s40", "-g100", "-2K50m", "--heap-sort=yes", "--secondary=yes", "-t", str(threads), str(work_assembly), str(read1)
            ]
            #mapper_cmd = [
            #    "minimap2", "-ax", "sr", "-t", str(threads), str(work_assembly), str(read1)
            #]
            #bwa_index_cmd = ["bwa-mem2", "index", str(work_assembly)]
            #subprocess.run(bwa_index_cmd, check=True)
            #mapper_cmd = ["bwa-mem2", "mem", "-t", str(threads), str(work_assembly), str(read1)]
            if read2:
                mapper_cmd.append(str(read2))

        sorted_bam = Path(tempfile.mkstemp(prefix=output_file.stem + "_sorted_", suffix=".bam")[1])
        #temp_files.append(sorted_bam)

        # Pipe: mapper | samtools view -bS -F 4 | samtools sort -o sorted_bam
        view_cmd = ["samtools", "view", "-@", str(threads), "-F", "4", "-bS", "-"]
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), "-"]
        print("COMMAND_MAP:", " ".join(mapper_cmd), flush=True)
        print("COMMAND_VIEW:", " ".join(view_cmd), flush=True)
        print("COMMAND_SORT:", " ".join(sort_cmd), flush=True)

        p1 = subprocess.Popen(mapper_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        # Close p1.stdout in parent to allow p1 to receive SIGPIPE if p2 exits
        p1.stdout.close()
        subprocess.run(sort_cmd, stdin=p2.stdout, check=True)
        p2.stdout.close()

        # Add MD tags and write final BAM (suppress warnings about secondary/supplementary alignments)
        with output_file.open("wb") as outfh:
            subprocess.run(
                ["samtools", "calmd", "-@", str(threads), "-b", str(sorted_bam), str(work_assembly)],
                stdout=outfh,
                stderr=subprocess.DEVNULL,
                check=True
            )

        # Index final BAM
        subprocess.run(["samtools", "index", "-@", str(threads), str(output_file)], check=True)

    finally:
        # cleanup temporary files
        for f in temp_files:
            try:
                Path(f).unlink()
            except Exception:
                pass

def add_mapping_all_args(parser):
    parser.add_argument('--csv', required=True, help='CSV file with comma-separated columns: read1,read2,sequencing_type,assembly_file')
    parser.add_argument('-a', '--assembly', help='Assembly file to use for all rows (overrides CSV field; optional)')
    parser.add_argument('--sequencing-type', choices=['long', 'paired-short', 'single-short'], help='Sequencing type to use for all rows (overrides CSV field; optional)')
    parser.add_argument('--circular', action='store_true', help='Concatenate each contig to itself during the mapping to circularize it')
    parser.add_argument('-o', '--output-dir', required=True, help='Directory to create and place outputs (must NOT exist)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads to pass to minimap2 and samtools (default: 4)')

def run_mapping_all(args):
    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    outdir = Path(args.output_dir)
    if outdir.exists():
        raise FileExistsError(f"Output directory already exists: {outdir}")
    outdir.mkdir(parents=True, exist_ok=False)

    global_seqtype = args.sequencing_type if getattr(args, 'sequencing_type', None) else None

    global_assembly = Path(args.assembly) if getattr(args, 'assembly', None) else None
    if global_assembly and not global_assembly.exists():
        raise FileNotFoundError(f"Provided assembly file not found: {global_assembly}")

    with csv_path.open(newline='') as fh:
        reader = csv.reader(fh)
        rows = list(reader)

    # If first row looks like header, skip it
    if rows:
        first = [c.lower() for c in rows[0]]
        if any("assembly" in c or "read" in c or "sequenc" in c for c in first):
            rows = rows[1:]

    if not rows:
        print("No rows to process in CSV")
        return 0

    for i, raw in enumerate(rows, start=1):
        # Expect assembly, read1, read2, sequencing_type
        row = [c.strip() for c in raw]
        while len(row) < 4:
            row.append("")
        read1, read2, seqtype, assembly = row[:4]
        read2 = read2 or None
        seqtype = seqtype.lower() if seqtype else None

        if not read1:
            raise ValueError(f"Row {i}: read1 is required")
        read1p = Path(read1)
        if not read1p.exists():
            raise FileNotFoundError(f"Row {i}: read1 file not found: {read1p}")

        read2p = Path(read2) if read2 else None
        if read2p and not read2p.exists():
            raise FileNotFoundError(f"Row {i}: read2 file not found: {read2p}")
        
        seqtype_to_use = global_seqtype if global_seqtype else seqtype
        if not seqtype_to_use:
            raise ValueError(f"Row {i}: sequencing_type is required")
        
        assembly_to_use = global_assembly if global_assembly else (Path(assembly) if assembly else None)
        if assembly_to_use is None:
            raise ValueError(f"Row {i}: no assembly provided (neither CSV nor --assembly)")
        if not assembly_to_use.exists():
            raise FileNotFoundError(f"Row {i}: assembly not found: {assembly_to_use}")

        path_read1 = Path(read1p)
        basename_read1 = path_read1.name.replace("".join(path_read1.suffixes), "")
        basename_assembly = Path(assembly_to_use).stem
        desired_bam = outdir / f"{basename_read1}_mapped_on_{basename_assembly}.bam"

        print(f"Processing row {i}: {read1p} -> assembly {assembly_to_use} (seqtype={seqtype_to_use}) -> {desired_bam}")

        # Create args namespace for run_mapping_per_sample
        sample_args = argparse.Namespace(
            read1=str(read1p),
            read2=str(read2p) if read2p else None,
            assembly=str(assembly_to_use),
            sequencing_type=seqtype_to_use,
            circular=getattr(args, 'circular', False),
            output=str(desired_bam),
            threads=args.threads,
        )
        run_mapping_per_sample(sample_args)

    print("All rows processed")
    return 0

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Map reads with minimap2 and produce MD-tagged indexed BAM")
    add_mapping_per_sample_args(parser)
    args = parser.parse_args()
    map_with_mapper(args.threads, Path(args.assembly), args.sequencing_type, Path(args.read1), Path(args.read2) if args.read2 else None, Path(args.output_file), args.circular)
