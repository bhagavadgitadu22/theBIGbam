from pathlib import Path
import subprocess
import shutil
import tempfile
import argparse
import csv
from typing import Optional
import sys
from mgfeatureviewer import slurm_utils

from Bio import SeqIO

def _double_fasta(in_path: Path, out_path: Path) -> None:
    records = list(SeqIO.parse(str(in_path), "fasta"))
    with out_path.open("w") as fh:
        for rec in records:
            rec.seq = rec.seq + rec.seq
            SeqIO.write(rec, fh, "fasta")

def add_mapping_per_sample_args(parser):
    parser.add_argument('--read1', required=True, help='Read1 fastq/fastq.gz')
    parser.add_argument('--read2', help='Read2 fastq/fastq.gz (optional)')
    parser.add_argument('--assembly', required=True, help='Assembly fasta to map against')
    parser.add_argument('--sequencing-type', required=True, choices=['long', 'short'], help='Sequencing type: long or short')
    parser.add_argument('--threads', type=int, default=4, help='Threads to pass to minimap2/samtools')
    parser.add_argument('--circular', action='store_true', help='Treat assembly as circular (double sequences)')
    parser.add_argument('--output', required=True, help='Output BAM path (will be written)')

def run_mapping_per_sample(args):
    read2 = Path(args.read2) if getattr(args, 'read2', None) else None
    return map_with_minimap2(
        args.threads,
        Path(args.assembly),
        args.sequencing_type,
        Path(args.read1),
        read2,
        Path(args.output),
        circular=bool(getattr(args, 'circular', False)),
    ) or 0

def map_with_minimap2(threads: int, assembly_file: Path, sequencing_type: str, read1: Path,
                      read2: Optional[Path], output_file: Path, circular: bool = False) -> None:
    """Run minimap2 + samtools pipeline and produce final indexed BAM at `output_file`.

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
        if sequencing_type == "long":
            preset = "map-ont"
        else:
            preset = "sr"

        minimap2_cmd = [
            "minimap2", "-ax", preset, "--secondary=no", "-t", str(threads), str(work_assembly), str(read1)
        ]
        if read2:
            minimap2_cmd.append(str(read2))

        sorted_bam = Path(tempfile.mkstemp(prefix=output_file.stem + "_sorted_", suffix=".bam")[1])
        temp_files.append(sorted_bam)

        # Pipe: minimap2 | samtools view -bS -F 2052 | samtools sort -o sorted_bam
        view_cmd = ["samtools", "view", "-@", str(threads), "-F", "2052", "-bS", "-"]
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), "-"]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        # Close p1.stdout in parent to allow p1 to receive SIGPIPE if p2 exits
        p1.stdout.close()
        subprocess.run(sort_cmd, stdin=p2.stdout, check=True)
        p2.stdout.close()

        # Add MD tags and write final BAM
        with output_file.open("wb") as outfh:
            subprocess.run(["samtools", "calmd", "-@", str(threads), "-b", str(sorted_bam), str(work_assembly)], stdout=outfh, check=True)

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
    parser.add_argument('--csv', required=True, help='CSV file with columns: assembly_file,read1,read2,sequencing_type')
    parser.add_argument('--assembly', help='Assembly file to use for all rows (overrides CSV field)')
    parser.add_argument('--circular', action='store_true', help='Treat assemblies as circular (double sequences)')
    parser.add_argument('--output-dir', required=True, help='Directory to create and place outputs (must NOT exist)')
    parser.add_argument('--threads', type=int, default=4, help='Threads to pass to minimap2/samtools')
    # Slurm options
    parser.add_argument('--use-slurm', action='store_true', help='Dispatch mapping jobs via Slurm (sbatch)')
    parser.add_argument('--max-concurrent', type=int, default=20, help='Max concurrent Slurm array tasks')
    parser.add_argument('--threads-per-job', type=int, default=4, help='CPUs per Slurm job')
    parser.add_argument('--sbatch-mem', default='8G', help='Memory request for sbatch (e.g. 8G)')
    parser.add_argument('--sbatch-time', default='02:00:00', help='Time request for sbatch (e.g. 02:00:00)')

def run_mapping_all(args):
    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    outdir = Path(args.output_dir)
    if outdir.exists():
        raise FileExistsError(f"Output directory already exists: {outdir}")
    outdir.mkdir(parents=True, exist_ok=False)

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
        assembly, read1, read2, seqtype = row[:4]
        read2 = read2 or None
        seqtype = seqtype.lower() if seqtype else None

        if not read1:
            raise ValueError(f"Row {i}: read1 is required")

        assembly_to_use = global_assembly if global_assembly else (Path(assembly) if assembly else None)
        if assembly_to_use is None:
            raise ValueError(f"Row {i}: no assembly provided (neither CSV nor --assembly)")
        if not assembly_to_use.exists():
            raise FileNotFoundError(f"Row {i}: assembly not found: {assembly_to_use}")

        if not seqtype:
            raise ValueError(f"Row {i}: sequencing_type is required (long or short)")
        if seqtype not in {"long", "short"}:
            raise ValueError(f"Row {i}: sequencing_type must be 'long' or 'short', got: {seqtype}")

        read1p = Path(read1)
        if not read1p.exists():
            raise FileNotFoundError(f"Row {i}: read1 file not found: {read1p}")

        read2p = Path(read2) if read2 else None
        if read2p and not read2p.exists():
            raise FileNotFoundError(f"Row {i}: read2 file not found: {read2p}")

        basename_read1 = read1p.stem
        basename_assembly = Path(assembly_to_use).stem
        desired_bam = outdir / f"{basename_read1}_mapped_on_{basename_assembly}.bam"

        print(f"Processing row {i}: {read1p} -> assembly {assembly_to_use} (seqtype={seqtype}) -> {desired_bam}")

        if getattr(args, 'use_slurm', False):
            # Submit as a Slurm job array (one array covering all rows). We will build a single array submission
            # outside the loop — so do nothing per-row here. Fallthrough to array submission below.
            continue

        # Local execution
        map_with_minimap2(args.threads, assembly_to_use, seqtype, read1p, read2p, desired_bam, circular=getattr(args, 'circular', False))

    # If Slurm requested, submit an array job for all rows
    if getattr(args, 'use_slurm', False):
        if not slurm_utils.has_slurm():
            raise EnvironmentError("Slurm requested via --use-slurm but 'sbatch' not found on PATH")

        csv_path = Path(args.csv)
        outdir = Path(args.output_dir)
        array_size = len(rows)
        concurrency = int(getattr(args, 'max_concurrent', 20))
        cpus = int(getattr(args, 'threads_per_job', 4))
        mem = getattr(args, 'sbatch_mem', '8G')
        time_req = getattr(args, 'sbatch_time', '02:00:00')

        # module_cli will be evaluated inside the bash script and should reference shell variables
        module_cli = f"{sys.executable} -m mgfeatureviewer.cli mapping-per-sample --read1 \"$read1\" --read2 \"$read2\" --assembly \"$assembly\" --sequencing-type \"$seqtype\" --threads {args.threads} --output \"$output\""
        if getattr(args, 'circular', False):
            module_cli += ' --circular'

        jobid = slurm_utils.submit_array(csv_path, outdir, module_cli, array_size=array_size, concurrency=concurrency, cpus_per_task=cpus, mem=mem, time=time_req)
        print(f"Submitted Slurm array job {jobid} for {array_size} mapping tasks")
        return 0

    print("All rows processed")
    return 0

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Map reads with minimap2 and produce MD-tagged indexed BAM")
    add_mapping_per_sample_args(parser)
    args = parser.parse_args()
    map_with_minimap2(args.threads, Path(args.assembly), args.sequencing_type, Path(args.read1), Path(args.read2) if args.read2 else None, Path(args.output_file), args.circular)
