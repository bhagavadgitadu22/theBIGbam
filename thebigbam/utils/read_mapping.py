from pathlib import Path
import subprocess
import shutil
import shlex
import sys
import tempfile
import warnings
from typing import Optional

from Bio import SeqIO

from thebigbam.utils.convert_circular_bam import convert_circular_bam

MAPPER_CHOICES = [
    'minimap2-sr', 'bwa-mem2',
    'minimap2-ont', 'minimap2-pb', 'minimap2-hifi',
    'minimap2-no-preset', 'minimap2-sr-secondary',
]

def _double_fasta(in_path: Path, out_path: Path) -> None:
    records = list(SeqIO.parse(str(in_path), "fasta"))
    with out_path.open("w") as fh:
        for rec in records:
            rec.seq = rec.seq + rec.seq
            SeqIO.write(rec, fh, "fasta")

def _validate_read_inputs(args):
    """Validate read-input constraints that argparse cannot express."""
    if getattr(args, 'read2', None) and not getattr(args, 'read1', None):
        raise ValueError("-r2/--read2 requires -r1/--read1")

def _validate_extra_params(args):
    """Warn if mapper-specific extra params don't match the chosen mapper."""
    mapper = args.mapper
    if getattr(args, 'bwa_params', None) and mapper.startswith('minimap2'):
        warnings.warn(f"--bwa-params ignored: mapper is '{mapper}' (not bwa-mem2)")
    if getattr(args, 'minimap2_params', None) and mapper == 'bwa-mem2':
        warnings.warn(f"--minimap2-params ignored: mapper is '{mapper}' (not minimap2)")

def add_mapping_per_sample_args(parser):
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads to pass to mapper and samtools (default: 4)')
    reads_group = parser.add_mutually_exclusive_group(required=True)
    reads_group.add_argument('-r1', '--read1', help='Read1 (fastq/fastq.gz)')
    reads_group.add_argument('--interleaved', help='Interleaved paired-end reads (fastq/fastq.gz; alternative to -r1/-r2)')
    parser.add_argument('-r2', '--read2', help='Read2 for paired-end reads (fastq/fastq.gz)')
    parser.add_argument('-a', '--assembly', required=True, help='Reference assembly to map against (fasta file)')
    parser.add_argument('-o', '--output', required=True, help='Output BAM path (will be written)')
    parser.add_argument('--mapper', choices=MAPPER_CHOICES, default='minimap2-sr-secondary', help='Mapper and preset to use. Default: short-read preset minimap2-sr-secondary that uses the same options as -ax sr but keeping secondary reads')
    parser.add_argument('--circular', action='store_true', help='Concatenate each contig to itself during the mapping to circularize it')
    parser.add_argument('--keep-unmapped', action='store_true', help='Keep unmapped reads in the output BAM (default: discard them)')
    parser.add_argument('--minimap2-params', dest='minimap2_params', default=None, help='Extra parameters for minimap2 (e.g., --minimap2-params "--secondary=no -N5")')
    parser.add_argument('--bwa-params', dest='bwa_params', default=None, help='Extra parameters for bwa-mem2 (e.g., --bwa-params "-M -B6")')

def run_mapping_per_sample(args):
    _validate_read_inputs(args)
    _validate_extra_params(args)

    read1 = Path(args.read1) if getattr(args, 'read1', None) else None
    read2 = Path(args.read2) if getattr(args, 'read2', None) else None
    interleaved = Path(args.interleaved) if getattr(args, 'interleaved', None) else None
    minimap2_params = getattr(args, 'minimap2_params', None)
    bwa_params = getattr(args, 'bwa_params', None)

    return map_with_mapper(
        args.threads,
        Path(args.assembly),
        args.mapper,
        read1,
        read2,
        Path(args.output),
        circular=bool(getattr(args, 'circular', False)),
        keep_unmapped=bool(getattr(args, 'keep_unmapped', False)),
        interleaved=interleaved,
        minimap2_params=minimap2_params,
        bwa_params=bwa_params,
    ) or 0

def _get_version() -> str:
    """Get theBIGbam version from pyproject.toml or fallback."""
    try:
        from importlib.metadata import version, PackageNotFoundError
        return version("thebigbam")
    except PackageNotFoundError:
        return "unknown"


def _inject_bam_headers(bam_path: Path, circular: bool, threads: int, command_line: str) -> None:
    """Inject @PG and @CO headers into a BAM file using samtools reheader."""
    version = _get_version()

    # Read existing header
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        capture_output=True, text=True, check=True,
    )
    header_lines = result.stdout.rstrip("\n").split("\n")

    # Append @PG and @CO lines
    pg_line = f"@PG\tID:theBIGbam\tPN:theBIGbam\tVN:{version}\tCL:{command_line}"
    circular_str = "true" if circular else "false"
    co_line = f"@CO\ttheBIGbam:circular={circular_str}"

    header_lines.append(pg_line)
    header_lines.append(co_line)
    new_header = "\n".join(header_lines) + "\n"

    # Write new header to temp file, reheader to temp BAM, then replace original
    header_file = Path(tempfile.mkstemp(suffix=".sam", prefix="header_")[1])
    reheadered_bam = Path(tempfile.mkstemp(suffix=".bam", prefix="reheader_")[1])
    try:
        header_file.write_text(new_header)
        with reheadered_bam.open("wb") as out:
            subprocess.run(
                ["samtools", "reheader", str(header_file), str(bam_path)],
                stdout=out, check=True,
            )
        shutil.move(str(reheadered_bam), str(bam_path))
    finally:
        for f in (header_file, reheadered_bam):
            try:
                f.unlink()
            except OSError:
                pass


def map_with_mapper(threads: int, assembly_file: Path, mapper: str, read1: Optional[Path],
                      read2: Optional[Path], output_file: Path, circular: bool = False,
                      keep_unmapped: bool = False, interleaved: Optional[Path] = None,
                      minimap2_params: Optional[str] = None, bwa_params: Optional[str] = None) -> None:
    """Run mapper + samtools pipeline and produce final indexed BAM at `output_file`."""

    # Executable check based on mapper
    if mapper.startswith('minimap2'):
        for exe in ("minimap2", "samtools"):
            if shutil.which(exe) is None:
                raise FileNotFoundError(f"Required executable not found on PATH: {exe}")
    elif mapper == 'bwa-mem2':
        for exe in ("bwa-mem2", "samtools"):
            if shutil.which(exe) is None:
                raise FileNotFoundError(f"Required executable not found on PATH: {exe}")

    assembly = Path(assembly_file)
    if not assembly.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly}")

    work_assembly = assembly
    temp_files = []
    bwa_index_sidecars = []
    try:
        if circular:
            doubled = Path(tempfile.mkstemp(prefix=assembly.stem + "_doubled_", suffix=".fasta")[1])
            _double_fasta(assembly, doubled)
            work_assembly = doubled
            temp_files.append(doubled)

        # bwa-mem2 indexing
        if mapper == 'bwa-mem2':
            bwa_index_cmd = ["bwa-mem2", "index", str(work_assembly)]
            print("COMMAND_INDEX:", " ".join(bwa_index_cmd), flush=True)
            subprocess.run(bwa_index_cmd, check=True)
            # Track sidecar files for cleanup in circular mode
            if circular:
                for ext in ('.amb', '.ann', '.bwt', '.pac', '.sa', '.0123', '.bwt.2bit.64'):
                    sidecar = Path(str(work_assembly) + ext)
                    if sidecar.exists():
                        bwa_index_sidecars.append(sidecar)

        # Build mapper command
        reads_args = []
        if interleaved:
            reads_args = [str(interleaved)]
        elif read1 and read2:
            reads_args = [str(read1), str(read2)]
        elif read1:
            reads_args = [str(read1)]

        extra_params = []
        if mapper.startswith('minimap2') and minimap2_params:
            extra_params = shlex.split(minimap2_params)
        elif mapper == 'bwa-mem2' and bwa_params:
            extra_params = shlex.split(bwa_params)

        if mapper == 'minimap2-sr':
            mapper_cmd = ["minimap2", "-ax", "sr", "-t", str(threads)] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'minimap2-sr-secondary':
            mapper_cmd = [
                "minimap2", "-a", "-k21", "-w11", "--sr", "--frag=yes", "-A2", "-B8",
                "-O12,32", "-E2,1", "-r100", "-p.5", "-N20", "-f1000,5000", "-n2",
                "-m25", "-s40", "-g100", "-2K50m", "--heap-sort=yes", "--secondary=yes",
                "-t", str(threads),
            ] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'minimap2-ont':
            mapper_cmd = ["minimap2", "-ax", "map-ont", "-t", str(threads)] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'minimap2-pb':
            mapper_cmd = ["minimap2", "-ax", "map-pb", "-t", str(threads)] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'minimap2-hifi':
            mapper_cmd = ["minimap2", "-ax", "map-hifi", "-t", str(threads)] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'minimap2-no-preset':
            mapper_cmd = ["minimap2", "-a", "-t", str(threads)] + extra_params + [str(work_assembly)] + reads_args
        elif mapper == 'bwa-mem2':
            mapper_cmd = ["bwa-mem2", "mem", "-t", str(threads)]
            if interleaved:
                mapper_cmd.append("-p")
            mapper_cmd += extra_params + [str(work_assembly)] + reads_args
        else:
            raise ValueError(f"Unknown mapper: {mapper}")

        sorted_bam = Path(tempfile.mkstemp(prefix=output_file.stem + "_sorted_", suffix=".bam")[1])

        # Pipe: mapper | samtools view -bS -F 4 | samtools sort -o sorted_bam
        view_cmd = ["samtools", "view", "-@", str(threads), "-bS", "-"]
        if not keep_unmapped:
            view_cmd[-1:-1] = ["-F", "4"]
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

        # Inject @PG and @CO metadata headers
        cmd_parts = [mapper, str(assembly)]
        if read1:
            cmd_parts.extend(["-r1", str(read1)])
        if read2:
            cmd_parts.extend(["-r2", str(read2)])
        if interleaved:
            cmd_parts.extend(["--interleaved", str(interleaved)])
        if circular:
            cmd_parts.append("--circular")
        command_line = "thebigbam mapping-per-sample " + " ".join(cmd_parts)
        _inject_bam_headers(output_file, circular, threads, command_line)

        # Re-index after header rewrite
        subprocess.run(["samtools", "index", "-@", str(threads), str(output_file)], check=True)

        # Convert doubled-reference circular BAM to SAM-spec circular BAM
        if circular:
            convert_circular_bam(output_file, output_file, threads=threads)

    finally:
        # cleanup temporary files
        for f in temp_files:
            try:
                Path(f).unlink()
            except OSError:
                pass
        for f in bwa_index_sidecars:
            try:
                f.unlink()
            except OSError:
                pass

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Map reads and produce MD-tagged indexed BAM")
    add_mapping_per_sample_args(parser)
    args = parser.parse_args()
    run_mapping_per_sample(args)
