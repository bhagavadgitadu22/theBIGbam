from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from thebigbam.utils.read_mapping import _resolve_assembly, _inject_bam_headers


def _has_md_tags(input_path: str, threads: int) -> bool:
    """Check whether the first aligned reads in the file carry MD tags."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", "260", "-@", str(threads), str(input_path)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    try:
        for i, line in enumerate(proc.stdout):
            if i >= 10:
                break
            fields = line.split("\t")
            if len(fields) > 11:
                for field in fields[11:]:
                    if field.startswith("MD:Z:"):
                        return True
        return False
    finally:
        proc.kill()
        proc.wait()


def add_format_mapping_args(parser):
    parser.add_argument('-i', '--input', required=True, dest='input_file',
        help='Input mapping file (SAM, BAM, or CRAM)')
    parser.add_argument('-o', '--output', required=True,
        help='Output BAM path (will be written with accompanying .bai index)')
    parser.add_argument('-t', '--threads', type=int, default=4,
        help='Threads for samtools (default: 4)')
    parser.add_argument('-a', '--assembly', default=None,
        help='Reference FASTA file or directory for samtools calmd (MD tag generation). '
             'Only needed if the input lacks MD tags.')
    parser.add_argument('--keep-unmapped', action='store_true',
        help='Keep unmapped reads in the output BAM (default: discard them)')
    parser.add_argument('--min-read-percent-identity', type=float, default=0.0,
        dest='min_read_percent_identity',
        help='Exclude reads by overall percent identity e.g. 95 for 95%%. [default: 0 (disabled)]')
    parser.add_argument('--min-read-aligned-percent', type=float, default=0.0,
        dest='min_read_aligned_percent',
        help='Exclude reads by percent aligned bases e.g. 95 means 95%% of the read\'s bases must be aligned. [default: 0 (disabled)]')


def run_format_mapping(args) -> int:
    input_path = Path(args.input_file)
    output_path = Path(args.output)
    threads = args.threads

    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", flush=True)
        return 2

    if input_path.resolve() == output_path.resolve():
        print("ERROR: Input and output paths must be different.", flush=True)
        return 2

    if shutil.which("samtools") is None:
        print("ERROR: samtools not found on PATH.", flush=True)
        return 2

    assembly_path = None
    assembly_tmp = None
    if args.assembly:
        assembly_path, assembly_tmp = _resolve_assembly(args.assembly)

    sorted_bam = Path(tempfile.mkstemp(prefix=output_path.stem + "_sorted_", suffix=".bam")[1])
    temp_files = [sorted_bam]

    try:
        # Step 1: convert to sorted BAM
        view_cmd = ["samtools", "view", "-@", str(threads), "-bS"]
        if not args.keep_unmapped:
            view_cmd.extend(["-F", "4"])
        view_cmd.append(str(input_path))
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), "-"]

        print("COMMAND_VIEW:", " ".join(view_cmd), flush=True)
        print("COMMAND_SORT:", " ".join(sort_cmd), flush=True)

        p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        subprocess.run(sort_cmd, stdin=p1.stdout, check=True)
        p1.stdout.close()
        p1.wait()

        # Step 2: apply read quality filters
        min_id = args.min_read_percent_identity
        min_al = args.min_read_aligned_percent
        identity_expr = f"(qlen-[NM])*100/qlen>={min_id}" if min_id > 0.0 else ""
        aligned_expr = f"(qlen-sclen)*100/qlen>={min_al}" if min_al > 0.0 else ""
        combined_expr = " && ".join(e for e in [identity_expr, aligned_expr] if e)

        if combined_expr:
            total_before = int(subprocess.run(
                ["samtools", "view", "-@", str(threads), "-c", str(sorted_bam)],
                capture_output=True, text=True, check=True,
            ).stdout.strip())

            filtered_bam = Path(tempfile.mkstemp(prefix=output_path.stem + "_filtered_", suffix=".bam")[1])
            temp_files.append(filtered_bam)
            filter_cmd = [
                "samtools", "view", "-@", str(threads), "-b",
                "-e", combined_expr, "-o", str(filtered_bam), str(sorted_bam),
            ]
            print("COMMAND_FILTER:", " ".join(filter_cmd), flush=True)
            subprocess.run(filter_cmd, check=True)

            total_after = int(subprocess.run(
                ["samtools", "view", "-@", str(threads), "-c", str(filtered_bam)],
                capture_output=True, text=True, check=True,
            ).stdout.strip())

            excluded = total_before - total_after
            filters_used = []
            if identity_expr:
                filters_used.append(f"--min-read-percent-identity {min_id}%")
            if aligned_expr:
                filters_used.append(f"--min-read-aligned-percent {min_al}%")
            print(f"Read filtering ({' + '.join(filters_used)}): "
                  f"{excluded}/{total_before} reads excluded ({total_after} remaining)", flush=True)

            shutil.move(str(filtered_bam), str(sorted_bam))

        # Step 3: MD tag handling
        has_md = _has_md_tags(str(sorted_bam), threads)

        if has_md:
            print("MD tags detected in input — skipping calmd.", flush=True)
            shutil.move(str(sorted_bam), str(output_path))
        elif assembly_path is not None:
            print("No MD tags detected — running samtools calmd.", flush=True)
            calmd_cmd = ["samtools", "calmd", "-@", str(threads), "-b", str(sorted_bam), str(assembly_path)]
            print("COMMAND_CALMD:", " ".join(calmd_cmd), flush=True)
            with output_path.open("wb") as outfh:
                subprocess.run(calmd_cmd, stdout=outfh, stderr=subprocess.DEVNULL, check=True)
        else:
            print("WARNING: No MD tags detected and no --assembly provided. "
                  "The output BAM may not work with all thebigbam calculate modules. "
                  "Provide -a/--assembly to generate MD tags with samtools calmd.", flush=True)
            shutil.move(str(sorted_bam), str(output_path))

        # Step 4: index
        subprocess.run(["samtools", "index", "-@", str(threads), str(output_path)], check=True)

        # Step 5: inject @PG and @CO headers
        cmd_parts = ["-i", str(input_path), "-o", str(output_path)]
        if args.assembly:
            cmd_parts.extend(["-a", str(args.assembly)])
        if args.keep_unmapped:
            cmd_parts.append("--keep-unmapped")
        if min_id > 0.0:
            cmd_parts.extend(["--min-read-percent-identity", str(min_id)])
        if min_al > 0.0:
            cmd_parts.extend(["--min-read-aligned-percent", str(min_al)])
        command_line = "thebigbam format-mapping " + " ".join(cmd_parts)
        _inject_bam_headers(output_path, circular=False, threads=threads, command_line=command_line)

        # Step 6: re-index after header rewrite
        subprocess.run(["samtools", "index", "-@", str(threads), str(output_path)], check=True)

        print(f"Done. Output: {output_path} (+{output_path}.bai)", flush=True)
        return 0

    finally:
        for f in temp_files:
            try:
                Path(f).unlink(missing_ok=True)
            except OSError:
                pass
        if assembly_tmp:
            try:
                assembly_tmp.unlink()
            except OSError:
                pass
