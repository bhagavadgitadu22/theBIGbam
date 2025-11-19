from pathlib import Path
import subprocess
import shutil
from typing import List

from Bio import SeqIO
from typing import Optional

def add_annotation_args(parser):
    parser.add_argument('--csv', required=True, help='CSV file with columns: assembly_file,read1,read2,sequencing_type')
    parser.add_argument('--assembly', help='Assembly file to use for all rows (overrides CSV field)')
    parser.add_argument('--annotation_tool', required=True, choices=['pharokka', 'bakta'], help='Annotation tool to run')
    parser.add_argument('--annotation_db', required=True, help='Annotation tool database path/identifier (passed to pharokka/bakta via --db)')
    parser.add_argument('--meta', action='store_true', help='Pass --meta flag to the annotator')
    parser.add_argument('--threads', type=int, default=4, help='Threads to pass to the annotation tool')
    parser.add_argument('--genbank', required=True, help='Final GenBank file path to write (gbk or gbff extension)')

def _collect_assembly_paths_from_csv(csv_path: Path) -> List[Path]:
    import csv

    with csv_path.open(newline='') as fh:
        reader = csv.reader(fh)
        rows = list(reader)

    # Skip header if present
    if rows:
        first = [c.lower() for c in rows[0]]
        if any('assembly' in c or 'read' in c for c in first):
            rows = rows[1:]

    assemblies = []
    for r in rows:
        if not r:
            continue
        cols = [c.strip() for c in r]
        if not cols:
            continue
        # assembly is first column
        a = cols[0]
        if a:
            assemblies.append(Path(a))
    return assemblies

def _dereplicate_by_basename(paths):
    seen = set()
    out = []
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Assembly not found: {p}")
        key = p.stem
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
    return out

def _combine_fastas(assemblies, out_path: Path):
    # Defensive coercion: accept strings or Path-like objects
    out_path = Path(out_path)
    assemblies = [Path(a) for a in assemblies]

    # Write combined fasta to out_path
    records_written = 0
    with out_path.open('w') as outfh:
        for asm in assemblies:
            for rec in SeqIO.parse(str(asm), 'fasta'):
                rec.name = rec.id
                rec.description = rec.description
                SeqIO.write(rec, outfh, 'fasta')
                records_written += 1
    if records_written == 0:
        raise RuntimeError('No sequences written to combined FASTA')

def _run_pharokka(combined_fasta: Path, out_dir: Path, threads: int, db: Optional[str] = None, meta: bool = False) -> Path:
    # pharokka invocation per requested format:
    # pharokka.py --threads <threads> --db <db> --infile <input_assembly> --outdir <output-dir> [--meta]
    out_dir = Path(out_dir)
    cmd = [
        'pharokka.py',
        '--threads', str(threads),
    ]
    if db:
        cmd.extend(['-d', str(db)])
    cmd.extend(['--infile', str(combined_fasta), '--outdir', str(out_dir)])
    if meta:
        cmd.append('--meta')
    subprocess.run(cmd, check=True)
    # search for .gbk files in out_dir
    gbks = list(out_dir.rglob('*.gbk'))
    if not gbks:
        raise RuntimeError('pharokka did not produce a .gbk file in output directory')
    # If multiple, choose the first
    return gbks[0]

def _run_bakta(combined_fasta: Path, out_dir: Path, threads: int, db: Optional[str] = None, meta: bool = False) -> Path:
    # bakta invocation per requested format:
    # bakta --threads <threads> --db <db> --output <output-dir> [--meta <input_assembly>]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / 'bakta_out'
    cmd = [
        'bakta',
        '--threads', str(threads),
    ]
    if db:
        cmd.extend(['--db', str(db)])
    cmd.extend(['--output', str(prefix)])
    if meta:
        # user requested pattern: --meta <input_assembly>
        cmd.extend(['--meta', str(combined_fasta)])
    cmd.append(str(combined_fasta))
    subprocess.run(cmd, check=True)
    # bakta may write <prefix>.gbk or outputs in out_dir; search for .gbk
    gbks = list(out_dir.rglob('*.gbff'))
    if not gbks:
        raise RuntimeError('bakta did not produce a .gbff file in output directory')
    return gbks[0]

def run_annotation(args):
    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f'CSV not found: {csv_path}')

    assemblies = _collect_assembly_paths_from_csv(csv_path)
    if getattr(args, 'assembly', None):
        assemblies = [Path(args.assembly)]

    assemblies = _dereplicate_by_basename(assemblies)

    out_gbk = Path(args.genbank)
    if out_gbk.exists():
        raise FileExistsError(f'Output GenBank file already exists: {out_gbk}')
    
    out_dir = out_gbk.parent
    combined = out_dir / 'combined_assemblies.fasta'
    _combine_fastas(assemblies, combined)

    # Deduplicate contigs by name using seqkit (rmdup -n)
    deduped = out_dir / 'combined_assemblies_deduped.fasta'
    if shutil.which('seqkit') is None:
        raise FileNotFoundError('Required executable not found on PATH: seqkit')
    subprocess.run([
        'seqkit', 'rmdup', '-n', '--threads', str(args.threads), '-o', str(deduped), str(combined)
    ], check=True)

    annotation_dir = out_dir / 'annotation_out'
    if args.annotation_tool == 'pharokka':
        gbk_path = _run_pharokka(deduped, annotation_dir, args.threads, db=getattr(args, 'annotation_db', None), meta=getattr(args, 'meta', False))
    else:
        gbk_path = _run_bakta(deduped, annotation_dir, args.threads, db=getattr(args, 'annotation_db', None), meta=getattr(args, 'meta', False))

    # copy result to desired location
    shutil.copy(gbk_path, out_gbk)

    return 0

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Combine assemblies and run annotator (pharokka or bakta)')
    add_annotation_args(parser)
    args = parser.parse_args()
    raise SystemExit(run_annotation(args))