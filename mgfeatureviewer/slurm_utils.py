import shutil
import subprocess
from pathlib import Path
from typing import Optional

def has_slurm() -> bool:
    return shutil.which("sbatch") is not None

def submit_sbatch(script_path: Path) -> str:
    """Submit an sbatch script and return the job id string."""
    res = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True, check=True)
    out = res.stdout.strip()
    # sbatch typically prints: "Submitted batch job 123456"
    parts = out.split()
    jobid = parts[-1] if parts else out
    return jobid

def submit_array(csv_path: Path, outdir: Path, module_cli: str, array_size: int, concurrency: int, cpus_per_task: int,
                 mem: str = "8G", time: str = "02:00:00", columns: Optional[list] = None) -> str:
    """
    Create and submit a job-array script that runs a Python CLI module for each row in `csv_path`.

    - `module_cli` is the Python module/command to run for a single task (e.g. "mgfeatureviewer.cli mapping-per-sample ...").
    - `array_size` is the number of tasks (rows).
    - `concurrency` is max concurrent tasks (slurm array %limit).
    Returns the sbatch job id.
    """
    script = outdir / "sbatch_array.sh"
    out_log = outdir / "sbatch_%A_%a.out"

    # The script: it computes the row index from SLURM_ARRAY_TASK_ID, extracts the CSV row using Python,
    # then exposes columns as shell variables and executes the provided module CLI line.
    cols = columns or ['assembly', 'read1', 'read2', 'seqtype']
    read_vars = ' '.join(cols)

    content = f"""#!/bin/bash
#SBATCH --job-name=mgfv_array
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={out_log}
#SBATCH --array=1-{array_size}%{concurrency}

IDX=$SLURM_ARRAY_TASK_ID
# Extract the CSV row for this task (1-based index)
ROW=$(python - <<'PY'
import csv,sys
idx=int(sys.argv[1])
csvp=sys.argv[2]
with open(csvp,newline='') as fh:
    r=list(csv.reader(fh))
print('|||'.join(r[idx-1]))
PY
 $IDX {csv_path})

IFS='|||' read -r {read_vars} <<< "$ROW"

# Build an 'output' variable if read1 and assembly exist (common case)
if [ -n "$read1" ] && [ -n "$assembly" ]; then
  readbase=$(basename "$read1")
  readbase=${{readbase%.*}}
  asmbase=$(basename "$assembly")
  asmbase=${{asmbase%.*}}
  output={outdir}/"${{readbase}}_mapped_on_${{asmbase}}.bam"
fi

# Run the requested CLI command; module_cli should be a string that may reference the shell variables above
CMD="{module_cli}"
eval $CMD
"""
    
    script.write_text(content)
    # Submit with sbatch
    jobid = submit_sbatch(script)
    return jobid
