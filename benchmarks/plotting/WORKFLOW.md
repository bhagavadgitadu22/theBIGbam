# Plotting performance benchmark

This benchmark is the repeatable end-to-end plotting scenario used when evaluating
server, filtering, rendering, and browser-performance changes. It is an internal
developer benchmark, not user documentation.

## Reference dataset and selections

- Database: `examples/thebigbam_nomis_mags_filtered_at_90_75_24_06_2026.db`
- MAG: `337R_metabat_172`
- Restricted sample: `GL34_UP`
- Filters: `Annotations.activity = Defense` AND `Coverage.Coverage_mean > 100`
- MAG ordering: `Coverage.Coverage_mean`, descending
- All-samples variable: `Primary alignments`
- Sample ordering: `Coverage.Coverage_mean`, descending

## Procedure

Record this reference workflow, then replay it once with a newly started server
(cold caches) and optionally again against the same server (warm caches):

1. Open the Defense-value lookup and select `Defense`.
2. Add the coverage row, open its histogram lookup, and set coverage above 100.
3. Apply filters and record the resulting MAG/contig/sample availability.
4. Enter MAG VIEW, select `337R_metabat_172`, and select `GL34_UP`.
5. Enable Coverage and Misalignment, then render the one-sample MAG.
6. Order contigs by Coverage mean and render again.
7. Add `activity = Defense` coloring rules to the gene map and MAG overview, then render.
8. Enter ALL SAMPLES and measure the transition.
9. APPLY with no variable selected and measure the validation response.
10. Select Primary alignments and render the restricted sample set.
11. Relax the coverage row, apply filters, and render the larger sample set.
12. Order samples by Coverage mean, enable a common y scale, and render again.
13. Select a Defense-associated contig, then manually verify double-click copy,
    zoom, and hover behavior.

Recorded settings changes, filter lookups, **APPLY FILTERS**, and **APPLY** each
produce elapsed browser time, matching server timing lines, a screenshot, and any
failure. Clipboard and hover gestures are not currently scenario actions and remain
manual checks. Results are written as JSON and TSV under the selected output
directory.

## Running

```bash
python benchmarks/plotting/run_workflow.py \
  --scenario benchmarks/plotting/scenarios/my_workflow.json \
  --db examples/thebigbam_nomis_mags_filtered_at_90_75_24_06_2026.db \
  --cold-and-warm \
  --output benchmarks/plotting/results
```

Use `--url http://localhost:5006` to benchmark an already-running server. In this
mode the runner does not start or stop the server and cannot capture its stdout.
The scenario is required. `--db` is required only when the runner starts the
server. Before opening a browser, the runner validates the scenario format,
actions, step numbers, and (for a locally started server) the recorded database
filename. Unsupported semantic actions fail clearly instead of being skipped.
Settings changes are replayed through the application's normal settings
restoration boundary, so dynamic AND/OR sections, variables, feature controls,
annotation/MAG coloring rules, and every plotting parameter use the same
ontology and compatibility behavior as **SAVE SETTINGS**. When the runner owns
the server, skipped database-specific options are promoted from restoration
warnings to failed benchmark steps.

Each run writes a scenario-shaped result (`results.cold.json` and, when
requested, `results.warm.json`). Its steps retain their original `sequence`,
`action`, `changes`, and `details`, augmented with `status`,
`duration_seconds`, `error`, `artifacts`, and matching server diagnostic lines.
`results.tsv` provides the same measurements for tabular comparison.

## Recording a reusable scenario

The hidden `serve --scenario` option records a manual plotting session as an
ordered, settings-shaped JSON document:

```bash
thebigbam serve \
  --db examples/thebigbam_nomis_mags_filtered_at_90_75_24_06_2026.db \
  --time \
  --scenario benchmarks/plotting/scenarios/nomis_mags.json
```

The document uses the same ontology as **SAVE SETTINGS** for view mode,
selection, filtering, variables, annotation configuration, and plotting
parameters. It additionally records semantic operations that settings alone do
not represent, including filter lookups, **APPLY FILTERS**, and **APPLY**.
Scenario files describe intent only: they contain no execution status or timing.
Benchmark result JSON copies this structure and augments its numbered steps with
status, duration, errors, and artifacts from that particular run.

The complete JSON document is replaced atomically after every recorded revision.
It is therefore safe to inspect after a Slurm cancellation or process kill: a
forced termination may omit the latest interaction, but it does not leave
truncated JSON. Correctness does not depend on final process cleanup.

If multiple browser sessions connect to the same recording server, the first
uses the requested filename and later sessions use `.session-2`, `.session-3`,
and so on, so their actions are never interleaved.

To identify a slow or problematic step without reading the JSON, run the hidden
scenario-description utility:

```bash
thebigbam describe-scenario benchmarks/plotting/scenarios/nomis_mags.json
```

It prints exactly one deterministic line per step. The number at the start of
each line is the stable `sequence` stored in the scenario, so it can be quoted
in benchmark reports and discussions. When the input is an augmented benchmark
result rather than a plain scenario, the same lines also show recorded status
and duration.
