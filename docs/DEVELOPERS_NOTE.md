Need module load llvm for HPC cluster
I work in base environment (where I also installed minimap2...)

To serve plots in cluster:
```bash
mgfeatureviewer serve \
  --db examples/outputs/AKIRA/akira.db \
  --port 5006
ssh -N -L 5006:localhost:5006 boutroux@jed.epfl.ch
```

To avoid recompiling everything everytime, I had to change the target dir for cargo builds:
```bash
export CARGO_TARGET_DIR=~/.cargo-target/mgfeatureviewer
# OR maybe this helped: export PATH='/home/boutroux/.duckdb/cli/latest':$PATH
maturin develop --release
```

To interrogate SQL databases, use commands like:
```bash
duckdb examples/outputs/phageterm_long_short_and_nextera/phageterm_long_short_and_nextera_duckdb.db "SELECT * FROM Explicit_phage_mechanisms;"
```

CI (best practice)
- run: cargo test
- run: maturin develop --release
- run: pytest tests/ -v

For detailed tests:
pytest tests/test_pipeline.py::test_calculate_linear_bams -v -s
  The -s flag will show the print output so we can see the SAM file size, unsorted BAM size, etc.
