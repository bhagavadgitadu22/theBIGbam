# Developer Notes

## Building for development

Install in editable mode with maturin (Rust changes take effect after rebuild, Python changes take effect immediately):

```bash
pip install maturin
maturin develop --release
```

To speed up repeated builds, set a persistent target directory so Cargo does not recompile from scratch:

```bash
export CARGO_TARGET_DIR=~/.cargo-target/thebigbam
maturin develop --release
```

For faster compilation without optimization (useful during development):

```bash
maturin develop  # Without --release
```

## HPC cluster setup

On HPC systems, you may need to load the LLVM module before compiling:

```bash
module load llvm
```

If compilation fails with clang/libclang errors:

```bash
sudo apt install -y clang libclang-dev  # Debian/Ubuntu
```

## Serving plots remotely

To access the Bokeh server running on a remote cluster, use SSH port forwarding:

```bash
# On the cluster
thebigbam serve --db my_database.db --port 5006

# On your local machine
ssh -N -L 5006:localhost:5006 user@cluster.example.com
# Then open http://localhost:5006 in your browser
```

## Querying the database directly

You can inspect the DuckDB database directly using the DuckDB CLI:

```bash
duckdb my_database.db "SELECT * FROM Explicit_phage_mechanisms;"
duckdb my_database.db "SELECT * FROM Explicit_coverage LIMIT 10;"
duckdb my_database.db "SHOW TABLES;"
```

## Running tests

```bash
# Rust tests
cargo test

# Build Python bindings
maturin develop --release

# Python tests
pytest tests/ -v

# Single test with output
pytest tests/test_pipeline.py::test_calculate_linear_bams -v -s
```
