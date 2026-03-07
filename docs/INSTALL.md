# Installation Guide

theBIGbam is a hybrid Python/Rust tool that combines fast Rust-based BAM processing with interactive Python visualization.

## Option 2: conda (includes mapping tools)

```bash
conda install -c bioconda thebigbam
```

This installs theBIGbam along with **samtools**, **minimap2**, and **bwa-mem2** — everything needed for both BAM processing and read mapping.

Option 1: pip (recommended)

Pre-built wheels are available for Linux and macOS — no Rust toolchain needed:

```bash
pip install thebigbam
```

**Note:** If you need read mapping capabilities (`thebigbam mapping-per-sample`), you must install **samtools**, **minimap2**, and/or **bwa-mem2** separately. See [External Tools](#external-tools) below.

## Option 3: from source

Requires the Rust toolchain and Python 3.9+.

### Step 1: Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### Step 2: Install theBIGbam

```bash
git clone https://github.com/bhagavadgitadu22/theBIGbam
cd theBIGbam
pip install .
```

This compiles the Rust code (5-10 minutes first time) and installs all Python dependencies.

### For developers

```bash
git clone https://github.com/bhagavadgitadu22/theBIGbam
cd theBIGbam
pip install maturin
maturin develop --release
```

Changes to Python code take effect immediately without reinstalling.

## External Tools

**Only needed if you want to use `thebigbam mapping-per-sample` to map reads.** If you already have BAM files, skip this.

Required tools for mapping:

- **samtools** — BAM file manipulation
- **minimap2** — Read aligner (long reads or short reads)
- **bwa-mem2** — Alternative read aligner (short reads)

Install via conda (easiest):

```bash
conda install -c bioconda samtools minimap2 bwa-mem2
```

Or install individually from their official sources.

## Verify Installation

```bash
# Check main command works
thebigbam -h

# Run quick test with example data
thebigbam calculate \
  -g tests/HK97/HK97_GCF_000848825.1_pharokka.gbk \
  -b tests/HK97/ \
  -o tests/HK97/test.db \
  -t 2

# Visualize interactively the test data
thebigbam serve --db tests/HK97/test.db --port 5006
# Open browser to http://localhost:5006
```

## Common Issues

**"rust-htslib" compilation errors:**

- On Linux: Install development headers: `sudo apt-get install libbz2-dev liblzma-dev zlib1g-dev clang libclang-dev`
- On macOS: `brew install xz bzip2 zlib`
- On HPC clusters: you may need to load the LLVM module first: `module load llvm`

**"maturin: command not found" or build fails:**

- Make sure Rust is installed: `rustc --version`
- Install maturin: `pip install maturin`
- Try: `maturin develop --release`

**Python import errors:**

- Ensure you're in the correct environment if using conda
- Try reinstalling: `pip uninstall thebigbam && pip install .`
- Or force rebuild: `pip install --force-reinstall --no-cache-dir .`

**Slow compilation:**

- First-time compilation takes 5-10 minutes
- Subsequent installs are faster (~1 minute) due to caching
- For less optimization but faster compile use: `maturin develop  # Without --release`

## Updating

**pip or conda:**

```bash
pip install --upgrade thebigbam
# or
conda update -c bioconda thebigbam
```

**From source:**

```bash
cd theBIGbam
git pull
pip install --force-reinstall .
```

## Uninstalling

```bash
pip uninstall thebigbam
```
