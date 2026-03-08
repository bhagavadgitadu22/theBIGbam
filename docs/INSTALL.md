# Installation Guide

## Option 1: conda (includes mapping tools)

```bash
conda install -c bioconda thebigbam
```

This installs theBIGbam along with **samtools**, **minimap2**, and **bwa-mem2** — everything needed for both BAM processing and read mapping.

---

## Option 2: pip

Pre-built wheels are available for Linux and macOS — no Rust toolchain needed:

```bash
pip install thebigbam
```

---

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

This compiles the Rust code and installs all Python dependencies.

### For developers

```bash
git clone https://github.com/bhagavadgitadu22/theBIGbam
cd theBIGbam
pip install maturin
# On HPC clusters: you may need to load the LLVM module first: `module load llvm`
maturin develop --release
```

With this installations:

- Changes to Python code take effect immediately without reinstalling

- The Rust code is still compiled, but you can recompile it quickly with `maturin develop` after making changes

---

## External mapping tools

3 mapping packages are needed if you want to use `thebigbam mapping-per-sample` to map reads:

- **samtools** — BAM file manipulation
- **minimap2** — Read aligner (long reads or short reads)
- **bwa-mem2** — Alternative read aligner (short reads)

They are installed directly with option 1. Otherwise they can be installed separately via conda:

```bash
conda install -c bioconda samtools minimap2 bwa-mem2
```

Or install individually from their official sources.

---

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
