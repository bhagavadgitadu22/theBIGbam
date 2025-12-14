# MGFeatureViewer Installation Guide

MGFeatureViewer is a hybrid Python/Rust tool that combines fast Rust-based BAM processing with interactive Python visualization. Installation requires three components:

1. **Rust toolchain** (to compile the fast calculation engine)
2. **Python environment** (for the CLI and visualization)
3. **External tools** (optional, only needed if you want to do mapping/annotation)

## Detailed Installation

### Step 1: Install Rust 1.7+

Rust is needed to compile the fast calculation engine. This is a one-time setup.

**Linux, macOS, or WSL:**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env  # Add Rust to your PATH
```

**Windows:**
1. Download [rustup-init.exe](https://rustup.rs/)
2. Run the installer and follow prompts
3. Restart your terminal

**Verify Rust installation:**
```bash
rustc --version
cargo --version
```

### Step 2: Install Python 3.9+

Check if you have Python installed:
```bash
python --version  # or python3 --version
```

If not installed, download from [python.org](https://www.python.org/downloads/) or use your system package manager.

### Step 3: Install MGFeatureViewer

**For users (recommended):**
```bash
git clone https://github.com/bhagavadgitadu22/MGFeatureViewer
cd MGFeatureViewer
pip install .
```

This will:
- Compile the Rust code (takes 5-10 minutes first time)
- Install Python dependencies (bokeh, biopython, pysam, etc.)
- Create the `mgfeatureviewer` command

**For developers:**
```bash
git clone https://github.com/bhagavadgitadu22/MGFeatureViewer
cd MGFeatureViewer
pip install maturin
maturin develop --release
```

This installs in "editable" mode - changes to Python code take effect immediately without reinstalling.

### Step 4: (Optional) Install External Tools

**Only needed if you want to:**
- Map reads to assemblies (minimap2, samtools)
- Annotate genomes (pharokka, bakta)

**If you already have BAM files and GenBank annotations, skip this step.**

**Using conda/mamba (easiest):**
```bash
mamba env create -f mgfeatureviewer_env.yaml
conda activate mgfeatureviewer
pip install .  # Re-install in conda environment
```

This installs:
- **minimap2** - Read aligner
- **samtools** - BAM file manipulation
- **seqkit** - Sequence file manipulation
- **pharokka** or **bakta** - Genome annotation tools

Alternatively, install these tools individually via your system package manager or from their official sources.

## Verify Installation

```bash
# Check main command works
mgfeatureviewer -h

# Run quick test with example data
mgfeatureviewer calculate \
  -g examples/inputs/HK97/HK97_GCF_000848825.1_pharokka.gbk \
  -b examples/inputs/HK97/ \
  -m coverage \
  -o test.db \
  --circular \
  -t 2

# Visualize interactively the test data
mgfeatureviewer serve --db examples/outputs/HK97/test.db --port 5006
# Open browser to http://localhost:5006
```

## Common Issues

**"rust-htslib" compilation errors:**
- On Linux: Install development headers: `sudo apt-get install libbz2-dev liblzma-dev zlib1g-dev`
- On macOS: `brew install xz bzip2 zlib`

**"maturin: command not found" or build fails:**
- Make sure Rust is installed: `rustc --version`
- Install maturin: `pip install maturin`
- Try: `maturin develop --release`

**Python import errors:**
- Ensure you're in the correct environment if using conda
- Try reinstalling: `pip uninstall mgfeatureviewer && pip install .`
- Or force rebuild: `pip install --force-reinstall --no-cache-dir .`

**Slow compilation:**
- First-time compilation takes 5-10 minutes
- Subsequent installs are faster (~1 minute) due to caching
- For less optimization but faster compile use: `maturin develop  # Without --release`

## Updating

```bash
cd MGFeatureViewer
git pull
pip install --force-reinstall .
```

For developers:
```bash
git pull
maturin develop --release
```

## Uninstalling

```bash
pip uninstall mgfeatureviewer
```

The Rust toolchain and external tools (conda environment) remain installed and can be removed separately if desired.
