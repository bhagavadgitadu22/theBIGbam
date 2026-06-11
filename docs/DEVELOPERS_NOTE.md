# Developer Notes

## Installing from GitHub

### Prerequisites

| Requirement       | Why                                     | Install                                                                                          |
| ----------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------ |
| **Rust** (stable) | Compiles the Rust extension via maturin | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` then `source $HOME/.cargo/env` |
| **Python ≥ 3.10** | Minimum supported version               | Check with `python --version`                                                                    |

A **virtual environment** is recommended to avoid conflicts with system packages:

```bash
python -m venv thebigbam-env
source thebigbam-env/bin/activate
```

Or with conda:

```bash
conda create -n thebigbam python=3.10
conda activate thebigbam
```

### Build and install

```bash
git clone https://github.com/bhagavadgitadu22/theBIGbam
cd theBIGbam
pip install .
```

`pip install .` automatically invokes [maturin](https://www.maturin.rs/) (the build backend declared in `pyproject.toml`), compiles the Rust extension with full optimizations, and installs all Python dependencies.

### Verify

Test your install as explained in [the main page](../README.md#check-installation-succeeded)

---

## Modifying theBIGbam

Thank you for your interest in contributing to theBIGbam! To start coding, you will need to:

- Fork the repository
- Create a feature branch (`git checkout -b my-feature`)
- Make your changes
- Compile the tool with maturin and test everything works
- Run tests: `pytest tests/`
- Submit a pull request from your fork to the main GitHub repository with a clear description of your changes

### Compiling with maturin

The standard `pip install .` uses maturin as a build backend behind the scenes, but does not install the `maturin` CLI itself. For development, you need it explicitly so you can rebuild the tool with `maturin develop`:

```bash
pip install "maturin>=1.4,<2.0"
```

After running a first compilation with maturin, you will not need to recompile when modifying Python files. You will need to do so after modifying Rust files though with:

```bash
maturin develop --release
```

Building the tool can be quite long if building from scratch. To speed up repeated builds, set a persistent target directory so Cargo does not recompile everything:

```bash
export CARGO_TARGET_DIR=~/.cargo-target/thebigbam
```

If working from an HPC cluster, you might need to load llvm (open-source compiler)  if available with `module load llvm`.

For faster compilation without optimization, use `maturin develop` without `--release`.
