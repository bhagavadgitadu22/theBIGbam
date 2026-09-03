# Developer notes

## Prerequisites

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

## Build and install (for users)

If you only want the latest unpackaged code for theBIGbam, you can install the tool from the GitHub repository:

```bash
git clone https://github.com/bhagavadgitadu22/theBIGbam
cd theBIGbam
pip install .
```

## Build and install (for users)

Thank you for considering helping us making theBIGbam a better place!

To start contributing, you first need to fork the [repository](https://github.com/bhagavadgitadu22/theBIGbam) on your GitHub account and create a feature branch (`git checkout -b my-feature`). 

To install theBIGbam with the development dependencies:

```bash
python -m pip install -e ".[dev]"
python -m playwright install chromium
```

The first command includes theBIGbam normal install plus dependencies like pytest, Ruff, maturin, and Playwright. 

The second command installs the Chromium browser used by plotting performance and browser-interaction tests. On systems where Chromium's shared libraries are missing, an administrator may need to install the corresponding OS packages.

## Coding (for users)

After making your changes to the code, you will need to:

- [Compile the tool](#compiling-with-maturin) with maturin and test everything works
- Run fast tests (the default): `pytest tests/` or `pytest -m fast tests/`
- Run mapping/database integration tests: `pytest -m integration tests/`
- Run real-browser tests: `pytest -m browser tests/`
- Run every test category: `pytest -o addopts="" tests/`
- [Monitor the time performance](#monitor-the-performance) after your changes 

Finally when happy with your code, you can submit a pull request from your fork to the main GitHub repository with a clear description of your changes.

### Compiling with maturin

After installing maturin, you need to compile the tool a first time:

```bash
maturin develop --release
```

Later on, you will only need to recompile theBIGbam when modifying the Rust files.

**Tips:**

- Building the tool can be quite long if building from scratch. To speed up repeated builds, set a persistent target directory so Cargo does not recompile everything: 
  
  ```bash
  export CARGO_TARGET_DIR=~/.cargo-target/thebigbam
  ```

- For faster compilation without optimization, use `maturin develop` without `--release`

- If working from an HPC cluster, you might need to load llvm (open-source compiler)  if available with `module load llvm`

### Monitor the performance

theBIGbam was designed to enable the exploration of large alignment files, making performance a critical consideration during development. Developer options are available to assess how code changes affect the performance of the `calculate` and `serve` operations.

The `--time` flag can be used with both `calculate` and `serve` to record the execution time of individual processing steps in a log file.

For `serve`, the `--scenario` parameter can be given the path to a new text file. During the localhost session, theBIGbam records the sequence of user actions, such as filtering contigs with a mean coverage >100, selecting a contig, or displaying its primary-read coverage. This recorded sequence of actions is called a **scenario**.

Scenarios can be inspected with the `describe-scenario` utility and replayed with `replay-scenario`. The latter reports timing statistics for the recorded actions, allowing `serve` performance to be compared before and after code changes using an identical and reproducible sequence of user interactions.
