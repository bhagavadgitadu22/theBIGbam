# Contributing to theBIGbam

Thank you for your interest in contributing to theBIGbam!

## Reporting Bugs

Please open an issue on [GitHub](https://github.com/bhagavadgitadu22/theBIGbam/issues) with:

- A clear description of the problem
- Steps to reproduce
- Your OS, Python version, and theBIGbam version (`thebigbam --version`)
- Any relevant error messages or log output

## Suggesting Features

Open an issue with the "feature request" label describing:

- What you want to achieve
- Why existing functionality doesn't cover your use case

## Development Setup

1. Clone the repository:
   
   ```bash
   git clone https://github.com/bhagavadgitadu22/theBIGbam.git
   cd theBIGbam
   ```

2. Install Rust: https://rust-lang.org/tools/install

3. Create a conda environment for external tools:
   
   ```bash
   conda env create -f thebigbam_env.yaml
   conda activate thebigbam
   ```

4. Install in development mode:
   
   ```bash
   pip install -e ".[dev]"
   ```

5. Build the Rust extension:
   
   ```bash
   maturin develop --features python
   ```

## Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b my-feature`)
3. Make your changes
4. Run tests: `pytest tests/`
5. Submit a pull request with a clear description of your changes

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
