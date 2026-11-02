# Changelog

## 0.8.0

- Default calculations to coverage and misalignment. Other mapping modules are opt-in.
- Add BLAST scopes: `none` (default), `termini`, `contig`, `mag`, and `all`.
- Add candidate terminal-repeat searches and reuse compatible stored evidence.
- Reduce memory use during annotation loading, BAM processing, and database writes.
- Add `calculate --diagnose` timing logs and size reports. Rename `serve --time` to `--diagnose`.
- Allow `--extend` to retry incomplete samples while preserving completed samples.
- Add the Empathi annotation color template.
- Fix repeat intervals, hit identities, partner tooltips, and selection handling.
- Upgrade note: rebuild legacy undirected BLAST hit tracks before using them.

## 0.7.0

- Add session history and navigation between previous application states.
- Save and restore complete sessions, including history.

## 0.6.0

- Add APPLY FILTERS to control when filter changes take effect.
- Add scenario recording and replay benchmarks.
- Improve filtering responsiveness and interface styling.

## 0.5.0

- Save plotting settings and restore them with `serve --json`.
- Reorganize plotting code and fix server startup issues.

## 0.4.0

- Add interactive filter histograms and treemaps.
- Add MAG track markers and contig ordering within MAGs.
- Add `format-mapping` and improve mapping performance and memory use.
- Improve calculation memory handling and database writes.
- Support MAGs in CDS mapping analysis.
- Improve annotation matching and add `--keep-multiple`.
- Preserve empty samples and correct counts after data removal.

## 0.3.0

- Add MAG calculation and visualization.
- Add the RNA module and splicing features.
- Store compressed blobs with zoom levels and base-resolution chunks.
- Add annotation coloring and CDS analysis commands.
- Add annotation-file extension and mapping-data purge utilities.
- Speed up plotting, exports, and database writes.

## 0.1.5

- Update dependency version guidance and release packaging.

## 0.1.4

- Simplify installation instructions and fix the logo.

## 0.1.3

- Fix Conda build dependencies and Python compatibility constraints.

## 0.1.2

- Require Python 3.10 or newer and update Rust bindings accordingly.

## 0.1.1

- Update installation documentation and release configuration.

## 0.1.0

- Initial release of BAM processing, genome visualization, and circular mapping tools.
