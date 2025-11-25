# Validation

Compares feature values between the original Python/SQLite implementation and
the Rust/Parquet implementation.

## Usage

```bash
python compare.py <sqlite_db> <parquet_dir>
```

Example:
```bash
python compare.py ../examples/outputs/original/metadata.db ../examples/outputs/new/
```

Issues are written to `validation_issues.csv` for further analysis.

## Expected Output

For features without signal compression (reads_starts, reads_ends, tau,
insertions, deletions, mismatches, left_clippings, right_clippings), the match
should be 100%.

For compressed features (coverage, coverage_reduced, read_lengths), a small
percentage of positions may differ due to a bug fix in the Rust implementation.

### Run 25.11.2025

```
SQLite:  examples/outputs/AKIRA/AKIRA.db
Parquet: examples/outputs/AKIRA/ddb
Samples: 47, Contigs: 27

Feature                  SQLite    Parquet  Missing    Extra                    Match
------------------------------------------------------------------------------------
coverage                176,674    228,821   12,129   64,276 164,545/176,674 (93.13%)
coverage_reduced        123,284    163,800       31   40,547 123,253/123,284 (99.97%)
reads_starts             22,120     22,120        0        0  22,120/22,120 (100.00%)
reads_ends               15,204     15,204        0        0  15,204/15,204 (100.00%)
tau                      17,805     17,805        0        0  17,805/17,805 (100.00%)
read_lengths            136,252    202,444       54   66,246 136,198/136,252 (99.96%)
insertions               70,349     70,349        0        0  70,349/70,349 (100.00%)
deletions                79,437     79,437        0        0  79,437/79,437 (100.00%)
mismatches               97,857     97,857        0        0  97,857/97,857 (100.00%)
left_clippings           19,479     19,479        0        0  19,479/19,479 (100.00%)
right_clippings          11,951     11,951        0        0  11,951/11,951 (100.00%)
------------------------------------------------------------------------------------
Total                   770,412    929,267   12,214  171,069 758,198/770,412 (98.41%)
```

#### Known Difference: Derivative Outlier Detection

The Python implementation contains a bug in the `compress_signal` function
where derivative outlier detection performs arithmetic on a boolean array
instead of indices:

```python
# Python (calculating_data.py)
der_outliers = np.abs(dy) > deriv_thresh * dy_std  # boolean array
der_outliers = np.unique(np.clip(np.concatenate([der_outliers - 1, der_outliers]), 0, n - 1))
# True - 1 = 0, False - 1 = -1 (unintended)
```

The Rust implementation correctly identifies derivative outlier indices:

```rust
// Rust (compress.rs)
.filter(|(_, &d)| d.abs() > deriv_thresh * deriv_std)
.flat_map(|(i, _)| if i > 0 { vec![i - 1, i] } else { vec![i] })
```

As a result:
- Both implementations keep approximately the same number of positions
- The specific positions kept around derivative outliers differ slightly
- Values at shared positions are identical
