#!/usr/bin/env python3
"""Compare SQLite (original) vs Parquet (new) feature outputs."""

import argparse
import csv
import sqlite3
from pathlib import Path

import duckdb

FEATURES = [
    "coverage", "coverage_reduced", "reads_starts", "reads_ends", "tau",
    "read_lengths", "insertions", "deletions", "mismatches",
    "left_clippings", "right_clippings"
]


def compare(sqlite_path: Path, parquet_dir: Path, verbose: bool = False) -> bool:
    sqlite_conn = sqlite3.connect(sqlite_path)
    duck = duckdb.connect()

    samples = dict(sqlite_conn.execute("SELECT Sample_id, Sample_name FROM Sample").fetchall())
    contigs = dict(sqlite_conn.execute("SELECT Contig_id, Contig_name FROM Contig").fetchall())

    print(f"SQLite:  {sqlite_path}")
    print(f"Parquet: {parquet_dir}")
    print(f"Samples: {len(samples)}, Contigs: {len(contigs)}")
    print()

    feature_stats = {f: {"sqlite": 0, "parquet": 0, "common": 0, "exact": 0} for f in FEATURES}
    value_mismatches = []
    missing_positions = []

    for sample_id, sample_name in samples.items():
        parquet_file = parquet_dir / "features" / f"{sample_name}.parquet"
        if not parquet_file.exists():
            parquet_file = parquet_dir / "features" / f"{sample_name.replace('_with_MD', '')}.parquet"
        if not parquet_file.exists():
            continue

        for feature in FEATURES:
            try:
                sqlite_data = sqlite_conn.execute(
                    f"SELECT Contig_id, Position, Value FROM Feature_{feature} WHERE Sample_id = ?",
                    (sample_id,)
                ).fetchall()
                parquet_data = duck.execute(
                    f"SELECT contig_name, position, value FROM '{parquet_file}' WHERE feature = '{feature}'"
                ).fetchall()
            except Exception:
                continue

            if not sqlite_data or not parquet_data:
                continue

            sqlite_vals = {(contigs.get(cid, ""), pos): val for cid, pos, val in sqlite_data}
            parquet_vals = {(cname, pos): val for cname, pos, val in parquet_data}

            stats = feature_stats[feature]
            stats["sqlite"] += len(sqlite_vals)
            stats["parquet"] += len(parquet_vals)

            common_keys = sqlite_vals.keys() & parquet_vals.keys()
            only_sqlite = sqlite_vals.keys() - parquet_vals.keys()

            stats["common"] += len(common_keys)

            for key in only_sqlite:
                missing_positions.append({
                    "sample": sample_name,
                    "feature": feature,
                    "contig": key[0],
                    "position": key[1],
                    "sqlite_value": sqlite_vals[key],
                })

            for key in common_keys:
                sqlite_val = sqlite_vals[key]
                parquet_val = parquet_vals[key]
                if sqlite_val == parquet_val:
                    stats["exact"] += 1
                else:
                    value_mismatches.append({
                        "sample": sample_name,
                        "feature": feature,
                        "contig": key[0],
                        "position": key[1],
                        "sqlite_value": sqlite_val,
                        "parquet_value": parquet_val,
                    })

    sqlite_conn.close()

    print(f"{'Feature':<20} {'SQLite':>10} {'Parquet':>10} {'Missing':>8} {'Extra':>8} {'Match':>24}")
    print("-" * 84)

    total_sqlite, total_parquet, total_common, total_exact = 0, 0, 0, 0
    for feature in FEATURES:
        s = feature_stats[feature]
        if s["sqlite"] == 0 and s["parquet"] == 0:
            continue
        total_sqlite += s["sqlite"]
        total_parquet += s["parquet"]
        total_common += s["common"]
        total_exact += s["exact"]

        missing = s["sqlite"] - s["common"]
        extra = s["parquet"] - s["common"]
        pct = s["exact"] / s["sqlite"] * 100 if s["sqlite"] > 0 else 0
        match_str = f"{s['exact']:,}/{s['sqlite']:,} ({pct:.2f}%)"
        print(f"{feature:<20} {s['sqlite']:>10,} {s['parquet']:>10,} {missing:>8,} {extra:>8,} {match_str:>24}")

    print("-" * 84)

    if total_common == 0:
        print("No common positions found")
        return False

    total_missing = total_sqlite - total_common
    total_extra = total_parquet - total_common
    total_pct = total_exact / total_sqlite * 100
    match_str = f"{total_exact:,}/{total_sqlite:,} ({total_pct:.2f}%)"
    print(f"{'Total':<20} {total_sqlite:>10,} {total_parquet:>10,} {total_missing:>8,} {total_extra:>8,} {match_str:>24}")

    # Write CSV if there are issues
    if missing_positions or value_mismatches:
        csv_path = Path("validation_issues.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["type", "sample", "feature", "contig", "position", "sqlite_value", "parquet_value"])
            for m in missing_positions:
                writer.writerow(["missing", m["sample"], m["feature"], m["contig"], m["position"], m["sqlite_value"], ""])
            for m in value_mismatches:
                writer.writerow(["mismatch", m["sample"], m["feature"], m["contig"], m["position"], m["sqlite_value"], m["parquet_value"]])
        print(f"\nWrote {len(missing_positions) + len(value_mismatches)} issues to {csv_path}")

    # Print summary of issues
    if missing_positions:
        print(f"\nMissing from Parquet ({len(missing_positions)}):")
        print(f"{'Sample':<35} {'Feature':<18} {'Contig':<30} {'Pos':>8} {'Value':>12}")
        print("-" * 108)
        for m in missing_positions[:10]:
            print(f"{m['sample']:<35} {m['feature']:<18} {m['contig']:<30} {m['position']:>8} {m['sqlite_value']:>12.2f}")
        if len(missing_positions) > 10:
            print(f"... and {len(missing_positions) - 10} more (see CSV)")

    if value_mismatches:
        print(f"\nValue mismatches ({len(value_mismatches)}):")
        print(f"{'Sample':<35} {'Feature':<18} {'Contig':<25} {'Pos':>8} {'SQLite':>12} {'Parquet':>12}")
        print("-" * 115)
        for m in value_mismatches[:10]:
            print(f"{m['sample']:<35} {m['feature']:<18} {m['contig']:<25} {m['position']:>8} {m['sqlite_value']:>12.4f} {m['parquet_value']:>12.4f}")
        if len(value_mismatches) > 10:
            print(f"... and {len(value_mismatches) - 10} more (see CSV)")

    return total_exact == total_sqlite


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare SQLite vs Parquet outputs")
    parser.add_argument("sqlite_db", type=Path)
    parser.add_argument("parquet_dir", type=Path)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if not args.sqlite_db.exists():
        raise SystemExit(f"Not found: {args.sqlite_db}")
    if not args.parquet_dir.exists():
        raise SystemExit(f"Not found: {args.parquet_dir}")

    compare(args.sqlite_db, args.parquet_dir, args.verbose)
