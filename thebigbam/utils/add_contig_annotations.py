import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

GENBANK_EXTS = ('.gbk', '.gbff', '.gb', '.genbank')
GFF_EXTS = ('.gff', '.gff3')
ANNOTATION_EXTS = GENBANK_EXTS + GFF_EXTS
CSV_EXTS = ('.csv',)


EPILOG = """\
CSV format (wide / spreadsheet style):
  Columns listed in --match-by are used to locate existing features.
  All other columns are new qualifiers to add. Empty cells are skipped.

  The special column name "feature_type" matches the feature type
  (e.g. CDS, gene, tRNA) which in GenBank files appears as the keyword
  before the location (not as a /qualifier).

Example GenBank feature:
      CDS             4028..5185
                     /locus_tag="TTVDVOOI_CDS_0005"
                     /product="major head protein"
                     /phrog="10"

  To add a "category" qualifier to this CDS, create a CSV file:

    feature_type,locus_tag,category
    CDS,TTVDVOOI_CDS_0005,structural

  Then run:
    thebigbam add-contig-annotations \\
        -g annotations.gbk \\
        --csv new_qualifiers.csv \\
        --match-by feature_type,locus_tag \\
        -o annotations_enriched.gbk

  Here --match-by feature_type,locus_tag means: find the feature whose
  type is "CDS" AND whose /locus_tag is "TTVDVOOI_CDS_0005", then add
  /category="structural" to it.

  If --match-by is omitted, the first CSV column is used for matching.
"""


def add_add_contig_annotations_args(parser):
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.epilog = EPILOG
    parser.add_argument('-g', '--genbank', required=True,
                        help='Path to annotation FILE or DIRECTORY (.gbk, .gbff, .gb, .genbank, .gff, .gff3)')
    parser.add_argument('--csv', dest='csv_file', required=True,
                        help='CSV file (or directory of CSV files) with new annotations in wide format. '
                             'Columns in --match-by identify features; all other columns are qualifiers to add. '
                             'Use "feature_type" as a column name to match by feature type (CDS, gene, etc.)')
    parser.add_argument('--match-by', dest='match_by', default=None,
                        help='Comma-separated column names used for feature matching. '
                             'Use "feature_type" to match the feature type (CDS, gene, tRNA, ...). '
                             'Other column names match /qualifier values in GenBank or attributes in GFF3. '
                             'Default: first CSV column')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file or directory for modified annotations')
    parser.add_argument('--force', action='store_true', default=False,
                        help='Replace existing qualifiers instead of erroring')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of annotation files to process in parallel (default: 1)')


def _list_dir(directory, exts):
    return sorted(
        os.path.join(directory, f) for f in os.listdir(directory)
        if f.lower().endswith(exts)
    )


def _resolve_input(path, exts, label):
    p = Path(path)
    if not p.exists():
        sys.exit(f"ERROR: {label} path not found: {path}")
    if p.is_file():
        if not p.name.lower().endswith(exts):
            sys.exit(f"ERROR: Unsupported {label} file format: {p.name}. Supported: {', '.join(exts)}")
        return [str(p)]
    files = _list_dir(str(p), exts)
    if not files:
        sys.exit(f"ERROR: No {label} files found in directory '{path}'. Supported: {', '.join(exts)}")
    return files


def _load_csv_files(csv_files, match_by_cols):
    """Parse CSV files into a list of annotation rows.

    Returns (rows, match_by_cols, value_cols) where each row is a dict with
    all CSV columns.  match_by_cols/value_cols are determined from the first
    file's header.
    """
    rows = []
    value_cols = None

    for csv_path in csv_files:
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames
            if header is None:
                sys.exit(f"ERROR: CSV file '{csv_path}' is empty or has no header")

            if match_by_cols is None:
                match_by_cols = [header[0]]

            missing = [c for c in match_by_cols if c not in header]
            if missing:
                sys.exit(f"ERROR: Match-by columns not found in '{csv_path}' header: {', '.join(missing)}")

            file_value_cols = [c for c in header if c not in match_by_cols]
            if value_cols is None:
                value_cols = file_value_cols
            elif file_value_cols != value_cols:
                print(f"WARNING: CSV file '{csv_path}' has different value columns than first file; "
                      f"using union of all columns", flush=True)
                for c in file_value_cols:
                    if c not in value_cols:
                        value_cols.append(c)

            for i, row in enumerate(reader, start=2):
                row['_source'] = f"{csv_path}:{i}"
                rows.append(row)

    return rows, match_by_cols, value_cols


def _get_feature_match_key(feature, match_by_cols, fmt):
    """Extract match-by values from a feature. Returns tuple or None if any key missing."""
    vals = []
    for col in match_by_cols:
        if fmt == "genbank":
            if col == "feature_type":
                vals.append(feature.type)
            elif col in feature.qualifiers:
                vals.append(feature.qualifiers[col][0])
            else:
                return None
        else:
            if col == "feature_type":
                vals.append(feature.get("feature_type", ""))
            elif col in feature.get("attributes", {}):
                vals.append(feature["attributes"][col])
            else:
                return None
    return tuple(vals)


def _process_genbank_file(annot_path, csv_lookup, match_by_cols, value_cols, force):
    """Process a single GenBank file. Returns (output_records, stats)."""
    records = list(SeqIO.parse(annot_path, "genbank"))
    stats = {
        "file": annot_path,
        "total_features": 0,
        "features_modified": 0,
        "qualifiers_added": {col: 0 for col in value_cols},
        "matched_keys": set(),
        "multi_match_rows": 0,
        "multi_match_features": 0,
        "errors": [],
        "modified": False,
    }

    key_to_features = {}
    for rec in records:
        for feat in rec.features:
            stats["total_features"] += 1
            key = _get_feature_match_key(feat, match_by_cols, "genbank")
            if key is not None:
                key_to_features.setdefault(key, []).append(feat)

    for key, row_list in csv_lookup.items():
        if key not in key_to_features:
            continue

        features = key_to_features[key]
        is_multi = len(features) > 1

        for row_source, row_values in row_list:
            stats["matched_keys"].add(key)
            if is_multi:
                stats["multi_match_rows"] += 1
                stats["multi_match_features"] += len(features)

            for feat in features:
                feat_modified = False
                for col in value_cols:
                    val = row_values.get(col, "")
                    if not val:
                        continue
                    if col in feat.qualifiers and not force:
                        key_desc = ", ".join(f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
                        stats["errors"].append(
                            f"{row_source}: qualifier '{col}' already exists on {key_desc} "
                            f"(use --force to overwrite)")
                        continue
                    feat.qualifiers[col] = [val]
                    stats["qualifiers_added"][col] += 1
                    feat_modified = True
                if feat_modified:
                    stats["features_modified"] += 1
                    stats["modified"] = True

    return records, stats


def _parse_gff_attributes(attr_str):
    """Parse GFF3 column 9 into an ordered list of (key, value) tuples."""
    if not attr_str or attr_str == ".":
        return []
    pairs = []
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            pairs.append((k, v))
        else:
            pairs.append((item, ""))
    return pairs


def _rebuild_gff_attributes(pairs):
    return ";".join(f"{k}={v}" if v else k for k, v in pairs)


def _process_gff_file(annot_path, csv_lookup, match_by_cols, value_cols, force):
    """Process a single GFF3 file. Returns (output_lines, stats)."""
    stats = {
        "file": annot_path,
        "total_features": 0,
        "features_modified": 0,
        "qualifiers_added": {col: 0 for col in value_cols},
        "matched_keys": set(),
        "multi_match_rows": 0,
        "multi_match_features": 0,
        "errors": [],
        "modified": False,
    }

    with open(annot_path) as fh:
        lines = fh.readlines()

    feature_lines = []
    for i, line in enumerate(lines):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feature_type = parts[2]
        attrs = _parse_gff_attributes(parts[8])
        attr_dict = {k: v for k, v in attrs}
        feature_lines.append((i, parts, feature_type, attrs, attr_dict))
        stats["total_features"] += 1

    key_to_feature_indices = {}
    for idx, (line_idx, parts, feature_type, attrs, attr_dict) in enumerate(feature_lines):
        vals = []
        valid = True
        for col in match_by_cols:
            if col == "feature_type":
                vals.append(feature_type)
            elif col in attr_dict:
                vals.append(attr_dict[col])
            else:
                valid = False
                break
        if valid:
            key = tuple(vals)
            key_to_feature_indices.setdefault(key, []).append(idx)

    for key, row_list in csv_lookup.items():
        if key not in key_to_feature_indices:
            continue

        feat_indices = key_to_feature_indices[key]
        is_multi = len(feat_indices) > 1

        for row_source, row_values in row_list:
            stats["matched_keys"].add(key)
            if is_multi:
                stats["multi_match_rows"] += 1
                stats["multi_match_features"] += len(feat_indices)

            for fidx in feat_indices:
                line_idx, parts, feature_type, attrs, attr_dict = feature_lines[fidx]
                feat_modified = False

                for col in value_cols:
                    val = row_values.get(col, "")
                    if not val:
                        continue
                    if col in attr_dict and not force:
                        key_desc = ", ".join(f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
                        stats["errors"].append(
                            f"{row_source}: qualifier '{col}' already exists on {key_desc} "
                            f"(use --force to overwrite)")
                        continue

                    existing_idx = next((i for i, (k, _) in enumerate(attrs) if k == col), None)
                    if existing_idx is not None:
                        attrs[existing_idx] = (col, val)
                    else:
                        attrs.append((col, val))
                    attr_dict[col] = val
                    stats["qualifiers_added"][col] += 1
                    feat_modified = True

                if feat_modified:
                    parts[8] = _rebuild_gff_attributes(attrs)
                    lines[line_idx] = "\t".join(parts) + "\n"
                    stats["features_modified"] += 1
                    stats["modified"] = True

    return lines, stats


def _process_one_file(annot_path, csv_lookup, match_by_cols, value_cols, force):
    lower = annot_path.lower()
    if lower.endswith(GENBANK_EXTS):
        records, stats = _process_genbank_file(annot_path, csv_lookup, match_by_cols, value_cols, force)
        return "genbank", records, stats
    elif lower.endswith(GFF_EXTS):
        lines, stats = _process_gff_file(annot_path, csv_lookup, match_by_cols, value_cols, force)
        return "gff", lines, stats
    else:
        sys.exit(f"ERROR: Unsupported annotation format: {annot_path}")


def _write_output(fmt, data, output_path):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    if fmt == "genbank":
        with open(output_path, "w") as fh:
            SeqIO.write(data, fh, "genbank")
    else:
        with open(output_path, "w") as fh:
            fh.writelines(data)


def run_add_contig_annotations(args):
    annot_files = _resolve_input(args.genbank, ANNOTATION_EXTS, "annotation")
    csv_files = _resolve_input(args.csv_file, CSV_EXTS, "CSV")

    match_by_cols = [c.strip() for c in args.match_by.split(",")] if args.match_by else None

    rows, match_by_cols, value_cols = _load_csv_files(csv_files, match_by_cols)

    if not value_cols:
        sys.exit("ERROR: No qualifier columns to add (all columns are match-by columns)")
    if not rows:
        sys.exit("ERROR: No data rows found in CSV file(s)")

    print(f"Matching features by: {match_by_cols}", flush=True)
    print(f"Qualifiers to add: {value_cols}", flush=True)
    print(f"CSV rows loaded: {len(rows)}", flush=True)
    print(f"Annotation files to process: {len(annot_files)}", flush=True)

    csv_lookup = {}
    for row in rows:
        key = tuple(row.get(c, "") for c in match_by_cols)
        row_values = {c: row.get(c, "") for c in value_cols}
        csv_lookup.setdefault(key, []).append((row['_source'], row_values))

    is_dir_input = Path(args.genbank).is_dir()
    output_path = Path(args.output)

    all_matched_keys = set()
    all_errors = []
    total_features = 0
    total_features_modified = 0
    total_multi_match_rows = 0
    total_multi_match_features = 0
    total_qualifiers_added = {col: 0 for col in value_cols}

    def _handle_result(annot_path, fmt, data, stats):
        nonlocal total_features, total_features_modified
        nonlocal total_multi_match_rows, total_multi_match_features

        total_features += stats["total_features"]
        total_features_modified += stats["features_modified"]
        total_multi_match_rows += stats["multi_match_rows"]
        total_multi_match_features += stats["multi_match_features"]
        all_matched_keys.update(stats["matched_keys"])
        all_errors.extend(stats["errors"])
        for col in value_cols:
            total_qualifiers_added[col] += stats["qualifiers_added"][col]

        if stats["modified"]:
            if is_dir_input:
                out = output_path / Path(annot_path).name
            else:
                out = output_path
            _write_output(fmt, data, str(out))
            print(f"  {Path(annot_path).name}: modified {stats['features_modified']} features, "
                  f"added {sum(stats['qualifiers_added'].values())} qualifiers", flush=True)
        else:
            print(f"  {Path(annot_path).name}: no modifications, skipped", flush=True)

    if args.threads > 1 and len(annot_files) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(_process_one_file, af, csv_lookup, match_by_cols, value_cols, args.force): af
                for af in annot_files
            }
            for future in as_completed(futures):
                af = futures[future]
                fmt, data, stats = future.result()
                _handle_result(af, fmt, data, stats)
    else:
        for af in annot_files:
            fmt, data, stats = _process_one_file(af, csv_lookup, match_by_cols, value_cols, args.force)
            _handle_result(af, fmt, data, stats)

    unmatched_keys = set(csv_lookup.keys()) - all_matched_keys
    for key in unmatched_keys:
        for row_source, _ in csv_lookup[key]:
            key_desc = ", ".join(f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
            all_errors.append(f"{row_source}: no feature found for {key_desc}")

    print(f"\nSummary:", flush=True)
    print(f"  Features annotated: {total_features_modified} / {total_features} total features", flush=True)
    if total_multi_match_rows > 0:
        print(f"  Rows with multiple matches: {total_multi_match_rows} "
              f"(affected {total_multi_match_features} features)", flush=True)
    print(f"  Qualifiers added:", flush=True)
    for col in value_cols:
        print(f"    {col}: {total_qualifiers_added[col]}", flush=True)

    if all_errors:
        print(f"\nErrors ({len(all_errors)}):", flush=True)
        for err in all_errors:
            print(f"  {err}", flush=True)
        return 1

    return 0
