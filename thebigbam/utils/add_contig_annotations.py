import argparse
import csv
import os
import shutil
import sys
import urllib.parse
from pathlib import Path

from Bio import SeqIO

GENBANK_EXTS = ('.gbk', '.gbff', '.gb', '.genbank')
GFF_EXTS = ('.gff', '.gff3')
ANNOTATION_EXTS = GENBANK_EXTS + GFF_EXTS
CSV_EXTS = ('.csv',)

MULTI_VALUE_SEP = "^"


def _sanitize_value(val, warn_cb):
    if "^" in val:
        warn_cb(f"value contains reserved separator '^' — replacing with '_': {val!r}")
        return val.replace("^", "_")
    return val


EPILOG = """\
CSV format (wide / spreadsheet style):
  Columns listed in --match-by are used to locate existing features.
  All other columns are new qualifiers to add. Empty cells are skipped.

Special column names for --match-by:
  "feature_type"  Matches the feature type keyword (CDS, gene, tRNA, ...).
                  In GenBank this is the keyword before the location
                  (not a /qualifier). In GFF3 this is column 3.
  "locus"         Matches the record/sequence identifier.
                  In GenBank this is the LOCUS name. In GFF3 this is
                  column 1 (seqid). Useful to restrict matching to a
                  specific contig.

  All other column names match /qualifier values in GenBank files or
  attribute keys in the GFF3 column 9.

Example GenBank feature:
  LOCUS       TTVDVOOI   10532 bp   DNA
  ...
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

  To restrict to a specific contig, add "locus" to --match-by:
    --match-by locus,feature_type,locus_tag
  and add a "locus" column to the CSV with the LOCUS/seqid name.

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
                             'Special column names: "feature_type" (CDS, gene, ...), '
                             '"locus" (LOCUS name in GenBank / seqid in GFF3)')
    parser.add_argument('--match-by', dest='match_by', default=None,
                        help='Comma-separated column names used for feature matching. '
                             'Special keywords: "feature_type" matches the feature type (CDS, gene, tRNA, ...); '
                             '"locus" matches the record/sequence name (LOCUS in GenBank, seqid in GFF3). '
                             'Other column names match /qualifier values in GenBank or attributes in GFF3. '
                             'Default: first CSV column')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file or directory for modified annotations')
    parser.add_argument('--force', action='store_true', default=False,
                        help='Replace existing qualifiers instead of keeping the first value')
    parser.add_argument('--keep-multiple', dest='keep_multiple', action='store_true', default=False,
                        help='When a qualifier already exists on a feature, append the new value '
                             'with "^" separator instead of keeping the first value. '
                             'Warning: links between qualifiers from the same CSV row are lost. '
                             'Use --force to overwrite instead.')
    parser.add_argument('--prefix', dest='prefix', default=None,
                        help='String prepended to every qualifier name written to the output '
                             '(e.g. --prefix mydb_ writes CSV column "category" as "mydb_category"). '
                             'Default: no prefix')


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


def _get_genbank_match_key(feature, match_by_cols, record_name):
    """Extract match-by values from a GenBank feature. Returns tuple or None if any key missing."""
    vals = []
    for col in match_by_cols:
        if col == "feature_type":
            vals.append(feature.type)
        elif col == "locus":
            vals.append(record_name)
        elif col in feature.qualifiers:
            vals.append(feature.qualifiers[col][0])
        else:
            return None
    return tuple(vals)


def _process_genbank_file(annot_path, csv_lookup, match_by_cols, value_cols, force, prefix=None, keep_multiple=False):
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
        "warnings": [],
        "modified": False,
    }

    key_to_features = {}
    for rec in records:
        for feat in rec.features:
            stats["total_features"] += 1
            key = _get_genbank_match_key(feat, match_by_cols, rec.name)
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
                    out_col = (prefix + col) if prefix else col
                    val = _sanitize_value(
                        val, lambda msg: stats["warnings"].append(f"{row_source}: {msg}"))
                    if out_col in feat.qualifiers:
                        key_desc = ", ".join(
                            f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
                        if force:
                            feat.qualifiers[out_col] = [val]
                        elif keep_multiple:
                            existing = feat.qualifiers[out_col][0]
                            feat.qualifiers[out_col] = [existing + MULTI_VALUE_SEP + val]
                            stats["warnings"].append(
                                f"{row_source}: qualifier '{out_col}' already exists on "
                                f"{key_desc} — appending value (use --force to overwrite)")
                        else:
                            stats["warnings"].append(
                                f"{row_source}: qualifier '{out_col}' already exists on "
                                f"{key_desc} — keeping first value "
                                f"(use --keep-multiple to append, --force to overwrite)")
                            continue
                    else:
                        feat.qualifiers[out_col] = [val]
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


def _process_gff_file(annot_path, csv_lookup, match_by_cols, value_cols, force, prefix=None, keep_multiple=False):
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
        "warnings": [],
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
        seqid = parts[0]
        feature_type = parts[2]
        attrs = _parse_gff_attributes(parts[8])
        attr_dict = {k: v for k, v in attrs}
        feature_lines.append((i, parts, seqid, feature_type, attrs, attr_dict))
        stats["total_features"] += 1

    key_to_feature_indices = {}
    for idx, (line_idx, parts, seqid, feature_type, attrs, attr_dict) in enumerate(feature_lines):
        vals = []
        valid = True
        for col in match_by_cols:
            if col == "feature_type":
                vals.append(feature_type)
            elif col == "locus":
                vals.append(seqid)
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
                line_idx, parts, seqid, feature_type, attrs, attr_dict = feature_lines[fidx]
                feat_modified = False

                for col in value_cols:
                    val = row_values.get(col, "")
                    if not val:
                        continue
                    out_col = (prefix + col) if prefix else col
                    val = _sanitize_value(
                        val, lambda msg: stats["warnings"].append(f"{row_source}: {msg}"))
                    existing_idx = next(
                        (i for i, (k, _) in enumerate(attrs) if k == out_col), None)
                    if out_col in attr_dict:
                        key_desc = ", ".join(
                            f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
                        if force:
                            encoded_val = urllib.parse.quote(val, safe=" /:.^*$@!+?")
                            if existing_idx is not None:
                                attrs[existing_idx] = (out_col, encoded_val)
                            attr_dict[out_col] = encoded_val
                        elif keep_multiple:
                            existing_decoded = urllib.parse.unquote(attr_dict[out_col])
                            combined = existing_decoded + MULTI_VALUE_SEP + val
                            encoded_combined = urllib.parse.quote(combined, safe=" /:.^*$@!+?")
                            if existing_idx is not None:
                                attrs[existing_idx] = (out_col, encoded_combined)
                            attr_dict[out_col] = encoded_combined
                            stats["warnings"].append(
                                f"{row_source}: qualifier '{out_col}' already exists on "
                                f"{key_desc} — appending value (use --force to overwrite)")
                        else:
                            stats["warnings"].append(
                                f"{row_source}: qualifier '{out_col}' already exists on "
                                f"{key_desc} — keeping first value "
                                f"(use --keep-multiple to append, --force to overwrite)")
                            continue
                    else:
                        encoded_val = urllib.parse.quote(val, safe=" /:.^*$@!+?")
                        attrs.append((out_col, encoded_val))
                        attr_dict[out_col] = encoded_val
                    stats["qualifiers_added"][col] += 1
                    feat_modified = True

                if feat_modified:
                    parts[8] = _rebuild_gff_attributes(attrs)
                    lines[line_idx] = "\t".join(parts) + "\n"
                    stats["features_modified"] += 1
                    stats["modified"] = True

    return lines, stats


def _process_one_file(annot_path, csv_lookup, match_by_cols, value_cols, force, prefix=None, keep_multiple=False):
    lower = annot_path.lower()
    if lower.endswith(GENBANK_EXTS):
        records, stats = _process_genbank_file(annot_path, csv_lookup, match_by_cols, value_cols, force, prefix, keep_multiple)
        return "genbank", records, stats
    elif lower.endswith(GFF_EXTS):
        lines, stats = _process_gff_file(annot_path, csv_lookup, match_by_cols, value_cols, force, prefix, keep_multiple)
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

    prefix = args.prefix or ""
    out_cols = [(prefix + col) for col in value_cols]

    print(f"Matching features by: {match_by_cols}", flush=True)
    print(f"Qualifiers to add: {out_cols}", flush=True)
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
    all_warnings = []
    total_features = 0
    total_features_modified = 0
    total_multi_match_rows = 0
    total_multi_match_features = 0
    total_qualifiers_added = {col: 0 for col in value_cols}
    total_files_copied = 0
    pending_copies = []  # input paths of unmodified dir-input files (strings only)

    def _handle_result(annot_path, fmt, data, stats):
        nonlocal total_features, total_features_modified, total_files_copied
        nonlocal total_multi_match_rows, total_multi_match_features

        total_features += stats["total_features"]
        total_features_modified += stats["features_modified"]
        total_multi_match_rows += stats["multi_match_rows"]
        total_multi_match_features += stats["multi_match_features"]
        all_matched_keys.update(stats["matched_keys"])
        all_errors.extend(stats["errors"])
        all_warnings.extend(stats["warnings"])
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
            if is_dir_input:
                pending_copies.append(annot_path)
                print(f"  {Path(annot_path).name}: no modifications, deferred", flush=True)
            else:
                print(f"  {Path(annot_path).name}: no modifications, not written", flush=True)

    for af in annot_files:
        fmt, data, stats = _process_one_file(af, csv_lookup, match_by_cols, value_cols, args.force, prefix or None, args.keep_multiple)
        _handle_result(af, fmt, data, stats)

    if pending_copies:
        if total_features_modified == 0:
            print("\nError: no modifications made in any file; output not written.", flush=True)
            return 1
        for src in pending_copies:
            dst = output_path / Path(src).name
            shutil.copy2(src, dst)
            total_files_copied += 1
            print(f"  {Path(src).name}: no modifications, copied", flush=True)

    unmatched_keys = set(csv_lookup.keys()) - all_matched_keys
    for key in unmatched_keys:
        for row_source, _ in csv_lookup[key]:
            key_desc = ", ".join(f"{mc}={kv}" for mc, kv in zip(match_by_cols, key))
            all_errors.append(f"{row_source}: no feature found for {key_desc}")

    print(f"\nSummary:", flush=True)
    print(f"  Features annotated: {total_features_modified} / {total_features} total features", flush=True)
    if total_files_copied > 0:
        print(f"  Files copied unchanged: {total_files_copied}", flush=True)
    if total_multi_match_rows > 0:
        print(f"  Rows with multiple matches: {total_multi_match_rows} "
              f"(affected {total_multi_match_features} features)", flush=True)
    print(f"  Qualifiers added:", flush=True)
    for col, out_col in zip(value_cols, out_cols):
        print(f"    {out_col}: {total_qualifiers_added[col]}", flush=True)

    if all_warnings:
        print(f"\nWarnings ({len(all_warnings)}):", flush=True)
        for w in all_warnings:
            print(f"  {w}", flush=True)

    if all_errors:
        print(f"\nErrors ({len(all_errors)}):", flush=True)
        for err in all_errors:
            print(f"  {err}", flush=True)
        return 1

    return 0
