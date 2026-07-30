#!/usr/bin/env python3
"""Check that the two benchmark arms produced equivalent results.

Compares, per library / library group:

  * .dedup.stats and merged group .stats -- as key -> value maps. Never byte-diffed:
    `pairtools stats --merge` has a non-reproducible line order upstream.
  * .cool                                -- total counts, nnz and per-pixel values.
  * .nodups.scaling.tsv                  -- numerically, within tolerance.
  * the deduplicated pairs themselves    -- contacts (chrom1/pos1/chrom2/pos2/
    strand1/strand2) present in one arm but not the other, plus a breakdown of
    same-contact rows that differ only in a read-dependent column (mapq1,
    mapq2, pair_type, ...).

With `pairtools.dedup_backend: duckdb` (the default) a small number of *contacts*
are *expected* to differ: where a chain of near-duplicates straddles a chunk
boundary, pairtools drops the link and reports the tail as unique while duckdb
keeps the chain. Upstream measures roughly 3 rows per million. Differences up to
--dedup-tolerance (as a fraction of the total) are therefore reported but not
treated as failures.

Separately, and much more commonly, duckdb and pairtools agree on every contact
but keep a different member of each duplicate family as the surviving row. That
row carries that read's own mapq1/mapq2 (and, in principle, pair_type) -- not a
disagreement about the contact, just about which read represents it. This is
reported but never counted toward --dedup-tolerance. It does mean that a filter
applied before binning (e.g. `mapq_30`) can let a different subset of contacts
through in each arm, which is why coolers built from a mapq-filtered stream are
not guaranteed to be pixel-identical even though the contacts underneath agree.

    python benchmarking/compare_outputs.py [bench_dir] [--dedup-tolerance 1e-4]

Exits non-zero if contacts differ beyond that.
"""

import argparse
import gzip
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

PAIRS_LIBRARY = "results/pairs_library"
STATS_GROUP = "results/stats_library_group"
COOLERS_LIBRARY = "results/coolers_library"

problems = []
notes = []


def fail(message):
    problems.append(message)
    print(f"  FAIL  {message}")


def ok(message):
    print(f"  ok    {message}")


def note(message):
    notes.append(message)
    print(f"  note  {message}")


# --------------------------------------------------------------------------- stats


def parse_stats(path):
    """A pairtools .stats file is `key<TAB>value` lines."""
    values = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        key, value = parts[0], parts[-1]
        try:
            values[key] = float(value)
        except ValueError:
            values[key] = value
    return values


def compare_stats(a_path, b_path, label):
    a, b = parse_stats(a_path), parse_stats(b_path)
    only_a = sorted(set(a) - set(b))
    only_b = sorted(set(b) - set(a))
    if only_a or only_b:
        fail(f"{label}: key sets differ (only baseline: {only_a[:5]}, only candidate: {only_b[:5]})")
        return
    differing = {k: (a[k], b[k]) for k in a if a[k] != b[k]}
    if not differing:
        ok(f"{label}: identical ({len(a)} keys)")
        return
    total = a.get("total") or 1.0
    headline = {k: v for k, v in differing.items() if k in ("total", "total_dups", "total_nodups", "cis", "trans")}
    worst = max(
        (abs(v[0] - v[1]) / total for v in differing.values() if isinstance(v[0], float) and isinstance(v[1], float)),
        default=0.0,
    )
    detail = ", ".join(f"{k}: {v[0]:g} vs {v[1]:g}" for k, v in list(headline.items() or differing.items())[:6])
    note(f"{label}: {len(differing)} of {len(a)} keys differ (max {worst:.2e} of total) -- {detail}")


# --------------------------------------------------------------------------- pairs


def read_pairs(path, converter_env=None):
    """Return (body rows, column names) of a .pairs.gz / .parquet file.

    Parquet is converted to text with `pairtools_parquet parquet-to-csv` so both arms
    are compared in the same representation.
    """
    if path.suffix == ".parquet":
        cmd = ["pairtools_parquet", "parquet-to-csv", str(path)]
        if converter_env:
            cmd = converter_env + cmd
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        text = out.stdout
    else:
        with gzip.open(path, "rt") as handle:
            text = handle.read()

    columns = []
    rows = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].split()
            continue
        rows.append(line)
    return rows, columns


# Columns that identify the read rather than the contact. Dedup is free to keep a
# different member of a duplicate family, which changes these and nothing else.
READ_ID_COLUMNS = {"readID", "parent_readID"}

# The contact itself: what dedup groups duplicates by. Everything else in a
# .pairs row (pair_type, mapq1, mapq2, ...) describes the specific read that was
# kept as a duplicate family's representative, which the two backends are free
# to disagree on -- see the module docstring.
CONTACT_COLUMNS = {"chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"}


def split_columns(columns):
    """Split column indices into (contact-identity, read-dependent-attribute)."""
    contact_idx = [i for i, name in enumerate(columns) if name in CONTACT_COLUMNS]
    attribute_idx = [
        i
        for i, name in enumerate(columns)
        if name not in CONTACT_COLUMNS and name not in READ_ID_COLUMNS
    ]
    return contact_idx, attribute_idx


def _select(fields, idx):
    return tuple(fields[i] for i in idx if i < len(fields))


def _index_by_contact(rows, contact_idx, attribute_idx):
    by_contact = defaultdict(list)
    for row in rows:
        fields = row.split("\t")
        by_contact[_select(fields, contact_idx)].append(_select(fields, attribute_idx))
    return by_contact


def compare_pairs(a_path, b_path, label, tolerance, converter_env):
    try:
        a_rows, a_cols = read_pairs(a_path, converter_env)
        b_rows, b_cols = read_pairs(b_path, converter_env)
    except (subprocess.CalledProcessError, OSError) as exc:
        fail(f"{label}: could not read pairs ({exc})")
        return

    if a_cols != b_cols:
        fail(f"{label}: columns differ ({a_cols} vs {b_cols})")
        return

    total = max(len(a_rows), 1)
    contact_idx, attribute_idx = split_columns(a_cols)
    attribute_names = [a_cols[i] for i in attribute_idx]

    a_by_contact = _index_by_contact(a_rows, contact_idx, attribute_idx)
    b_by_contact = _index_by_contact(b_rows, contact_idx, attribute_idx)

    # Contacts missing on one side entirely -- the real disagreement, and what
    # --dedup-tolerance is measured against.
    only_a = 0
    only_b = 0
    # Contacts present in both. `identical` rows agree on every column;
    # `differing` rows agree on the contact but not on some read-dependent
    # column -- expected when dedup kept a different duplicate-family member.
    identical = 0
    differing = 0
    per_column_diffs = Counter()

    for key in set(a_by_contact) | set(b_by_contact):
        a_attrs = sorted(a_by_contact.get(key, []))
        b_attrs = sorted(b_by_contact.get(key, []))
        matched = min(len(a_attrs), len(b_attrs))
        only_a += len(a_attrs) - matched
        only_b += len(b_attrs) - matched
        for a_row, b_row in zip(a_attrs[:matched], b_attrs[:matched]):
            if a_row == b_row:
                identical += 1
            else:
                differing += 1
                for name, a_val, b_val in zip(attribute_names, a_row, b_row):
                    if a_val != b_val:
                        per_column_diffs[name] += 1

    fraction = (only_a + only_b) / total

    if only_a == 0 and only_b == 0:
        if differing:
            detail = ", ".join(f"{name}: {n}" for name, n in per_column_diffs.most_common())
            ok(
                f"{label}: {len(a_rows)} rows, contacts identical -- "
                f"{identical} rows fully identical, {differing} same contact but "
                f"differ in a read-dependent column ({detail}) -- dedup kept a "
                f"different member of the duplicate family"
            )
        else:
            ok(f"{label}: {len(a_rows)} rows, identical")
        return

    message = (
        f"{label}: {len(a_rows)} vs {len(b_rows)} rows; contacts differ -- "
        f"{only_a} only in baseline, {only_b} only in candidate "
        f"({fraction:.2e} of total; of the rest, {differing} same contact but "
        f"differ in a read-dependent column)"
    )
    if fraction <= tolerance:
        note(message + " -- within the documented dedup divergence")
    else:
        fail(message)


# -------------------------------------------------------------------------- coolers


def compare_coolers(a_path, b_path, label):
    try:
        import cooler
        import numpy as np
    except ImportError:
        note(f"{label}: cooler/numpy not importable here, skipping cooler comparison")
        return

    a, b = cooler.Cooler(str(a_path)), cooler.Cooler(str(b_path))
    a_info, b_info = a.info, b.info
    if a.binsize != b.binsize:
        fail(f"{label}: binsize {a.binsize} vs {b.binsize}")
        return
    if a_info["nnz"] != b_info["nnz"] or a_info["sum"] != b_info["sum"]:
        note(
            f"{label}: nnz {a_info['nnz']} vs {b_info['nnz']}, "
            f"sum {a_info['sum']} vs {b_info['sum']}"
        )
    a_pixels = a.pixels()[:]
    b_pixels = b.pixels()[:]
    merged = a_pixels.merge(
        b_pixels, on=["bin1_id", "bin2_id"], how="outer", suffixes=("_a", "_b")
    ).fillna(0)
    diff = np.abs(merged["count_a"] - merged["count_b"])
    if diff.max() == 0:
        ok(f"{label}: {a_info['nnz']} pixels, identical")
    else:
        note(
            f"{label}: {int((diff > 0).sum())} of {len(merged)} pixels differ, "
            f"max delta {int(diff.max())}, total delta {int(diff.sum())}"
        )


# -------------------------------------------------------------------------- scaling


def compare_scaling(a_path, b_path, label):
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        note(f"{label}: pandas/numpy not importable here, skipping scaling comparison")
        return

    a = pd.read_table(a_path)
    b = pd.read_table(b_path)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        fail(f"{label}: shape/columns differ ({a.shape} vs {b.shape})")
        return
    numeric = a.select_dtypes("number").columns
    close = np.allclose(a[numeric].fillna(0), b[numeric].fillna(0), rtol=1e-9, atol=0)
    if close:
        ok(f"{label}: {len(a)} rows, identical")
    else:
        delta = (a[numeric].fillna(0) - b[numeric].fillna(0)).abs().max().max()
        note(f"{label}: numeric values differ, max delta {delta:g}")


# ------------------------------------------------------------------------------ main


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bench_dir", nargs="?", default="bench", type=Path)
    parser.add_argument("--baseline", default="pairtools")
    parser.add_argument("--candidate", default="parquet")
    parser.add_argument(
        "--dedup-tolerance",
        type=float,
        default=1e-4,
        help="fraction of rows allowed to differ before it counts as a failure",
    )
    parser.add_argument(
        "--converter-prefix",
        default="",
        help="command prefix used to run pairtools_parquet (e.g. a conda run wrapper)",
    )
    args = parser.parse_args()

    a_dir = args.bench_dir / args.baseline
    b_dir = args.bench_dir / args.candidate
    for d in (a_dir, b_dir):
        if not d.exists():
            sys.exit(f"No such arm directory: {d}")
    converter_env = args.converter_prefix.split() if args.converter_prefix else None

    print("\n== dedup stats ==")
    for a_path in sorted((a_dir / PAIRS_LIBRARY).glob("*.dedup.stats")):
        b_path = b_dir / PAIRS_LIBRARY / a_path.name
        if b_path.exists():
            compare_stats(a_path, b_path, a_path.name)
        else:
            fail(f"{a_path.name}: missing in {args.candidate}")

    print("\n== merged group stats ==")
    for a_path in sorted((a_dir / STATS_GROUP).glob("*.stats")):
        b_path = b_dir / STATS_GROUP / a_path.name
        if b_path.exists():
            compare_stats(a_path, b_path, a_path.name)
        else:
            fail(f"{a_path.name}: missing in {args.candidate}")

    print("\n== deduplicated pairs ==")
    for a_path in sorted((a_dir / PAIRS_LIBRARY).glob("*.nodups.pairs.gz")):
        stem = a_path.name[: -len(".nodups.pairs.gz")]
        b_path = b_dir / PAIRS_LIBRARY / f"{stem}.nodups.parquet"
        if not b_path.exists():
            b_path = b_dir / PAIRS_LIBRARY / a_path.name
        if b_path.exists():
            compare_pairs(a_path, b_path, f"{stem}.nodups", args.dedup_tolerance, converter_env)
        else:
            fail(f"{stem}.nodups: missing in {args.candidate}")

    print("\n== coolers ==")
    for a_path in sorted((a_dir / COOLERS_LIBRARY).glob("*.cool")):
        b_path = b_dir / COOLERS_LIBRARY / a_path.name
        if b_path.exists():
            compare_coolers(a_path, b_path, a_path.name)
        else:
            fail(f"{a_path.name}: missing in {args.candidate}")

    print("\n== scaling ==")
    for a_path in sorted((a_dir / PAIRS_LIBRARY).glob("*.scaling.tsv")):
        b_path = b_dir / PAIRS_LIBRARY / a_path.name
        if b_path.exists():
            compare_scaling(a_path, b_path, a_path.name)
        else:
            fail(f"{a_path.name}: missing in {args.candidate}")

    print()
    print("=" * 70)
    if problems:
        print(f"{len(problems)} problem(s), {len(notes)} note(s)")
        for message in problems:
            print(f"  FAIL  {message}")
        sys.exit(1)
    print(f"No problems. {len(notes)} note(s) -- differences within expected bounds.")


if __name__ == "__main__":
    main()
