#!/usr/bin/env python3
"""Check that the two benchmark arms produced equivalent results.

Compares, per library / library group:

  * .dedup.stats and merged group .stats -- as key -> value maps. Never byte-diffed:
    `pairtools stats --merge` has a non-reproducible line order upstream.
  * .cool                                -- total counts, nnz and per-pixel values.
  * .nodups.scaling.tsv                  -- numerically, within tolerance.
  * the deduplicated pairs themselves    -- row counts, and the number of rows
    present in one arm but not the other.

With `pairtools.dedup_backend: duckdb` (the default) a small number of rows are
*expected* to differ: where a chain of near-duplicates straddles a chunk boundary,
pairtools drops the link and reports the tail as unique while duckdb keeps the chain.
Upstream measures roughly 3 rows per million. Differences up to --dedup-tolerance
(as a fraction of the total) are therefore reported but not treated as failures.

    python benchmarking/compare_outputs.py [bench_dir] [--dedup-tolerance 1e-4]

Exits non-zero if anything differs beyond that.
"""

import argparse
import gzip
import subprocess
import sys
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


def read_pairs_rows(path, converter_env=None):
    """Return the set of body rows of a .pairs.gz / .parquet file, and the row count.

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
    rows = [line for line in text.splitlines() if line and not line.startswith("#")]
    return rows


def compare_pairs(a_path, b_path, label, tolerance, converter_env):
    try:
        a_rows = read_pairs_rows(a_path, converter_env)
        b_rows = read_pairs_rows(b_path, converter_env)
    except (subprocess.CalledProcessError, OSError) as exc:
        fail(f"{label}: could not read pairs ({exc})")
        return

    a_set, b_set = set(a_rows), set(b_rows)
    only_a = len(a_set - b_set)
    only_b = len(b_set - a_set)
    total = max(len(a_rows), 1)

    if only_a == 0 and only_b == 0 and len(a_rows) == len(b_rows):
        ok(f"{label}: {len(a_rows)} rows, identical")
        return

    fraction = (only_a + only_b) / total
    message = (
        f"{label}: {len(a_rows)} vs {len(b_rows)} rows; "
        f"{only_a} only in baseline, {only_b} only in candidate "
        f"({fraction:.2e} of total)"
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
