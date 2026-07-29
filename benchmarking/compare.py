#!/usr/bin/env python3
"""Aggregate and compare the two benchmark arms produced by run_benchmark.sh.

Snakemake already writes a per-job TSV for every rule carrying a `benchmark:`
directive (wall seconds, max_rss, cpu_time, ...). This script sums those per rule
for each arm and prints the speedup, plus the on-disk size of the pairs
intermediates -- storage is a real secondary benefit of the parquet backend.

    python benchmarking/compare.py [bench_dir] [--json out.json]
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

PAIRS_DIRS = (
    "results/mapped_parsed_sorted_chunks",
    "results/pairs_runs",
    "results/pairs_library",
)


def read_benchmark_tsv(path):
    """Return the first data row of a snakemake benchmark TSV as a dict."""
    lines = path.read_text().strip().splitlines()
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    row = {}
    for key, value in zip(header, values):
        try:
            row[key] = float(value)
        except ValueError:
            row[key] = value
    return row


def collect_arm(workdir):
    """Aggregate every benchmark TSV under `workdir` by rule name."""
    rules = {}
    bench_root = workdir / "benchmarks"
    for tsv in sorted(bench_root.rglob("*.tsv")):
        row = read_benchmark_tsv(tsv)
        if row is None:
            continue
        rule = tsv.parent.relative_to(bench_root).as_posix()
        agg = rules.setdefault(rule, {"jobs": 0, "s": 0.0, "cpu_time": 0.0, "max_rss": 0.0})
        agg["jobs"] += 1
        agg["s"] += row.get("s", 0.0) or 0.0
        cpu = row.get("cpu_time")
        agg["cpu_time"] += cpu if isinstance(cpu, float) else 0.0
        rss = row.get("max_rss")
        if isinstance(rss, float):
            agg["max_rss"] = max(agg["max_rss"], rss)
    return rules


def total_wall_clock(workdir):
    """Parse the elapsed wall time recorded by `/usr/bin/time -v`, in seconds."""
    path = workdir / "total_time.txt"
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        if "Elapsed (wall clock) time" in line:
            stamp = line.split(": ", 1)[1].strip()
            parts = [float(p) for p in stamp.split(":")]
            seconds = 0.0
            for part in parts:
                seconds = seconds * 60 + part
            return seconds
    return None


def directory_bytes(workdir):
    sizes = {}
    for rel in PAIRS_DIRS:
        target = workdir / rel
        if not target.exists():
            continue
        out = subprocess.run(
            ["du", "-sb", str(target)], capture_output=True, text=True, check=False
        )
        if out.returncode == 0:
            sizes[rel] = int(out.stdout.split()[0])
    return sizes


def human_bytes(n):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.1f}{unit}" if unit != "B" else f"{int(n)}B"
        n /= 1024
    return f"{n:.1f}TiB"


def speedup(baseline, other):
    if not other:
        return "n/a"
    return f"{baseline / other:.2f}x"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bench_dir", nargs="?", default="bench", type=Path)
    parser.add_argument("--baseline", default="pairtools")
    parser.add_argument("--candidate", default="parquet")
    parser.add_argument("--json", type=Path, help="also write raw numbers here")
    args = parser.parse_args()

    arms = {}
    for name in (args.baseline, args.candidate):
        workdir = args.bench_dir / name
        if not workdir.exists():
            sys.exit(f"No such arm directory: {workdir}")
        arms[name] = {
            "rules": collect_arm(workdir),
            "total_wall_clock_s": total_wall_clock(workdir),
            "sizes": directory_bytes(workdir),
        }

    base, cand = arms[args.baseline]["rules"], arms[args.candidate]["rules"]
    all_rules = sorted(set(base) | set(cand))

    print()
    print(f"{'rule':<38} {'jobs':>5} {args.baseline+' s':>13} {args.candidate+' s':>13} {'speedup':>9}")
    print("-" * 82)
    base_total = cand_total = 0.0
    for rule in all_rules:
        b = base.get(rule, {})
        c = cand.get(rule, {})
        b_s, c_s = b.get("s", 0.0), c.get("s", 0.0)
        base_total += b_s
        cand_total += c_s
        jobs = max(b.get("jobs", 0), c.get("jobs", 0))
        print(f"{rule:<38} {jobs:>5} {b_s:>13.2f} {c_s:>13.2f} {speedup(b_s, c_s):>9}")
    print("-" * 82)
    print(f"{'TOTAL (sum of per-rule wall time)':<38} {'':>5} {base_total:>13.2f} {cand_total:>13.2f} {speedup(base_total, cand_total):>9}")

    b_wall = arms[args.baseline]["total_wall_clock_s"]
    c_wall = arms[args.candidate]["total_wall_clock_s"]
    if b_wall and c_wall:
        print(f"{'TOTAL (snakemake wall clock)':<38} {'':>5} {b_wall:>13.2f} {c_wall:>13.2f} {speedup(b_wall, c_wall):>9}")

    print()
    print(f"{'pairs on disk':<38} {'':>5} {args.baseline:>13} {args.candidate:>13} {'ratio':>9}")
    print("-" * 82)
    for rel in PAIRS_DIRS:
        b_sz = arms[args.baseline]["sizes"].get(rel)
        c_sz = arms[args.candidate]["sizes"].get(rel)
        if b_sz is None and c_sz is None:
            continue
        ratio = f"{b_sz / c_sz:.2f}x" if b_sz and c_sz else "n/a"
        print(f"{rel:<38} {'':>5} {human_bytes(b_sz or 0):>13} {human_bytes(c_sz or 0):>13} {ratio:>9}")
    print()

    if args.json:
        args.json.write_text(json.dumps(arms, indent=2))
        print(f"Raw numbers written to {args.json}")


if __name__ == "__main__":
    main()
