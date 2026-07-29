# Benchmarking the pairs backends

`distiller-sm` can run every pairs-processing step (parse, sort, merge, dedup, select,
scaling, stats) either with [`pairtools`](https://github.com/open2c/pairtools) over
gzipped text `.pairs`, or with
[`pairtools_parquet`](https://github.com/Phlya/pairtools_parquet) over `.parquet`.
Pick one with `pairtools.backend` in the config:

```yaml
pairtools:
    backend: 'parquet'   # or 'pairtools'
```

The tooling here runs the workflow both ways on the same input and compares the
results — both how fast they were and whether they agree.

## Running a benchmark

```sh
benchmarking/run_benchmark.sh -c 16 -f config/benchmark_large.yml
```

This runs the workflow twice, in `bench/pairtools/` and `bench/parquet/`, using
`snakemake --directory` so the two arms cannot overwrite each other's `results/`,
`logs/` or `benchmarks/`. Inputs are symlinked in rather than copied, and the bwa
index is built once and shared through snakemake's between-workflow cache (rule
`bwaindex` is marked `cache: True`).

Useful flags:

| flag | meaning |
| --- | --- |
| `-c, --cores N` | cores per arm (default 4) |
| `-f, --configfile FILE` | which config to run (default `config/config.yml`) |
| `-o, --outdir DIR` | where the arms go (default `bench/`) |
| `-a, --arm NAME` | run only one arm |
| `-- ...` | everything after `--` is passed to snakemake |

Then:

```sh
python benchmarking/compare.py bench           # timings and on-disk sizes
python benchmarking/compare_outputs.py bench   # do the two arms agree?
```

## Which config to benchmark with

**`config/config.yml`** (the default sacCer3 test, 4 × 100k reads) is a correctness
check, not a benchmark. Every rule finishes in seconds and the numbers are dominated
by process startup. Use it to confirm the two backends agree.

**`config/benchmark_large.yml`** pulls much deeper slices of the same SRA accessions
through the workflow's existing `sra:...?start=&end=` download path, so the pairs
steps get enough work to be measurable. sacCer3 is a 12 Mb genome, so deep
sequencing produces a high duplicate rate — which is what stresses `dedup`, the step
with the largest expected win. It also sets `map.chunksize`, which multiplies the
number of `parse | sort` pipes and gives `merge_runs` real work.

**`config/benchmark_synthetic.yml`** needs no network at all: generate the fixture
first with

```sh
python benchmarking/make_synthetic_test_data.py --outdir test
```

which writes a small random genome and paired-end reads sampled from it, with a
controlled PCR-duplicate rate. Useful for checking the wiring in a sandbox with no
access to UCSC/SRA. The reads carry no real Hi-C structure, so do not read anything
into its timings.

**Your own data is the only representative benchmark.** Point
`input.raw_reads_paths` at a real library and run the harness against that.

## Reading the output

`compare.py` sums the per-job `benchmark:` TSVs snakemake already writes for every
heavy rule, and prints wall time per rule for both arms with the speedup. It also
prints the total from `/usr/bin/time -v` around each `snakemake` invocation, which
includes scheduling overhead the per-rule TSVs miss, and the on-disk size of the
pairs intermediates — parquet's storage footprint is a real secondary benefit.

Expect the win to be concentrated in `merge_dedup` and `parse_sort_chunks`.
`parse` itself is not faster (it is pysam-bound); what the parquet backend removes
is the text serialisation between stages.

## Expect a small difference in dedup output

`compare_outputs.py` compares dedup stats, merged group stats, coolers, scaling
tables and the deduplicated pairs themselves. It reports differences rather than
asserting equality, and only exits non-zero when they exceed `--dedup-tolerance`
(default 1e-4 of rows).

That is because `pairtools_parquet dedup --backend duckdb` — the default, and the
fast one — deliberately differs from `pairtools dedup`. Where a chain of
near-duplicates straddles a chunk boundary, pairtools drops the link and reports the
tail as unique; duckdb keeps every row in its lookback window and resolves the chain.
Upstream measures roughly 3 rows per million.

To remove that difference, set:

```yaml
pairtools:
    dedup_backend: 'scipy'   # pairtools' own KD-tree, bit-exact
```

at the cost of most of dedup's speedup. `--max-mismatch 0` (this repo's default
`dedup.max_mismatch_bp`) takes a separate, order-independent code path that upstream
measures at 8.8x.

Note that merged group stats are compared as key → value maps, never byte-diffed:
`pairtools stats --merge` has a non-reproducible line order upstream, so the two arms
would not compare equal even against itself.
