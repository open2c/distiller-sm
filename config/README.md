# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

## Input

For each sample the input files are typically a pairs of .fastq.gz files, one with forward and one with reverse reads.
Additionally, the input can be specified as an accession in the SRA database, and reads will be downloaded automatically.

For each biological sample multiple technical replicates ("lanes") can be provided, they are then merged at the stage of pairs.

Biological samples (e.g. biological replicates) can also be grouped into "library groups", so they are merged at the level of coolers.

You need to provide the name of the genome assembly, the path to the bwa index with a wildcard, and to the chromsizes file.
The index doesn't need to already exist, as long as provided path matches exactly the fasta file with the reference genome (e.g. sequence in mm10.fa.gz, provide mm10.fa.gz*). If the index doesn't exist, it will be created.

## Mapping

Mapping can be done with bwa-mem, bwa-mem2, bwa-meme (all produce identical or near-identical results), or chromap.
Chromap outputs .pairs directly and works very fast, but you lose the flexibility of custom parising options.

## Pairs processing (`pairtools:` section)

Every pairs-processing step -- parse, sort, merge, dedup, select, scaling, stats -- runs
with one of two interchangeable implementations, chosen by `pairtools.backend`:

- `pairtools`: the classic pipeline. Intermediates are bgzipped text (`.pairs.gz`), and
  the deduplicated library pairs are indexed with `pairix`.
- `parquet`: [`pairtools_parquet`](https://github.com/Phlya/pairtools_parquet).
  Intermediates are `.parquet` and are streamed between steps as Arrow IPC, so the data
  is never serialised to text in between.

Only the intermediate format changes; the coolers, `.dedup.stats`, scaling tables and
MultiQC report are the same either way, and the `bin.filters` expressions are the same
Python syntax in both.

The remaining options only apply when `backend: parquet`:

| option | meaning |
| --- | --- |
| `dedup_backend` | `duckdb` (default, fastest), or `scipy`/`sklearn` for bit-exact parity with `pairtools dedup`. pairtools' `cython` backend does not exist here. |
| `memory` | memory budget handed to DuckDB for sort/merge/dedup. |
| `tmpdir` | where DuckDB spills, and where `dedup` spools its input when piped into. Point this at fast local scratch on a cluster. |
| `export_text_pairs` | additionally write `.nodups.pairs.gz` + pairix index. Off by default; nothing in the workflow needs it. |

`duckdb` deliberately differs from pairtools by roughly 3 rows per million: where a
chain of near-duplicates straddles a chunk boundary, pairtools drops the link and
reports the tail as unique, while duckdb keeps the chain intact. Use `scipy` if you
need byte-identical output. See `benchmarking/README.md`.

Note that with `backend: parquet` no `pairix` index is produced -- parquet cannot be
bgzf-indexed. Nothing downstream in this workflow reads one; enable
`export_text_pairs` if an external tool needs `.pairs.gz`.

### Merging chunks and runs

Each run's chunks normally go through two merge passes before dedup: chunks merge
into a per-run file (`{run}.pairs.gz`/`.parquet`, kept under
`output.dirs.pairs_runs`), then `merge_dedup` merges across a library's runs. Set
`pairtools.merge_runs_before_dedup: False` to skip the per-run merge and feed every
chunk from every run of a library into `merge_dedup` directly -- one merge pass
instead of `n_runs + 1`, and no per-run intermediate written. The final library
pairs are identical either way; this only changes how many merge passes happen and
whether a per-run file exists.

If you still want the per-run files as a separate artifact while using the direct
path, set `pairtools.keep_run_pairs: True` too -- they'll be produced as an extra
workflow target rather than as merge_dedup's input.