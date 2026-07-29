#!/usr/bin/env python3
"""Generate a small synthetic Hi-C dataset for exercising the workflow offline.

`setup_test.sh` downloads sacCer3 and real SRA reads, which needs network access to
UCSC and NCBI. This script builds an equivalent-shaped fixture locally instead: a
small random genome plus paired-end reads sampled from it, including a controlled
fraction of PCR duplicates so that dedup has something to do.

It is a substitute for `setup_test.sh` when the real fixtures are unavailable, not a
replacement for benchmarking on real data -- the reads carry no real Hi-C structure.

    python benchmarking/make_synthetic_test_data.py --outdir test

Produces the layout config/config.yml already expects:

    test/genome/sacCer3.fa.gz
    test/genome/sacCer3.chrom.sizes
    test/fastq/<library>/<run>/<accession>_{1,2}.fastq.gz
"""

import argparse
import gzip
import random
from pathlib import Path

# Mirrors the library/run layout in config/config.yml.
LIBRARIES = {
    "MATalpha_R1": {"lane1": "SRR2601842"},
    "MATalpha_R2": {"lane1": "SRR2601845"},
    "MATa_R1": {"lane1": "SRR2601848"},
    "MATa_R2": {"lane1": "SRR2601851"},
}

BASES = "ACGT"


def make_genome(outdir, chrom_sizes, seed):
    rng = random.Random(seed)
    genome_dir = outdir / "genome"
    genome_dir.mkdir(parents=True, exist_ok=True)

    sequences = {}
    with gzip.open(genome_dir / "sacCer3.fa.gz", "wt") as fasta:
        for name, size in chrom_sizes.items():
            seq = "".join(rng.choice(BASES) for _ in range(size))
            sequences[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fasta.write(seq[i : i + 60] + "\n")

    with open(genome_dir / "sacCer3.chrom.sizes", "w") as sizes:
        for name, size in chrom_sizes.items():
            sizes.write(f"{name}\t{size}\n")

    return sequences


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sample_pair(rng, sequences, read_length, cis_fraction):
    """One Hi-C-like contact: two read ends, mostly on the same chromosome."""
    names = list(sequences)
    chrom1 = rng.choice(names)
    if rng.random() < cis_fraction:
        chrom2 = chrom1
    else:
        chrom2 = rng.choice(names)

    seq1, seq2 = sequences[chrom1], sequences[chrom2]
    pos1 = rng.randrange(0, len(seq1) - read_length)
    if chrom1 == chrom2:
        # Log-uniform separation, so contacts concentrate at short distances.
        span = len(seq2) - read_length
        offset = int(10 ** rng.uniform(2, 6))
        pos2 = min(max(pos1 + rng.choice((-1, 1)) * offset, 0), span)
    else:
        pos2 = rng.randrange(0, len(seq2) - read_length)

    read1 = seq1[pos1 : pos1 + read_length]
    read2 = revcomp(seq2[pos2 : pos2 + read_length])
    return read1, read2


def mutate(rng, read, error_rate):
    if error_rate <= 0:
        return read
    out = []
    for base in read:
        if rng.random() < error_rate:
            out.append(rng.choice(BASES.replace(base, "")))
        else:
            out.append(base)
    return "".join(out)


def write_fastqs(outdir, sequences, n_reads, read_length, dup_rate, cis_fraction,
                 error_rate, seed):
    for lib_index, (library, runs) in enumerate(LIBRARIES.items()):
        for run_index, (run, accession) in enumerate(runs.items()):
            rng = random.Random(seed + 1000 * lib_index + run_index)
            run_dir = outdir / "fastq" / library / run
            run_dir.mkdir(parents=True, exist_ok=True)

            path1 = run_dir / f"{accession}_1.fastq.gz"
            path2 = run_dir / f"{accession}_2.fastq.gz"
            quality = "I" * read_length

            with gzip.open(path1, "wt") as f1, gzip.open(path2, "wt") as f2:
                emitted = 0
                previous = None
                while emitted < n_reads:
                    if previous is not None and rng.random() < dup_rate:
                        # A PCR duplicate: same fragment, sequenced again, so it gets
                        # a different read ID but the same mapped position.
                        read1, read2 = previous
                    else:
                        read1, read2 = sample_pair(
                            rng, sequences, read_length, cis_fraction
                        )
                        previous = (read1, read2)

                    # Illumina-style read ID, so dedup's by-tile statistics have
                    # coordinates to work with.
                    tile = rng.randrange(1101, 1105)
                    x = rng.randrange(1000, 20000)
                    y = rng.randrange(1000, 20000)
                    name = f"SIM:1:FLOWCELL:1:{tile}:{x}:{y}"

                    f1.write(f"@{name} 1:N:0:1\n{mutate(rng, read1, error_rate)}\n+\n{quality}\n")
                    f2.write(f"@{name} 2:N:0:1\n{mutate(rng, read2, error_rate)}\n+\n{quality}\n")
                    emitted += 1

            print(f"  {path1}  ({n_reads} read pairs)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--outdir", type=Path, default=Path("test"))
    parser.add_argument("--n-reads", type=int, default=20000,
                        help="read pairs per run (default: 20000)")
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--chrom-size", type=int, default=400000,
                        help="length of each synthetic chromosome")
    parser.add_argument("--n-chroms", type=int, default=4)
    parser.add_argument("--dup-rate", type=float, default=0.15,
                        help="fraction of read pairs that are PCR duplicates")
    parser.add_argument("--cis-fraction", type=float, default=0.7)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
    chrom_sizes = {
        f"chr{roman[i]}": args.chrom_size for i in range(min(args.n_chroms, len(roman)))
    }

    print(f"Writing synthetic genome to {args.outdir}/genome ...")
    sequences = make_genome(args.outdir, chrom_sizes, args.seed)
    print(f"Writing synthetic fastqs to {args.outdir}/fastq ...")
    write_fastqs(
        args.outdir,
        sequences,
        args.n_reads,
        args.read_length,
        args.dup_rate,
        args.cis_fraction,
        args.error_rate,
        args.seed,
    )
    print("\nDone. Note config/config.yml lists MATalpha_R1/lane2 as an SRA accession;")
    print("comment that run out (or use config/benchmark_synthetic.yml) to stay offline.")


if __name__ == "__main__":
    main()
