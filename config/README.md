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