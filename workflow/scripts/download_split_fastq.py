#!/usr/bin/env python

from urllib.parse import urlparse

# from common import organize_fastqs
from snakemake.shell import shell


def get_start_end(query):
    start, end = 0, None
    if query:
        for kv_pair in query.split("&"):
            k, v = kv_pair.split("=")
            if k == "start":
                start = v
            if k == "end":
                end = v
    return start, end


fastq_files = snakemake.params.fastq_files
if (len(fastq_files) == 1) and (fastq_files[0].startswith("sra:")):
    parsed = urlparse(fastq_files[0])
    srr, query = parsed.path, parsed.query
    start, end = get_start_end(query)
    if snakemake.config["map"]["use_custom_split"]:
        shell(
            f"""
            fastq-dump {srr} -Z --split-spot \
            {f'--minSpotId {start}' if start else ''} \
            {f'--maxSpotId {end}' if end else ''} \
                | python workflow/scripts/pyfilesplit --lines 4 \
                    >(bgzip -c -@{snakemake.params.bgzip_threads} > {snakemake.output.fastq1}) \
                    >(bgzip -c -@{snakemake.params.bgzip_threads} > {snakemake.output.fastq2})
            """
        )
    else:
        shell(
            f"""
            fastq-dump --origfmt --split-files --gzip \
                -O {snakemake.params.downloaded_fastqs_folder} {srr} \
                {f'--minSpotId {start}' if start else ''} \
                {f'--maxSpotId {end}' if end else ''}

            mv {snakemake.params.downloaded_fastqs_folder}/{srr}_1.fastq.gz {snakemake.output.fastq1}
            mv {snakemake.params.downloaded_fastqs_folder}/{srr}_2.fastq.gz {snakemake.output.fastq2}
            """
        )
