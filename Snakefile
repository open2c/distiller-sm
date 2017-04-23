from urllib.parse import urlparse
from _distiller_common import organize_fastqs, needs_downloading

configfile: "config.yml"

workdir: config['project_folder']

LIBRARY_RUN_FASTQS = organize_fastqs(config)


rule default:
    input: 
        expand(
            "pairs/libraries/{library}.nodups.pairs.gz", 
            library=LIBRARY_RUN_FASTQS.keys())


rule download_fastqs:
    params:
        library=lambda wc: wc.library,
        run=lambda wc: wc.run,
    output:
        fastq1='downloaded_fastqs/{library}.{run}.1.fastq.gz',
        fastq2='downloaded_fastqs/{library}.{run}.2.fastq.gz',
    run:
        fastq_files = LIBRARY_RUN_FASTQS[params.library][params.run]
        if (len(fastq_files) == 1) and (fastq_files[0].startswith('sra:')):
            parsed = urlparse(fastq_files[0])
            srr, query = parsed.path, parsed.query
            start, end = 0, None
            if query:
                for kv_pair in query.split('&'):
                    k,v = kv_pair.split('=')
                    if k == 'start':
                        start = v
                    if k == 'end':
                        end = v

            shell(
                ('fastq-dump --origfmt --split-files --gzip '
                 '-O downloaded_fastqs {srr}').format(srr=srr)
                + (' --minSpotId {}'.format(start) if start else '')
                + (' --maxSpotId {}'.format(end) if end else '')
                )
            shell(
                ('mv downloaded_fastqs/{srr}_1.fastq.gz '
                 'downloaded_fastqs/{library}.{run}.1.fastq.gz').format(
                     srr=srr, library=params.library, run=params.run))
            shell(
                ('mv downloaded_fastqs/{srr}_2.fastq.gz '
                 'downloaded_fastqs/{library}.{run}.2.fastq.gz').format(
                     srr=srr, library=params.library, run=params.run))

        
rule chunk_runs:
    input:
        fastq1=lambda wc: (
            'downloaded_fastqs/{}.{}.1.fastq.gz'.format(wc.library, wc.run)
            if needs_downloading(LIBRARY_RUN_FASTQS[wc.library][wc.run], 0)
            else LIBRARY_RUN_FASTQS[wc.library][wc.run][0]),
        fastq2=lambda wc: (
            'downloaded_fastqs/{}.{}.2.fastq.gz'.format(wc.library, wc.run)
            if needs_downloading(LIBRARY_RUN_FASTQS[wc.library][wc.run], 1)
            else LIBRARY_RUN_FASTQS[wc.library][wc.run][1]),
    params:
        chunksize=expand("{chunksize}", chunksize=4*config['chunksize']),
        library=lambda wc: wc.library,
        run=lambda wc: wc.run,
    output:
        chunk1=temp(dynamic('fastq_chunks/{library}.{run}.{chunk_id}.1.fastq.gz')),
        chunk2=temp(dynamic('fastq_chunks/{library}.{run}.{chunk_id}.2.fastq.gz')),
    shell:
        """
        zcat {input.fastq1} | split -l {params.chunksize} -d \
            --filter 'gzip > $FILE.1.fastq.gz' - \
            fastq_chunks/{params.library}.{params.run}.
        zcat {input.fastq2} | split -l {params.chunksize} -d \
            --filter 'gzip > $FILE.2.fastq.gz' - \
            fastq_chunks/{params.library}.{params.run}.
        """


rule map_chunks:
    input:
        fastq1='fastq_chunks/{library}.{run}.{chunk_id}.1.fastq.gz',
        fastq2='fastq_chunks/{library}.{run}.{chunk_id}.2.fastq.gz',
        index_fasta=expand('{index}', index=config['genome']['fasta_path']),
        index_other=expand('{index}.{res}', 
                           index=config['genome']['fasta_path'],
                           res=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        "sam/chunks/{library}.{run}.{chunk_id}.bam"
    shell: 
        """
        bwa mem -SP {input.index_fasta} {input.fastq1} {input.fastq2} \
            | samtools view -bS > {output}
        """


rule map_runs:
    input:
        fastq1=lambda wc: (
            'downloaded_fastqs/{}.{}.1.fastq.gz'.format(wc.library, wc.run)
            if needs_downloading(LIBRARY_RUN_FASTQS[wc.library][wc.run], 0)
            else LIBRARY_RUN_FASTQS[wc.library][wc.run][0]),
        fastq2=lambda wc: (
            'downloaded_fastqs/{}.{}.2.fastq.gz'.format(wc.library, wc.run)
            if needs_downloading(LIBRARY_RUN_FASTQS[wc.library][wc.run], 1)
            else LIBRARY_RUN_FASTQS[wc.library][wc.run][1]),
        index_fasta=expand('{index}', index=config['genome']['fasta_path']),
        index_other=expand('{index}.{res}', 
                           index=config['genome']['fasta_path'],
                           res=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        "sam/runs/{library}.{run}.bam"
    shell: 
        "bwa mem -SP {input.index_fasta} {input.fastq1} {input.fastq2} | samtools view -bS > {output}"


rule parse_runs:
    input:
        (dynamic("sam/chunks/{library}.{run}.{chunk_id}.bam")
         if config.get('chunksize', 0) 
         else "sam/runs/{library}.{run}.bam")
    params:
        dropsam_flag='--drop-sam' if config.get('drop_sam',False) else '',
        dropreadid_flag='--drop-readid' if config.get('drop_readid', False) else '',
        assembly=config['genome']['assembly']
    output:
        "pairsam/runs/{library}.{run}.pairsam.gz"
    run: 
        if config.get('chunksize', 0):
            shell(
                """
                cat <( samtools merge - {input} | samtools view -H ) \
                    <( samtools cat {input} | samtools view ) \
                    | pairsamtools parse {params.dropsam_flag} {params.dropreadid_flag} \
                    | pairsamtools sort -o {output}
                """
            )
        else:
            shell("""
                pairsamtools parse {input} {params.dropsam_flag} {params.dropreadid_flag} \
                | pairsamtools sort -o {output}
            """)


rule make_run_stats:
    input:
        "pairsam/runs/{library}.{run}.pairsam.gz"
    output:
        "stats/runs/{library}.{run}.stats.tsv"
    shell:
        "pairsamtools stats {input} -o {output}"


rule merge_runs_into_libraries:
    input:
        pairsams=lambda wc: expand(
            "pairsam/runs/{library}.{run}.pairsam.gz", 
            library=wc.library,
            run=LIBRARY_RUN_FASTQS[wc.library].keys(),
        ),
        stats=lambda wc: expand(
            "stats/runs/{library}.{run}.stats.tsv", 
            library=wc.library,
            run=LIBRARY_RUN_FASTQS[wc.library].keys(),
        ),
    output:
        pairsam="pairsam/libraries/{library}.pairsam.gz",
        stats="stats/libraries/{library}.stats.tsv",
    shell:
        """
        pairsamtools merge {input.pairsams} -o {output.pairsam}
        pairsamtools stats --merge {input.stats} -o {output.stats}
        """


rule make_pairs_bams:
    input:
        pairsam="pairsam/libraries/{library}.pairsam.gz"

    output:
        stats="stats/libraries/{library}.dedup.stats.tsv",
        nodups_pairs=   "pairs/libraries/{library}.nodups.pairs.gz",
        unmapped_pairs= "pairs/libraries/{library}.unmapped.pairs.gz",
        dups_pairs=     "pairs/libraries/{library}.dups.pairs.gz",
        nodups_sam=   "sam/libraries/{library}.nodups.bam",
        unmapped_sam= "sam/libraries/{library}.unmapped.bam",
        dups_sam=     "sam/libraries/{library}.dups.bam",
    run:
        if config.get('drop_sam', False):
            shell("""
            pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
                {input} \
                --output-rest >( pairsamtools split \
                    --output-pairs {output.unmapped_pairs} \
                    ) | \
            pairsamtools dedup \
                --output \
                    >( pairsamtools split \
                        --output-pairs {output.nodups_pairs} \
                     ) \
                --output-dups \
                    >( pairsamtools markasdup \
                        | pairsamtools split \
                            --output-pairs {output.dups_pairs} \
                     ) \
                --stats-file {output.stats}
            touch {output.dups_sam}
            touch {output.unmapped_sam}
            touch {output.nodups_sam}
            """)

        else:
            shell("""
            pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
                {input} \
                --output-rest >( pairsamtools split \
                    --output-pairs {output.unmapped_pairs} \
                    --output-sam {output.unmapped_sam} \
                    ) | \
            pairsamtools dedup \
                --output \
                    >( pairsamtools split \
                        --output-pairs {output.nodups_pairs} \
                        --output-sam {output.nodups_sam} \
                     ) \
                --output-dups \
                    >( pairsamtools markasdup \
                        | pairsamtools split \
                            --output-pairs {output.dups_pairs} \
                            --output-sam {output.dups_sam} \
                     ) \
                --stats-file {output.stats}
            """)


rule index_pairs:
    input:
        "pairs/libraries/{library}.nodups.pairs.gz"
    output:
        "pairs/libraries/{library}.nodups.pairs.gz.px2"
    shell: 
        "pairix {input}"


rule make_library_coolers:
    input:
        pairs="pairs/libraries/{library}.nodups.pairs.gz",
        pairix_index="pairs/libraries/{library}.nodups.pairs.gz.px2",
        chrom_sizes=expand(
            '{chrom_sizes}', chrom_sizes=config['genome']['chrom_sizes_path']),
    params:
        res=lambda wc: wc.res,
        assembly=expand("{assembly}", assembly=config['genome']['assembly'])
    output:
        "coolers/libraries/{library}.{res}.cool"
    shell:
        """
        cooler cload pairix \
            --assembly {params.assembly} \
            {input.chrom_sizes}:{params.res} {input.pairs} {output}
        """


rule merge_library_group_stats:
    input:
        stats=lambda wc: expand(
            "stats/libraries/{library}.stats.tsv",
            library=config['library_groups'][wc.library_group],
            ),
        stats_dedup=lambda wc: expand(
            "stats/libraries/{library}.dedup.stats.tsv",
            library=config['library_groups'][wc.library_group],
            )

    output:
        stats="stats/library_groups/{library_group}.stats.tsv"
    shell:
        """
        pairsamtools stats --merge {input.stats} {input.stats_dedup} -o {output.stats}
        """


rule make_library_group_coolers:
    input:
        coolers=lambda wc: expand(
            "coolers/libraries/{library}.{res}.cool", 
            library=config['library_groups'][wc.library_group],
            res=wc.res),
    output:
        cooler="coolers/library_groups/{library_group}.{res}.cool",
    shell:
        """
        cooler merge {output.cooler} {input.coolers}
        """


rule make_all_coolers:
    input: 
        library_coolers = expand(
            "coolers/libraries/{library}.{res}.cool",
            library=LIBRARY_RUN_FASTQS.keys(),
            res=config['cooler_resolutions']),
        library_group_coolers = expand(
            "coolers/library_groups/{library_group}.{res}.cool",
            library_group=config['library_groups'].keys(), 
            res=config['cooler_resolutions']),
        library_group_stats = expand(
            "stats/library_groups/{library_group}.stats.tsv",
            library_group=config['library_groups'].keys(), 
            )

        
