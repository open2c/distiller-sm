from _distiller_common import organize_fastqs

configfile: "config.yaml"

workdir: config['project_folder']


RUN_TO_FASTQS, RUN_FULL_NAMES, LIBRARY_TO_FASTQS, LIBRARY_NAMES = organize_fastqs(
    config)


rule default:
    input: expand("libraries/pairs/{library}.nodups.pairs.gz", library=LIBRARY_NAMES)


rule chunk_runs:
    input:
        fastq1=lambda wildcards: RUN_TO_FASTQS[wildcards.run][0],
        fastq2=lambda wildcards: RUN_TO_FASTQS[wildcards.run][1],
    params:
        chunksize=expand("{chunksize}", chunksize=4*config['chunksize']),
        run=lambda wildcards: wildcards.run,
    output:
        chunk1=dynamic('chunks/fastq/{run}/{run}.1.fastq.chunk.{chunk_id}.gz'),
        chunk2=dynamic('chunks/fastq/{run}/{run}.2.fastq.chunk.{chunk_id}.gz'),
    shell:
        """
        zcat {input.fastq1} | split -l {params.chunksize} -d \
            --filter 'gzip > $FILE.gz' - chunks/fastq/{params.run}/{params.run}.1.fastq.chunk.
        zcat {input.fastq2} | split -l {params.chunksize} -d \
            --filter 'gzip > $FILE.gz' - chunks/fastq/{params.run}/{params.run}.2.fastq.chunk.
        """


rule map_chunks:
    input:
        fastq1='chunks/fastq/{run}/{run}.1.fastq.chunk.{chunk_id}.gz',
        fastq2='chunks/fastq/{run}/{run}.2.fastq.chunk.{chunk_id}.gz',
        index_fasta=expand('{index}', index=config['genome']['fasta_path']),
        index_other=expand('{index}.{res}', 
                           index=config['genome']['fasta_path'],
                           res=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        "chunks/sam/{run}/{run}.chunk.{chunk_id}.bam"
    shell: 
        """
        bwa mem -SP {input.index_fasta} {input.fastq1} {input.fastq2} \
            | samtools view -bS > {output}
        """

rule map_runs:
    input:
        fastq1=lambda wildcards: RUN_TO_FASTQS[wildcards.run][0],
        fastq2=lambda wildcards: RUN_TO_FASTQS[wildcards.run][1],
        index_fasta=expand('{index}', index=config['genome']['fasta_path']),
        index_other=expand('{index}.{res}', 
                           index=config['genome']['fasta_path'],
                           res=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        "runs/sam/{run}.bam"
    shell: 
        "bwa mem -SP {input.index_fasta} {input.fastq1} {input.fastq2} | samtools view -bS > {output}"


rule parse_runs:
    input:
        dynamic("chunks/sam/{run}/{run}.chunk.{chunk_id}.bam") if config.get('chunksize', 0) else "runs/sam/{run}.bam" 
    output:
        "runs/pairsam/{run}.pairsam.gz"
    run: 
        if config.get('chunksize', 0):
            shell(
                """
                cat <( samtools merge - {input} | samtools view -H ) \
                    <( samtools cat {input} | samtools view ) \
                    | pairsamtools parse | pairsamtools sort -o {output}
                """
            )
        else:
            shell("pairsamtools parse {input} | pairsamtools sort -o {output}")


rule make_run_stats:
    input:
        "runs/pairsam/{run}.pairsam.gz"
    output:
        "runs/stats/{run}.stats.tsv"
    shell:
        "pairsamtools stats {input} -o {output}"


rule merge_runs_into_libraries:
    input:
        pairsams=lambda wildcards: expand("runs/pairsam/{run}.pairsam.gz", 
                                 run=LIBRARY_TO_FASTQS[wildcards.library]),
        stats=lambda wildcards: expand("runs/stats/{run}.stats.tsv", 
                                 run=LIBRARY_TO_FASTQS[wildcards.library]),
    output:
        pairsam="libraries/pairsam/{library}.pairsam.gz",
        stats="libraries/stats/{library}.stats.tsv",
    shell:
        """
        pairsamtools merge {input.pairsams} -o {output.pairsam}
        pairsamtools stats --merge {input.stats} -o {output.stats}
        """


rule make_pairs_bams:
    input:
        pairsam="libraries/pairsam/{library}.pairsam.gz"

    output:
        pairs="libraries/pairs/{library}.nodups.pairs.gz",
        stats="libraries/stats/{library}.dedup.stats.tsv"
    shell: """
        mkdir -p libraries/sam ;
        pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
            {input} \
            --output-rest >( pairsamtools split \
                --output-pairs libraries/pairs/{wildcards.library}.unmapped.pairs.gz \
                --output-sam libraries/sam/{wildcards.library}.unmapped.bam \
                ) | \
        pairsamtools dedup \
            --output \
                >( pairsamtools split \
                    --output-pairs libraries/pairs/{wildcards.library}.nodups.pairs.gz \
                    --output-sam libraries/sam/{wildcards.library}.nodups.bam \
                 ) \
            --output-dups \
                >( pairsamtools markasdup \
                    | pairsamtools split \
                        --output-pairs libraries/pairs/{wildcards.library}.dups.pairs.gz \
                        --output-sam libraries/sam/{wildcards.library}.dups.bam \
                 ) \
            --stats-file {output.stats}
        """


rule index_pairs:
    input:
        "libraries/pairs/{library}.nodups.pairs.gz"
    output:
        "libraries/pairs/{library}.nodups.pairs.gz.px2"
    shell: 
        "pairix {input}"


rule make_library_coolers:
    input:
        pairs="libraries/pairs/{library}.nodups.pairs.gz",
        pairix_index="libraries/pairs/{library}.nodups.pairs.gz.px2",
        chrom_sizes=expand(
            '{chrom_sizes}', chrom_sizes=config['genome']['chrom_sizes_path']),
    params:
        res=lambda wildcards: wildcards['res'],
        assembly=expand("{assembly}", assembly=config['genome']['name'])
    output:
        "libraries/coolers/{library}.{res}.cool"
    shell:
        """
        cooler cload pairix \
            --assembly {params.assembly} \
            {input.chrom_sizes}:{params.res} {input.pairs} {output}
        """


rule merge_library_group_stats:
    input:
        stats=lambda wildcard: expand(
            "libraries/stats/{library}.stats.tsv",
            library=config['library_groups'][wildcard.library_group],
            ),
        stats_dedup=lambda wildcard: expand(
            "libraries/stats/{library}.dedup.stats.tsv",
            library=config['library_groups'][wildcard.library_group],
            )

    output:
        stats="library_groups/stats/{library_group}.stats.tsv"
    shell:
        """
        pairsamtools stats --merge {input.stats} {input.stats_dedup} -o {output.stats}
        """


rule make_library_group_coolers:
    input:
        coolers=lambda wildcard: expand(
            "libraries/coolers/{library}.{res}.cool", 
            library=config['library_groups'][wildcard.library_group],
            res=wildcard.res),
    output:
        cooler="library_groups/coolers/{library_group}.{res}.cool",
    shell:
        """
        cooler merge {output.cooler} {input.coolers}
        """


rule make_all_coolers:
    input: 
        library_coolers = expand(
            "libraries/coolers/{library}.{res}.cool",
            library=LIBRARY_NAMES,
            res=config['cooler_resolutions']),
        library_group_coolers = expand(
            "library_groups/coolers/{library_group}.{res}.cool",
            library_group=config['library_groups'].keys(), 
            res=config['cooler_resolutions']),
        library_group_stats = expand(
            "library_groups/stats/{library_group}.stats.tsv",
            library_group=config['library_groups'].keys(), 
            )

        
