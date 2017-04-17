from _distiller_common import organize_fastqs

configfile: "config.yaml"

workdir: config['project_folder']


RUN_TO_FASTQS, RUN_FULL_NAMES, LIBRARY_TO_FASTQS, LIBRARY_NAMES = organize_fastqs(
    config)


rule all:
    input: expand(
        "libraries/pairs/{library}.nodups.pairs.gz.px2", 
        library=LIBRARY_NAMES)


rule map:
    input:
        fastq1=lambda wildcards: RUN_TO_FASTQS[wildcards.chunk][0],
        fastq2=lambda wildcards: RUN_TO_FASTQS[wildcards.chunk][1],
        index_fasta=expand('{index}', index=config['genome']['fasta_path']),
        index_other=expand('{index}.{res}', 
                     index=config['genome']['fasta_path'],
                     res=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        "chunks/sam/{chunk}.bam"
    shell: 
        "bwa mem -SP {input.index_fasta} {input.fastq1} {input.fastq2} | samtools view -bS > {output}"


rule parse:
    input:
        "chunks/sam/{chunk}.bam"
    output:
        "chunks/pairsam/{chunk}.pairsam.gz"
    shell: 
        "pairsamtools parse {input} | pairsamtools sort -o {output}"


rule merge_pairsam:
    input:
        lambda wildcards: expand("chunks/pairsam/{chunk}.pairsam.gz", 
                                 chunk=LIBRARY_TO_FASTQS[wildcards.library])
    output:
        "libraries/pairsam/{library}.pairsam.gz"
    shell:
        "pairsamtools merge {input} -o {output}"


rule make_pairs_bams:
    input:
        "libraries/pairsam/{library}.pairsam.gz"
    output:
        "libraries/pairs/{library}.nodups.pairs.gz"
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
                 )
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
        library="libraries/pairs/{library}.nodups.pairs.gz",
        library_index="libraries/pairs/{library}.nodups.pairs.gz.px2",
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
            {input.chrom_sizes}:{params.res} {input.library} {output}
        """
       
rule test_coolers:
    input: 
        expand("libraries/coolers/{library}.{res}.cool",
               library=LIBRARY_NAMES, res=config['cooler_resolutions'])
