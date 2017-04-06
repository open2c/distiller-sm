from _distiller_common import organize_fastqs

configfile: "config.yaml"

workdir: config['project_folder']

RUN_TO_FASTQS, RUN_FULL_NAMES, EXPERIMENT_TO_FASTQS, EXPERIMENT_NAMES = organize_fastqs(
    config)


rule all:
    input: expand("exps/pairs/{experiment}.nodups.pairs.gz", experiment=EXPERIMENT_NAMES)


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


rule merge:
    input:
        lambda wildcards: expand("chunks/pairsam/{chunk}.pairsam.gz", 
                                 chunk=EXPERIMENT_TO_FASTQS[wildcards.experiment])
    output:
        "exps/pairsam/{experiment}.pairsam.gz"
    shell:
        "pairsamtools merge {input} -o {output}"


rule make_pairs_bams:
    input:
        "exps/pairsam/{experiment}.pairsam.gz"
    output:
        "exps/pairs/{experiment}.nodups.pairs.gz"
    shell: """
        mkdir -p exps/sam ;
        pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
            {input} \
            --output-rest >( pairsamtools split \
                --output-pairs exps/pairs/{wildcards.experiment}.unmapped.pairs.gz \
                --output-sam exps/sam/{wildcards.experiment}.unmapped.bam \
                ) | \
        pairsamtools dedup \
            --output \
                >( pairsamtools split \
                    --output-pairs exps/pairs/{wildcards.experiment}.nodups.pairs.gz \
                    --output-sam exps/sam/{wildcards.experiment}.nodups.bam \
                 ) \
            --output-dups \
                >( pairsamtools markasdup \
                    | pairsamtools split \
                        --output-pairs exps/pairs/{wildcards.experiment}.dups.pairs.gz \
                        --output-sam exps/sam/{wildcards.experiment}.dups.bam \
                 )
        """

