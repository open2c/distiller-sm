from _distiller_common import flatten_tree

configfile: "config.yaml"

CHUNK_FASTQS = {config['exp_name_sep'].join(k):v 
                for k,v in flatten_tree(config['experiments']).items()}

CHUNK_NAMES = list(CHUNK_FASTQS.keys())

EXPERIMENT_CHUNKS = {exp:[config['exp_name_sep'].join(k)
                          for k in flatten_tree({exp:config['experiments'][exp]})]
                     for exp in config['experiments'] }
EXPERIMENT_NAMES = list(EXPERIMENT_CHUNKS.keys())

rule all:
    input: expand("exps/pairs/{experiment}.nodups.pairs.gz", experiment=EXPERIMENT_NAMES)

rule map:
    params:
        fastq1=lambda wildcards: CHUNK_FASTQS[wildcards.chunk][0],
        fastq2=lambda wildcards: CHUNK_FASTQS[wildcards.chunk][1],
        index=config['genome']['fasta_path']
    output:
        "chunks/sam/{chunk}.bam"
    shell: 
        'touch {output};'
        'echo "bwa mem -SP {params.index} {params.fastq1} {params.fastq2}"'
        '" | samtools view -bS > {output}"'

rule parse:
    input:
        "chunks/sam/{chunk}.bam"
    output:
        "chunks/pairsam/{chunk}.pairsam.gz"
    shell: 
        'touch {output};'
        'echo "pairsamtools parse {input} | pairsamtools sort -o {output}"'

rule merge:
    input:
        lambda wildcards: expand("chunks/pairsam/{chunk}.pairsam.gz", 
                                 chunk=EXPERIMENT_CHUNKS[wildcards.experiment])
    output:
        "exps/pairsam/{experiment}.pairsam.gz"
    shell:
        'touch {output};'
        'echo "pairsamtools merge {input} -o {output}"'

rule make_pairs:
    output:
        "exps/pairs/{experiment}.nodups.pairs.gz"
    shell: """
        pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
            --output-rest >( pairsamtools split \
                --output-pairs exps/pairs/{wildcard.experiment}.unmapped.pairs.gz \
                --output-sam exps/sam/{wildcard.experiment}.unmapped.bam \

        pairsamtools dedup \
            --output \
                >( pairsamtools split \
                    --output-pairs exps/pairs/{wildcard.experiment}.nodups.pairs.gz \
                    --output-sam exps/sam/{wildcard.experiment}.nodups.bam \
            --output-dups \
                >( pairsamtools markasdup \
                    | pairsamtools split \
                        --output-pairs exps/pairs/{wildcard.experiment}.dups.pairs.gz \
                        --output-sam exps/sam/{wildcard.experiment}.dups.bam \
        """
