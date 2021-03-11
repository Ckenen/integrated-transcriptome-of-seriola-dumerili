#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = config["samples"]
threads: 20


rule all:
    input:
        expand("results/tama/{sample}.bed", sample=samples),


rule TAMA:
    input:
        bam = "results/isoseq/aligned/{sample}.bam",
        bai = "results/isoseq/aligned/{sample}.bam.bai",
        fas = "data/genome/genome.fasta"
    output:
        "results/tama/{sample}.bed"
    log:
        "results/tama/{sample}.log"
    params:
        prefix = "results/tama/{sample}"
    shell:
        """
        set +u
        source activate py27
        python /home/chenzonggui/software/tama/tama_collapse.py \
            -b BAM -s {input.bam} -f {input.fas} -p {params.prefix} -x no_cap &> {log}
        conda deactivate
        """
