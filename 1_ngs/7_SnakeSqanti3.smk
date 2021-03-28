#!/usr/bin/env snakemake
threads = 8
names = [
    "stringtie_stringtie", 
    "stringtie_taco",
    "cufflink_cuffmerge",
    "cufflink_taco"
]

rule all:
    input:
        expand("results/sqanti3/{name}", name=names),

rule sqanti3:
    input:
        gtf = "results/collect/{name}.gtf",
        ref = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta"
    output:
        directory("results/sqanti3/{name}")
    log:
        log = "results/sqanti3/{name}.log"
    params:
        prefix = "results/sqanti3/{name}"
    threads:
        threads
    shell:
        """
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py -d {output} -t {threads} \
            --gtf {input.gtf} {input.ref} {input.fsa} &> {log}
        conda deactivate
        """