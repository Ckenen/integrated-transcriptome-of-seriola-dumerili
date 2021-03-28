#!/usr/bin/env snakemake
names = [
    "stringtie_stringtie", 
    "stringtie_taco",
    "cufflink_cuffmerge",
    "cufflink_taco"
]

rule all:
    input:
        expand("results/gffcompare/{name}.tracking", name=names),

rule gffcompare:
    input:
        ref = "data/genome/annotation.gtf",
        gtf = "results/collect/{name}.gtf"
    output:
        out1 = "results/gffcompare/{name}.loci",
        out2 = "results/gffcompare/{name}.tracking",
        out3 = "results/gffcompare/{name}.annotated.gtf",
        out4 = "results/gffcompare/{name}.stats",
    log:
        log = "results/gffcompare/{name}.log"
    params:
        prefix = "results/gffcompare/{name}"
    shell:
        """
        gffcompare -r {input.ref} -R -o {params.prefix} {input.gtf} &> {log}
        """