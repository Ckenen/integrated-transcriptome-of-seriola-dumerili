#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/prepare"
rs = ["R1", "R2"]


rule all:
    input:
        expand(outdir + "/fastqc/{sample}_{r}_fastqc.html", sample=samples, r=rs),

# FastQC

rule fastqc:
    input:
        fq = indir + "/{name}.fastq.gz"
    output:
        html = outdir + "/fastqc/{name}_fastqc.html"
    log:
        log = outdir + "/fastqc/{name}.log"
    shell:
        """
        fastqc -o `dirname {output.html}` {input} &> {log}
        """
        