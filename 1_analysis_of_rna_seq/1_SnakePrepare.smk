#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_ALL
INDIR = "data/datasets"
OUTDIR = "results/prepare"
RS = ["R1", "R2"]


rule all:
    input:
        expand(OUTDIR + "/fastqc/{sample}_{r}_fastqc.html", sample=SAMPLES, r=RS),

# FastQC

rule fastqc:
    input:
        fq = INDIR + "/{name}.fastq.gz"
    output:
        html = OUTDIR + "/fastqc/{name}_fastqc.html"
    log:
        log = OUTDIR + "/fastqc/{name}.log"
    shell:
        """
        fastqc -o `dirname {output.html}` {input} &> {log}
        """
        