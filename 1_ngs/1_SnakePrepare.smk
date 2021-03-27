#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "data/datasets"
outdir = "results/prepare"
reads = ["R1", "R2"]

rule all:
    input:
        expand(outdir + "/fastqc/{sample}_{read}_fastqc.html", sample=samples, read=reads),

# FastQC

rule fastqc:
    input:
        indir + "/{name}.fastq.gz"
    output:
        outdir + "/fastqc/{name}_fastqc.html"
    log:
        outdir + "/fastqc/{name}.log"
    params:
        odir = outdir + "/fastqc"
    shell:
        """
        fastqc -o {params.odir} {input} &> {log}
        """