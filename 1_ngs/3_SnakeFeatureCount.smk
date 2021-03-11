#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/expression"


rule all:
    input:
        expand(outdir + "/featureCount/{sample}.txt", sample=samples),


rule featureCount:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
        gtf = "data/genome/annotation.gtf"
    output:
        txt = outdir + "/featureCount/{sample}.txt"
    log:
        outdir + "/featureCount/{sample}.log"
    threads:
        threads
    shell:
        """
        featureCounts -T {threads} -s 2 -p -B -a {input.gtf} -o {output.txt} {input.bam} &> {log}
        """