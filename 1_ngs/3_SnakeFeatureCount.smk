#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/expression"


rule all:
    input:
        expand(outdir + "/featureCount/{sample}.txt", sample=samples),
        outdir + "/featureCount.merged.tsv",


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

rule merge_feature_count:
    input:
        count = expand(outdir + "/featureCount/{sample}.txt", sample=samples)
    output:
        outdir + "/featureCount.merged.tsv"
    params:
        paths = ",".join([outdir + "/featureCount/{sample}.txt".format(sample=sample) for sample in samples]),
        samples = ",".join(samples),
    shell:
        """
        ./scripts/merge_feature_count.py {params.paths} {params.samples} {output}
        """
