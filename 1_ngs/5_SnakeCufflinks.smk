#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/cufflinks"


rule all:
    input:
        expand("results/cufflinks/clout/{sample}", sample=samples),
        # "results/cufflinks/merged",


rule cufflinks:
    input:
        bam = indir + "/{sample}.bam"
    output:
        directory("results/cufflinks/clout/{sample}")
    log:
        "results/cufflinks/clout/{sample}.log"
    threads:
        threads
    shell:
        """
        cufflinks -p {threads} -o {output} {input.bam} &> {log}
        """

rule cuffmerge:
    input:
        gtf = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta",
        lis = ["results/cufflinks/clout/{sample}".format(sample=sample) for sample in samples]
    output:
        txt = "results/cufflinks/merged.txt",
        out = directory("results/cufflinks/merged")
    log:
        "results/cufflinks/merged.log"
    shell:
        """
        for f in {input.lis}; do echo $f; done > {output.txt}
        cuffmerge -g {input.gtf} -s {input.fsa} -p {threads} {output.txt} -o {output.out} &> {log}
        """
        