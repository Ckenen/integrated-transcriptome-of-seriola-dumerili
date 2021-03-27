#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/cufflinks"


rule all:
    input:
        expand(outdir + "/clout/{sample}", sample=samples),
        "results/cufflinks/merged",
        outdir + "/taco/merged"


# Cufflinks

rule cufflinks:
    input:
        bam = indir + "/{sample}.bam",
        gtf = "data/genome/annotation.gtf"
    output:
        directory(outdir + "/clout/{sample}")
    log:
        outdir + "/clout/{sample}.log"
    threads:
        threads
    shell:
        """
        cufflinks -p {threads} -o {output} -g {input.gtf} --library-type fr-firststrand \
            -L {wildcards.sample} {input.bam} &> {log}
        """

rule cuffmerge:
    input:
        gtf = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta",
        lis = [outdir + "/clout/{sample}".format(sample=sample) for sample in samples]
    output:
        txt = outdir + "/merged.gtf.txt",
        out = directory(outdir + "/merged")
    log:
        outdir + "/merged.log"
    threads:
        threads
    shell:
        """
        for f in {input.lis}; do echo "$f/transcripts.gtf"; done > {output.txt}
        cuffmerge -g {input.gtf} -s {input.fsa} -o {output.out} -p {threads} {output.txt}  &> {log}
        """
        
# TACO

rule taco_merge:
    input:
        gtfs = expand(rules.cufflinks.output, sample=samples)
    output:
        txt = outdir + "/taco/merged.gtf.txt",
        out = directory(outdir + "/taco/merged")
    log:
        outdir + "/taco/merged.log"
    params:
        gtf = outdir + "/taco/merged/assembly.gtf",
        bed = outdir + "/taco/merged/assembly.bed"
    threads:
        threads
    shell:
        """
        for f in {input.gtfs}; do echo "$f/transcripts.gtf"; done > {output.txt}
        taco_run -p {threads} -o {output.out} --gtf-expr-attr FPKM --filter-min-length 200 \
            --filter-min-expr 0.5 --isoform-frac 0.05 {output.txt} &> {log}
        """