#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/rmdup"
outdir = "results/expression"

rule all:
    input:
        expand(outdir + "/stringtie/{sample}", sample=samples),

rule stringtie:
    input:
        bam = indir + "/{sample}.bam",
        gtf = ANNOTATION_GTF
    output:
        out = directory(outdir + "/stringtie/{sample}")
    log:
        outdir + "/stringtie/{sample}.log"
    threads:
        8
    shell:
        """
        mkdir {output.out}
        stringtie {input.bam} --rf -e -A {output.out}/gene_abund.tab \
            -C {output.out}/cov_refs.gtf -p {threads} -G {input.gtf} \
            -o {output.out}/transcripts.gtf &> {log}
        """
