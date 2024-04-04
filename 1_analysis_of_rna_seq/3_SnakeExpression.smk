#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_ALL
INDIR = "results/mapping/rmdup"
OUTDIR = "results/expression"

rule all:
    input:
        expand(OUTDIR + "/stringtie/{sample}", sample=SAMPLES),

rule stringtie:
    input:
        bam = INDIR + "/{sample}.bam",
        gtf = ANNOTATION_GTF
    output:
        out = directory(OUTDIR + "/stringtie/{sample}")
    log:
        OUTDIR + "/stringtie/{sample}.log"
    conda:
        "stringtie"
    threads:
        4
    shell:
        """
        mkdir {output.out}
        stringtie {input.bam} \
            --rf \
            -e \
            -A {output.out}/gene_abund.tab \
            -C {output.out}/cov_refs.gtf \
            -p {threads} \
            -G {input.gtf} \
            -o {output.out}/transcripts.gtf &> {log}
        """
