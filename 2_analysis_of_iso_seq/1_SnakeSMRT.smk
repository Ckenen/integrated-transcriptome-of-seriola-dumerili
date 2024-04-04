#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "data/datasets"
OUTDIR = "results/smrt"

PRIMER = "data/primers.fasta"

rule all:
    input:
        expand(OUTDIR + "/ccs/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/demux/{sample}.primer_5p--primer_3p.bam", sample=SAMPLES),
        expand(OUTDIR + "/flnc/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/clustered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/polished/{sample}.bam", sample=SAMPLES),

# Generate circular consensus sequences (ccs) from subreads.
# time-consuming

rule ccs:
    input:
        bam = INDIR + "/{sample}.subreads.bam",
    output:
        bam = OUTDIR + "/ccs/{sample}.bam"
    log:
        OUTDIR + "/ccs/{sample}.log"
    conda:
        "smrt"
    threads:
        THREADS
    shell:
        """
        ccs --skip-polish \
            --min-passes {MIN_PASSES} \
            --min-length 200 \
            --min-rq 0.99 \
            --num-threads {threads} \
            {input.bam} {output.bam} &> {log}
        """

# Lima, Demultiplex Barcoded PacBio Data and Clip Barcodes

rule lima:
    input:
        bam = rules.ccs.output.bam,
        primer = PRIMER
    output:
        bam = OUTDIR + "/demux/{sample}.primer_5p--primer_3p.bam"
    log:
        OUTDIR + "/demux/{sample}.log"
    conda:
        "smrt"
    params:
        bam = OUTDIR + "/demux/{sample}.bam"
    threads:
        THREADS
    shell:
        """
        lima --isoseq -j {threads} {input.bam} \
            {input.primer} {params.bam} &> {log}
        """

# Remove polyA and concatemers from FL reads and generate FLNC transcripts (FL to FLNC)
# full-length, non-concatemer (FLNC) reads

rule refine:
    input:
        bam = rules.lima.output.bam,
        primer = PRIMER
    output:
        bam = OUTDIR + "/flnc/{sample}.bam"
    log:
        OUTDIR + "/flnc/{sample}.log"
    conda:
        "smrt"
    threads:
        THREADS
    shell:
        """
        isoseq3 refine --require-polya --min-polya-length 20 \
            -j {threads} {input.bam} {input.primer} {output.bam} &> {log}
        """

# Cluster FLNC reads and generate unpolished transcripts (FLNC to UNPOLISHED)
# time-consuming

rule cluster:
    input:
        bam = rules.refine.output.bam
    output:
        bam = OUTDIR + "/clustered/{sample}.bam"
    log:
        OUTDIR + "/clustered/{sample}.log"
    conda:
        "smrt"
    threads:
        THREADS
    shell:
        """
        isoseq3 cluster -j {threads} {input.bam} {output.bam} &> {log}
        """

# Polish transcripts using subreads (UNPOLISHED to POLISHED)
# time-consuming

rule polish:
    input:
        bam1 = rules.cluster.output.bam,
        bam2 = INDIR + "/{sample}.subreads.bam"
    output:
        bam = OUTDIR + "/polished/{sample}.bam"
    log:
        OUTDIR + "/polished/{sample}.log"
    conda:
        "smrt"
    threads:
        THREADS
    shell:
        """
        isoseq3 polish -j {threads} \
            {input.bam1} {input.bam2} {output.bam} &> {log}
        """
