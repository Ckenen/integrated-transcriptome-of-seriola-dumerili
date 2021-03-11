#!/usr/bin/env python
configfile: "config.yaml"
samples = config["samples"]
threads = 20
indir = "data/datasets"
outdir = "results/isoseq"


rule all:
    input:
        expand(outdir + "/ccs/{sample}.bam", sample=samples),
        expand(outdir + "/demux/{sample}.primer_5p--primer_3p.bam", sample=samples),
        expand(outdir + "/flnc/{sample}.bam", sample=samples),
        expand(outdir + "/clustered/{sample}.bam", sample=samples),
        expand(outdir + "/polished/{sample}.bam", sample=samples),
        outdir + "/genome.isoseq.mmi",
        expand(outdir + "/aligned/{sample}.bam", sample=samples),
        expand(outdir + "/collapsed/{sample}.gff", sample=samples),
        expand(outdir + "/collapsed/{sample}.sorted.gtf.gz", sample=samples),
        # expand(outdir + "/collapsed/{sample}.bed.gz", sample=samples),

# Generate circular consensus sequences (ccs) from subreads.

rule ccs:
    input:
        indir + "/{sample}.bam",
    output:
        outdir + "/ccs/{sample}.bam"
    log:
        outdir + "/ccs/{sample}.log"
    threads:
        threads
    shell:
        """
        ccs --skip-polish --min-passes 1 --min-rq 0.99 -j {threads} {input} {output} &> {log}
        """

# Lima, Demultiplex Barcoded PacBio Data and Clip Barcodes

rule lima:
    input:
        bam = rules.ccs.output,
        primer = "data/primers.fasta"
    output:
        outdir + "/demux/{sample}.primer_5p--primer_3p.bam"
    log:
        outdir + "/demux/{sample}.log"
    params:
        bam = outdir + "/demux/{sample}.bam"
    threads:
        threads
    shell:
        """
        lima --isoseq -j {threads} {input} {params.bam} &> {log}
        """

# Remove polyA and concatemers from FL reads and generate FLNC transcripts (FL to FLNC)
# full-length, non-concatemer (FLNC) reads

rule refine:
    input:
        bam = rules.lima.output,
        primer = "data/primers.fasta"
    output:
        bam = outdir + "/flnc/{sample}.bam"
    log:
        outdir + "/flnc/{sample}.log"
    threads:
        threads
    shell:
        """
        isoseq3 refine --require-polya --min-polya-length 20 -j {threads} {input} {output} &> {log}
        """

# Cluster FLNC reads and generate unpolished transcripts (FLNC to UNPOLISHED)

rule cluster:
    input:
        rules.refine.output.bam
    output:
        outdir + "/clustered/{sample}.bam"
    log:
        outdir + "/clustered/{sample}.log"
    threads:
        threads
    shell:
        """
        isoseq3 cluster -j {threads} {input} {output} &> {log}
        """

# Polish transcripts using subreads (UNPOLISHED to POLISHED)

rule polish:
    input:
        bam1 = rules.cluster.output,
        bam2 = indir + "/{sample}.bam",
        pbi2 = indir + "/{sample}.bam.pbi",
    output:
        outdir + "/polished/{sample}.bam"
    log:
        outdir + "/polished/{sample}.log"
    threads:
        threads
    shell:
        """
        isoseq3 polish -j {threads} {input.bam1} {input.bam2} {output} &> {log}
        """

# Index reference and store as .mmi file

rule mmindex_isoseq:
    input:
        "data/genome/genome.fasta"
    output:
        outdir + "/genome.isoseq.mmi"
    threads:
        threads
    shell:
        """
        pbmm2 index -j {threads} --preset ISOSEQ {input} {output}
        """

# Align PacBio reads to reference sequences

rule align:
    input:
        rules.mmindex_isoseq.output,
        rules.polish.output
    output:
        outdir + "/aligned/{sample}.bam"
    log:
        outdir + "/aligned/{sample}.log"
    threads:
        threads
    shell:
        """
        pbmm2 align -j {threads} --preset ISOSEQ --sort {input} {output} &> {log}
        """

# Collapse transcripts based on genomic mapping

rule collapse:
    input:
        bam = rules.align.output,
        ccs = rules.ccs.output
    output:
        gff = outdir + "/collapsed/{sample}.gff",
        fa = outdir + "/collapsed/{sample}.fasta",
        fq = outdir + "/collapsed/{sample}.fastq"
    threads:
        threads
    shell:
        """
        isoseq3 collapse -j {threads} {input.bam} {input.ccs} {output.gff}
        """

rule bgzipGtf:
    input:
        gff = outdir + "/collapsed/{sample}.gff"
    output:
        gtf = outdir + "/collapsed/{sample}.sorted.gtf.gz",
        tbi = outdir + "/collapsed/{sample}.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule gtf_to_bed:
    input:
        gtf = rules.bgzipGtf.output.gtf
    output:
        bed = outdir + "/collapsed/{sample}.bed.gz",
        tbi = outdir + "/collapsed/{sample}.bed.gz.tbi"
    shell:
        """
        ./scripts/pacbio_gtf_to_bed.py {input} | bedtools sort -i - | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

# Common 

rule pbindex:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.pbi"
    shell:
        """
        pbindex {input}
        """

rule bam_stats:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """
