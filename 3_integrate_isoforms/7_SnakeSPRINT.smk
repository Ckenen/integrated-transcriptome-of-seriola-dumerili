#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
indir = "../1_analysis_of_rna_seq/data/datasets"
outdir = "results/sprint"

rule all:
    input:
        outdir + "/prepare/ref.fasta",
        expand(outdir + "/fastqs/{sample}_{r}.fastq", sample=ngs_samples, r=["R1", "R2"]),
        expand(outdir + "/outputs/{sample}", sample=ngs_samples),
        expand(outdir + "/A2I/{sample}.tsv", sample=ngs_samples),
        expand(outdir + "/read/{sample}.txt", sample=ngs_samples),

rule prepare:
    input:
        fa = GENOME_FASTA,
        gtf = GTFS["asm"]
    output:
        fa = outdir + "/prepare/ref.fasta",
        gtf = outdir + "/prepare/ref.gtf",
    log:
        outdir + "/prepare/run.log"
    shell:
        """
        cp {input.fa} {output.fa}
        gzip -d -c {input.gtf} > {output.gtf}
        set +u; source activate SPRINT
        sprint prepare -t {output.gtf} {output.fa} `which bwa` &> {log}
        """

rule decompress_fastq:
    input:
        fq = indir + "/{name}.fastq.gz",
    output:
        fq = outdir + "/fastqs/{name}.fastq"
    shell:
        """
        gzip -d -c {input.fq} > {output.fq}
        """

rule sprint: # SPRINT_threads is a modified SPRINT version that support multi-thread samtools
    input:
        fa = rules.prepare.output.fa,
        fq1 = outdir + "/fastqs/{sample}_R1.fastq",
        fq2 = outdir + "/fastqs/{sample}_R2.fastq"
    output:
        out = directory(outdir + "/outputs/{sample}")
    log:
        outdir + "/outputs/{sample}.log"
    threads:
        16
    shell:
        """
        set +u
        source activate SPRINT
        # source activate SPRINT_threads
        sprint main -ss 0 -c 6 -p {threads} -1 {input.fq1} -2 {input.fq2} {input.fa} {output.out} `which bwa` `which samtools` &> {log}
        """

rule getA2I:
    input:
        rules.sprint.output.out
    output:
        txt = outdir + "/A2I/{sample}.tsv"
    log:
        outdir + "/A2I/{sample}.log"
    shell:
        """
        set +u; source activate SPRINT
        python ../common/SPRINT/utilities/getA2I.py 1 {input} {output} &> {log}
        """

rule get_read_number:
    input:
        bamdir = rules.sprint.output.out
    output:
        txt = outdir + "/read/{sample}.txt"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bamdir}/tmp/genome/all.bam | grep 'primary mapped' | awk '{{print $1}}' > {output.txt}
        """