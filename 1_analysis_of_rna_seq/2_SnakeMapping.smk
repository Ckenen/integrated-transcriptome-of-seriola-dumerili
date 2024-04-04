#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_ALL
INDIR = "data/datasets"
OUTDIR = "results/mapping"

rule all:
    input:
        OUTDIR + "/star/index",
        expand(OUTDIR + "/star/mapped/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/infered/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/rmdup/{sample}.bam", sample=SAMPLES),

rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf = ANNOTATION_GTF
    output:
        out = directory(OUTDIR + "/star/index")
    log:
        log = OUTDIR + "/star/index.log"
    conda:
        "star"
    threads:
        THREADS
    shell:
        """
        mkdir {output.out}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.out} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAindexNbases 13 \
            --genomeFastaFiles {input.fasta} &> {log}
        """

rule star_mapping:
    input:
        fq1 = INDIR + "/{sample}_R1.fastq.gz",
        fq2 = INDIR + "/{sample}_R2.fastq.gz", 
        idx = rules.star_index.output.out
    output:
        out = directory(OUTDIR + "/star/mapped/{sample}")
    log:
        log = OUTDIR + "/star/mapped/{sample}.log"
    conda:
        "star"
    threads:
        THREADS
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 10000000000 \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule filter_bam:
    input:
        bamdir = rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/filtered/{sample}.bam"
    log:
        OUTDIR + "/filtered/{sample}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} \
            -F 2308 \
            -f 3 \
            -d 'NH:1' \
            --expr 'rname =~ "^NW_"' \
            -o {output.bam} \
            {input.bamdir}/Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule infer_experiment:
    input:
        bam = rules.filter_bam.output.bam,
        bed = ANNOTATION_BED
    output:
        txt = OUTDIR + "/infered/{sample}.txt"
    log:
        OUTDIR + "/infered/{sample}.log"
    conda:
        "rseqc"
    shell:
        """
        infer_experiment.py \
            -i {input.bam} \
            -r {input.bed} > {output.txt} 2> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/rmdup/{sample}.bam",
        txt = OUTDIR + "/rmdup/{sample}_metrics.txt"
    log:
        OUTDIR + "/rmdup/{sample}.log"
    conda:
        "picard"
    shell:
        """
        picard MarkDuplicates -REMOVE_DUPLICATES true \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.txt} &> {log}
        samtools index -@ {threads} {output.bam}
        """
