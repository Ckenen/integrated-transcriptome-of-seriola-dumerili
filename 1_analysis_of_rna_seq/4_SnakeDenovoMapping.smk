#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_FINAL
INDIR = "data/datasets"
OUTDIR = "results/denovo_mapping"

rule all:
    input:
        OUTDIR + "/star/index",
        expand(OUTDIR + "/star/mapped.1st/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/star/mapped.2nd/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.flagstat.txt", sample=SAMPLES),
        expand(OUTDIR + "/rmdup/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/rmdup/{sample}.flagstat.txt", sample=SAMPLES),

rule star_index:
    input:
        fasta = GENOME_FASTA
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
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} &> {log}
        """

rule star_mapping_1st:
    input:
        fq1 = INDIR + "/{sample}_R1.fastq.gz",
        fq2 = INDIR + "/{sample}_R2.fastq.gz",
        idx = rules.star_index.output.out
    output:
        out = directory(OUTDIR + "/star/mapped.1st/{sample}")
    log:
        log = OUTDIR + "/star/mapped.1st/{sample}.log"
    conda:
        "star"
    threads:
        THREADS
    shell:
        """
        mkdir {output}
        STAR --outSAMtype None \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

# STAR --genomeLoad Remove --genomeDir results/denovo_mapping/star/index

rule star_mapping_2nd:
    input:
        fq1 = INDIR + "/{sample}_R1.fastq.gz",
        fq2 = INDIR + "/{sample}_R2.fastq.gz", 
        idx = rules.star_index.output.out,
        sjdb = expand(OUTDIR + "/star/mapped.1st/{sample}", sample=SAMPLES)
    output:
        out = directory(OUTDIR + "/star/mapped.2nd/{sample}")
    log:
        log = OUTDIR + "/star/mapped.2nd/{sample}.log"
    conda:
        "star"
    threads:
        THREADS
    params:
        sjdb = " ".join([OUTDIR + "/star/mapped.1st/%s/SJ.out.tab" % s for s in SAMPLES])
    shell:
        """
        mkdir {output}
        STAR --genomeDir {input.idx} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS \
            --limitBAMsortRAM 10000000000 \
            --limitSjdbInsertNsj 2000000 \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --sjdbFileChrStartEnd {params.sjdb} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule filter_bam:
    input:
        bamdir = rules.star_mapping_2nd.output.out
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

rule remove_duplicate:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/rmdup/{sample}.bam",
        txt = OUTDIR + "/rmdup/{sample}_metrics.txt"
    log:
        log = OUTDIR + "/rmdup/{sample}.log"
    conda:
        "picard"
    threads:
        THREADS
    shell:
        """
        picard MarkDuplicates \
            -REMOVE_DUPLICATES true \
            -TAGGING_POLICY All \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.txt} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat.txt"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
