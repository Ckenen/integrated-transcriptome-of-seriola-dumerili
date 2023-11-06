#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
samples = final_samples
indir = "data/datasets"
outdir = "results/denovo_mapping"

rule all:
    input:
        outdir + "/star/index",
        expand(outdir + "/star/mapped.1st/{sample}", sample=samples),
        expand(outdir + "/star/mapped.2nd/{sample}", sample=samples),
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/rmdup/{sample}.bam", sample=samples),

# Mapping to genome

rule star_index:
    input:
        fasta = GENOME_FASTA
    output:
        out = directory(outdir + "/star/index")
    log:
        log = outdir + "/star/index.log"
    threads:
        20
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} &> {log}
        """

rule star_mapping_1st:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz",
        idx = rules.star_index.output.out
    output:
        out = directory(outdir + "/star/mapped.1st/{sample}")
    log:
        log = outdir + "/star/mapped.1st/{sample}.log"
    threads:
        20
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
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz", 
        idx = rules.star_index.output.out,
        tabs = expand(outdir + "/star/mapped.1st/{sample}", sample=final_samples)
    output:
        out = directory(outdir + "/star/mapped.2nd/{sample}")
    log:
        log = outdir + "/star/mapped.2nd/{sample}.log"
    threads:
        20
    params:
        tabs = " ".join([outdir + "/star/mapped.1st/%s/SJ.out.tab" % s for s in samples])
    shell:
        """
        mkdir {output}
        STAR --genomeDir {input.idx} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS \
            --limitBAMsortRAM 10000000000 \
            --limitSjdbInsertNsj=2000000 \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --sjdbFileChrStartEnd {params.tabs} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

# Filtered

rule filter_bam:
    input:
        bamdir = rules.star_mapping_2nd.output.out
    output:
        bam = outdir + "/filtered/{sample}.bam"
    threads:
        8
    shell:
        """
        samtools view -@ {threads} -F 2308 -f 3 -d 'NH:1' --expr 'rname =~ "^NW_"' -o {output.bam} {input.bamdir}/Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """

# rmdup

rule remove_duplicate:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/rmdup/{sample}.bam",
        txt = outdir + "/rmdup/{sample}_metrics.txt"
    log:
        log = outdir + "/rmdup/{sample}.log"
    threads:
        8
    shell:
        """
        picard MarkDuplicates -REMOVE_DUPLICATES true -TAGGING_POLICY All \
            -I {input.bam} -O {output.bam} -M {output.txt} &> {log}
        samtools index -@ {threads} {output.bam}
        """


# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """