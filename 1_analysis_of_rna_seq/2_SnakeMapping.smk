#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/mapping"
#samples = samples[3:]

rule all:
    input:
        outdir + "/star/index",
        expand(outdir + "/star/mapped/{sample}", sample=samples),
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/infered/{sample}.txt", sample=samples),
        expand(outdir + "/rmdup/{sample}.bam", sample=samples),

# Build index and align reads to genome.

rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf = ANNOTATION_GTF
    output:
        out = directory(outdir + "/star/index")
    log:
        log = outdir + "/star/index.log"
    threads:
        20
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

# Run the following command when finished:
# STAR --genomeLoad Remove --genomeDir  results/mapping/star/index

rule star_mapping:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz", 
        idx = rules.star_index.output.out,
    output:
        out = directory(outdir + "/star/mapped/{sample}")
    log:
        log = outdir + "/star/mapped/{sample}.log"
    threads:
        20
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
        bam = outdir + "/filtered/{sample}.bam"
    threads:
        8
    shell:
        """
        samtools view -@ {threads} -F 2308 -f 3 -d 'NH:1' -o {output.bam} {input.bamdir}/Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """

rule infer_experiment:
    input:
        bam = rules.filter_bam.output.bam,
        bed = ANNOTATION_BED
    output:
        txt = outdir + "/infered/{sample}.txt"
    shell:
        """
        set +u; source activate py27
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        conda deactivate 
        """

rule rmdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/rmdup/{sample}.bam",
        txt = outdir + "/rmdup/{sample}_metrics.txt"
    log:
        outdir + "/rmdup/{sample}.log"
    threads:
        20
    shell:
        """
        picard MarkDuplicates -REMOVE_DUPLICATES true \
            -I {input.bam} -O {output.bam} -M {output.txt} &> {log}
        samtools index -@ {threads} {output.bam}
        """
