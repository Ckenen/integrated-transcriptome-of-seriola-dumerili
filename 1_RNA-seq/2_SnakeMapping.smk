#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = config["threads"]
indir = "data/datasets"
outdir = "results/mapping"

# 流程：构建索引、第一次比对、第二次比对、去除PCR重复

rule all:
    input:
        outdir + "/star/index",
        expand(outdir + "/star/mapped.1st/{sample}", sample=samples),
        expand(outdir + "/star/mapped.2nd/{sample}", sample=samples),
        
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/filtered/{sample}.bam.bai", sample=samples),
        expand(outdir + "/filtered/{sample}.flagstat", sample=samples),
        expand(outdir + "/filtered/{sample}_mt.bam", sample=samples),
        expand(outdir + "/filtered/{sample}_mt.bam.bai", sample=samples),
        expand(outdir + "/filtered/{sample}_mt.flagstat", sample=samples),

        expand(outdir + "/infered/{sample}.txt", sample=samples),
        
        expand(outdir + "/markdup/{sample}.bam", sample=samples),
        expand(outdir + "/markdup/{sample}.bam.bai", sample=samples),
        expand(outdir + "/markdup/{sample}.flagstat", sample=samples),
        
        expand(outdir + "/rmdup/{sample}.bam", sample=samples),
        expand(outdir + "/rmdup/{sample}.bam.bai", sample=samples),
        expand(outdir + "/rmdup/{sample}.flagstat", sample=samples),

        expand(outdir + "/rmdup/{sample}.fam", sample=samples),

        # expand(outdir + "/rmdup/{sample}.chrom_read_count.tsv", sample=samples),
        # expand(outdir + "/uniq/{sample}.SJ.bed.gz", sample=samples),
        # expand(outdir + "/uniq/{sample}.SJ.summary.txt", sample=samples),
        # expand(outdir + "/uniq/{sample}.fam", sample=samples),
        # expand(outdir + "/uniq/{sample}.lengths.txt", sample=samples),
        # expand(outdir + "/uniq/{sample}.bw", sample=samples),

# Mapping to genome

rule star_index:
    input:
        fasta = config["genome"]
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
        10
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

# STAR --genomeLoad Remove --genomeDir results/mapping/star/index

rule star_mapping_2nd:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz", 
        idx = rules.star_index.output.out,
        tabs = expand(outdir + "/star/mapped.1st/{sample}", sample=samples)
    output:
        out = directory(outdir + "/star/mapped.2nd/{sample}")
    log:
        log = outdir + "/star/mapped.2nd/{sample}.log"
    threads:
        20
    params:
        tabs = " ".join([outdir + "/star/mapped.1st/%s/SJ.out.tab" % s for s in samples]),
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
            --outFilterMultimapNmax 1 \
            --sjdbFileChrStartEnd {params.tabs} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

# Filtered

rule filter_alignment:
    input:
        rules.star_mapping_2nd.output.out
    output:
        bam = outdir + "/filtered/{sample}.bam",
        bam_mt = outdir + "/filtered/{sample}_mt.bam"
    shell:
        """
        ../common/scripts/mapping/filter_alignments.py \
            {input}/Aligned.sortedByCoord.out.bam \
            {output.bam}
        """

# Infer experiment

rule uncompress_bed:
    input:
        bed = config["ncbi_bed"]
    output:
        bed = outdir + "/infered/ref.bed"
    shell:
        """
        gzip -d -c {input.bed} > {output.bed}
        """

rule infer_experiment:
    input:
        bam = rules.filter_alignment.output.bam,
        bai = rules.filter_alignment.output.bam + ".bai",
        bed = rules.uncompress_bed.output.bed
    output:
        txt = outdir + "/infered/{sample}.txt"
    shell:
        """
        set +u; source activate py27
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        conda deactivate 
        """

# MarkDup

rule mark_duplicates:
    input:
        bam = rules.filter_alignment.output.bam,
        bai = rules.filter_alignment.output.bam + ".bai"
    output:
        bam = outdir + "/markdup/{sample}.bam",
        txt = outdir + "/markdup/{sample}_metrics.txt"
    log:
        log = outdir + "/markdup/{sample}.log"
    threads:
        8
    shell:
        """(
        set +u; source activate picard
        picard MarkDuplicates \
            -REMOVE_DUPLICATES false \
            -TAGGING_POLICY All \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.txt} 
        conda deactivate ) &> {log}
        """

# Uniq

rule remove_duplicate:
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bam + ".bai"
    output:
        bam = outdir + "/rmdup/{sample}.bam"
    log:
        log = outdir + "/rmdup/{sample}.log"
    threads:
        8
    shell:
        """
        samtools view -F 1024 -b -@ {threads} -o {output.bam} {input.bam}
        """
        # """
        # ./scripts/mapping/remove_duplicates.py {input.bam} {output.bam}
        # """

# Common rules

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        samtools index {input}
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

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """

rule bam_to_fam:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai"
    output:
        fam = "{prefix}.fam"
    threads:
        8
    shell:
        """
        famtools.py build {input.bam} {output.fam} &> /dev/null
        """

rule extract_splice_junction:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai",
        fsa = config["genome"]
    output:
        tmp = temp("{prefix}.SJ.unsorted.bed"),
        bed = "{prefix}.SJ.bed.gz",
        tbi = "{prefix}.SJ.bed.gz.tbi"
    shell:
        """
        ./scripts/extract_splice_junction.py {input.bam} {input.fsa} -rf {output.tmp}
        bedtools sort -i {output.tmp} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule splice_junction_motif:
    input:
        bed = "{prefix}.SJ.bed.gz"
    output:
        txt = "{prefix}.SJ.summary.txt"
    shell:
        """
        zcat {input.bed} | awk '{{print $4}}' | sort | uniq -c | sort -k1,1nr | awk '{{print $2"\\t"$1}}' > {output.txt}
        """

rule get_chrom_read_count:
    input:
        bam = "{prefix}.bam"
    output:
        tsv = "{prefix}.chrom_read_count.tsv"
    shell:
        """
        ./scripts/stat_chrom_read_count.py {input.bam} > {output.tsv}
        """
        # """
        # samtools view {input.bam} | awk '{{print $3}}' | sort | uniq -c | awk '{{print $2"\\t"$1}}' > {output.txt}
        # """

rule infer_fragment_length:
    input:
        fam = "{prefix}.fam"
    output:
        txt = "{prefix}.lengths.txt"
    log:
        log = "{prefix}.lengths.log"
    shell:
        """
        ./scripts/infer_fragment_length.py {input.fam} {output.txt} &> {log}
        """

rule fam_to_bigwig:
    input:
        fam = "{prefix}.fam"
    output:
        bw1 = "{prefix}.bw",
        bw2 = "{prefix}.+.bw",
        bw3 = "{prefix}.-.bw",
        bw4 = "{prefix}.norm.bw",
        bw5 = "{prefix}.norm.+.bw",
        bw6 = "{prefix}.norm.-.bw",
    shell:
        """
        ./scripts/fam_to_bigwig.sh --rf {input.fam} {wildcards.prefix} &> /dev/null
        """
