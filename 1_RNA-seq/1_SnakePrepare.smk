#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import glob
samples = [os.path.basename(p)[:-12] for p in glob.glob("data/datasets/*_R1.fastq.gz")]
indir = "data/datasets"
outdir = "results/prepare"
rs = ["R1", "R2"]


rule all:
    input:
        # FastQC
        expand(outdir + "/fastqc/{sample}_{r}_fastqc.html", sample=samples, r=rs),
        # STAR
        outdir + "/star/index",
        expand(outdir + "/star/mapped/{sample}", sample=samples),
        # StringTie
        expand(outdir + "/stringtie/fpkm/{sample}", sample=samples),
        # Analysis
        outdir + "/stringtie/gene_abund.tsv",
        outdir + "/stringtie/gene_abund.corr.heatmap.pdf",
        outdir + "/stringtie/gene_abund.corr.clustermap.pdf",
        outdir + "/stringtie/gene_abund.clustermap.pdf"

# FastQC
# Perform quality control on the raw data.

rule fastqc:
    input:
        fq = indir + "/{name}.fastq.gz"
    output:
        html = outdir + "/fastqc/{name}_fastqc.html"
    log:
        log = outdir + "/fastqc/{name}.log"
    shell:
        """
        fastqc -o `dirname {output.html}` {input} &> {log}
        """

# GTF
# Uncompress gtf file in advance for STAR and StringTie.

rule uncompress_gtf:
    input:
        gtf = config["annotation"]
    output:
        gtf = outdir + "/annotation.gtf"
    shell:
        """
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output.gtf}
        """

# STAR
# Build index and align reads to genome.
# Run the following command when finished:
# STAR --genomeLoad Remove --genomeDir  results/prepare/star/index

rule star_index:
    input:
        fasta = config["genome"],
        gtf = rules.uncompress_gtf.output.gtf
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

rule star_mapping:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz", 
        idx = outdir + "/star/index",
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
            --alignEndsType Local \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 10000000000 \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --outFilterMultimapNmax 1 \
            --genomeLoad LoadAndKeep \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

# StringTie
# Calculating FPKM.

rule stringtie_fpkm:
    input:
        bam = rules.star_mapping.output.out,
        gtf = rules.uncompress_gtf.output.gtf
    output:
        out = directory(outdir + "/stringtie/fpkm/{sample}")
    threads:
        8
    shell:
        """
        mkdir {output.out}
        stringtie {input.bam}/Aligned.sortedByCoord.out.bam --rf \
            -e -A {output.out}/gene_abund.tab -C {output.out}/cov_refs.gtf \
            -B -p {threads} -G {input.gtf} -o {output.out}/transcripts.gtf
        """

# Merge and plot

rule merge_stringtie_gene_abund:
    input:
        expand(outdir + "/stringtie/fpkm/{sample}", sample=samples),
    output:
        tsv = outdir + "/stringtie/gene_abund.tsv"
    shell:
        """
        ../common/scripts/expression/merge_stringtie_gene_abund.py {input} {output.tsv}
        """

rule plot_fpkm_corr_heatmap:
    input:
        tsv = rules.merge_stringtie_gene_abund.output.tsv
    output:
        pdf = outdir + "/stringtie/gene_abund.corr.heatmap.pdf"
    shell:
        """
        ../common/scripts/expression/plot_fpkm_corr_heatmap.py {input.tsv} {output.pdf}
        """

rule plot_fpkm_corr_clustermap:
    input:
        tsv = rules.merge_stringtie_gene_abund.output.tsv
    output:
        pdf = outdir + "/stringtie/gene_abund.corr.clustermap.pdf"
    shell:
        """
        ../common/scripts/expression/plot_fpkm_corr_clustermap.py {input.tsv} {output.pdf}
        """

rule plot_fpkm_clustermap:
    input:
        tsv = rules.merge_stringtie_gene_abund.output.tsv
    output:
        pdf = outdir + "/stringtie/gene_abund.clustermap.pdf"
    shell:
        """
        ../common/scripts/expression/plot_fpkm_clustermap.py {input.tsv} {output.pdf}
        """
