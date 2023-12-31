#!/usr/bin/env python
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/smrt"

rule all:
    input:
        # SMRT-link pipeline
        expand(outdir + "/min_passes_{p}/ccs/{sample}.bam", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/demux/{sample}.primer_5p--primer_3p.bam", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/flnc/{sample}.bam", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/clustered/{sample}.bam", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/polished/{sample}.bam", p=min_passes_list, sample=samples),
        # Mapping and collapse
        outdir + "/minimap2/genome.isoseq.mmi",
        expand(outdir + "/min_passes_{p}/minimap2/aligned/{sample}.bam", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/minimap2/aligned/{sample}.flagstat.txt", p=min_passes_list, sample=samples),
        expand(outdir + "/min_passes_{p}/collapsed/{sample}.gtf.gz", p=min_passes_list, sample=samples),

# Generate circular consensus sequences (ccs) from subreads.
# time-consuming

rule ccs:
    input:
        bam = indir + "/{sample}.subreads.bam",
    output:
        bam = outdir + "/min_passes_{p}/ccs/{sample}.bam"
    threads:
        threads
    shell:
        """
        ccs --skip-polish --min-passes {wildcards.p} --min-length 200 --min-rq 0.99 \
            --num-threads {threads} {input.bam} {output.bam}
        """

# Lima, Demultiplex Barcoded PacBio Data and Clip Barcodes

rule lima:
    input:
        bam = rules.ccs.output.bam,
        primer = "data/primers.fasta"
    output:
        outdir + "/min_passes_{p}/demux/{sample}.primer_5p--primer_3p.bam"
    params:
        bam = outdir + "/min_passes_{p}/demux/{sample}.bam"
    threads:
        threads
    shell:
        """
        lima --isoseq -j {threads} {input} {params.bam}
        """

# Remove polyA and concatemers from FL reads and generate FLNC transcripts (FL to FLNC)
# full-length, non-concatemer (FLNC) reads

rule refine:
    input:
        bam = rules.lima.output,
        primer = "data/primers.fasta"
    output:
        bam = outdir + "/min_passes_{p}/flnc/{sample}.bam"
    threads:
        threads
    shell:
        """
        isoseq3 refine --require-polya --min-polya-length 20 \
            -j {threads} {input} {output}
        """

# Cluster FLNC reads and generate unpolished transcripts (FLNC to UNPOLISHED)
# time-consuming

rule cluster:
    input:
        bam = rules.refine.output.bam
    output:
        bam = outdir + "/min_passes_{p}/clustered/{sample}.bam"
    threads:
        threads
    shell:
        """
        isoseq3 cluster -j {threads} {input.bam} {output.bam}
        """

# Polish transcripts using subreads (UNPOLISHED to POLISHED)
# time-consuming

rule polish:
    input:
        bam1 = rules.cluster.output.bam,
        bam2 = indir + "/{sample}.subreads.bam",
        pbi2 = indir + "/{sample}.subreads.bam.pbi",
    output:
        bam = outdir + "/min_passes_{p}/polished/{sample}.bam"
    log:
        outdir + "/min_passes_{p}/polished/{sample}.log"
    threads:
        threads
    shell:
        """
        isoseq3 polish -j {threads} {input.bam1} {input.bam2} {output} &> {log}
        """

# Index reference and store as .mmi file

rule mmindex_isoseq:
    input:
        fasta = GENOME_FASTA
    output:
        mmi = outdir + "/minimap2/genome.isoseq.mmi"
    threads:
        threads
    shell:
        """
        pbmm2 index -j {threads} --preset ISOSEQ {input.fasta} {output.mmi}
        """

# Align PacBio reads to reference sequences

rule align:
    input:
        mmi = rules.mmindex_isoseq.output.mmi,
        bam = rules.polish.output.bam
    output:
        bam = outdir + "/min_passes_{p}/minimap2/aligned/{sample}.bam"
    threads:
        threads
    shell:
        """
        pbmm2 align -j {threads} --preset ISOSEQ --sort {input.mmi} {input.bam} {output.bam}
        """

# Collapse transcripts based on genomic mapping

rule collapse:
    input:
        bam = rules.align.output.bam,
        ccs = rules.ccs.output.bam
    output:
        gff = outdir + "/min_passes_{p}/collapsed/{sample}.gff",
        gtf = outdir + "/min_passes_{p}/collapsed/{sample}.gtf.gz"
    threads:
        threads
    shell:
        """
        isoseq3 collapse -j {threads} {input.bam} {input.ccs} {output.gff}
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# Common 

# rule pacbio_index:
#     input:
#         "{prefix}.bam"
#     output:
#         "{prefix}.bam.pbi"
#     shell:
#         """
#         pbindex {input}
#         """

# rule bam_index:
#     input:
#         "{prefix}.bam"
#     output:
#         "{prefix}.bam.bai"
#     shell:
#         """
#         bamtools index -in {input}
#         """

rule bam_stats:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """
