#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 20
outdir = "results/mapping"

rule all:
    input:
        outdir + "/genome.splice.mmi",
        expand(outdir + "/flnc/fasta/{sample}.flnc.fasta", sample=samples),
        expand(outdir + "/flnc/mapped/{sample}.bam", sample=samples),
        # expand(outdir + "/mapped/{sample}.bam", sample=samples),
        # expand(outdir + "/mapped/{sample}.bam.bai", sample=samples),
        expand(outdir + "/markInternalPriming/{sample}.bam", sample=samples),
        expand(outdir + "/markInternalPriming/{sample}.bam.bai", sample=samples),
        expand(outdir + "/markInternalPriming/{sample}.stats", sample=samples),
        expand(outdir + "/clean/{sample}.bam", sample=samples),
        expand(outdir + "/clean/{sample}.bam.bai", sample=samples),
        expand(outdir + "/clean/{sample}.stats", sample=samples),

rule minimap2_index:
    input:
        fasta = config["genome_fasta"]
    output:
        mmi = outdir + "/genome.splice.mmi"
    threads:
        threads
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fasta}
        """

rule minimap2_align:
    input:
        fasta = "results/isoseq/polished/{sample}.hq.fasta.gz",
        mmi = rules.minimap2_index.output.mmi
    output:
        fasta = temp(outdir + "/mapped/{sample}.hq.fasta"),
        sam = temp(outdir + "/mapped/{sample}.tmp.sam"),
        tmp = temp(outdir + "/mapped/{sample}.tmp.bam"),
        bam = outdir + "/mapped/{sample}.bam",
    threads:
        threads
    shell:
        """
        pigz -d -c -p {threads} {input.fasta} > {output.fasta}
        minimap2 -ax splice -t {threads} --secondary=no -C5 {input.mmi} {output.fasta} > {output.sam} 2>/dev/null
        samtools view -@ {threads} -b {output.sam} > {output.tmp}
        samtools sort -@ {threads} {output.tmp} > {output.bam}
        """

rule flnc_to_fasta:
    input:
        bam = "results/isoseq/flnc/{sample}.bam"
    output:
        fasta = outdir + "/flnc/fasta/{sample}.flnc.fasta"
    shell:
        """
        samtools view {input.bam} | awk '{{print ">"$1"\\n"$10}}' > {output.fasta}
        """

rule flnc_align:
    input:
        fasta = rules.flnc_to_fasta.output.fasta,
        mmi = rules.minimap2_index.output.mmi
    output:
        sam = temp(outdir + "/flnc/mapped/{sample}.tmp.sam"),
        tmp = temp(outdir + "/flnc/mapped/{sample}.tmp.bam"),
        bam = outdir + "/flnc/mapped/{sample}.bam"
    log:
        outdir + "/flnc/mapped/{sample}.log"
    threads:
        threads
    shell:
        """(
        minimap2 -ax splice -t {threads} --secondary=no -C5 {input.mmi} {input.fasta} > {output.sam}
        samtools view -@ {threads} -b {output.sam} > {output.tmp}
        samtools sort -@ {threads} {output.tmp} > {output.bam} ) &> {log}
        """

# mark internal priming

rule mark_internal_priming:
    input:
        bam = rules.minimap2_align.output.bam,
        fasta = config["genome_fasta"]
    output:
        bam = outdir + "/markInternalPriming/{sample}.bam"
    shell:
        """
        ./scripts/mark_internal_priming.py {input.bam} {input.fasta} {output.bam}
        """

rule remove_internal_priming_alignments:
    input:
        bam = rules.mark_internal_priming.output.bam,
    output:
        bam = outdir + "/clean/{sample}.bam"
    shell:
        """
        bamtools filter -in {input.bam} -out {output.bam} -tag XP:0
        """

# common rules


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        bamtools index -in {input}
        """

rule bam_stat:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """