#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 20
outdir = "results/mapping"

rule all:
    input:
        outdir + "/genome.splice.mmi",
        expand(outdir + "/mapped/{sample}.bam", sample=samples),
        # expand(outdir + "/mapped/{sample}.filtered.bam", sample=samples),

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
        bai = outdir + "/mapped/{sample}.bam.bai",
    threads:
        threads
    shell:
        """
        pigz -d -c -p {threads} {input.fasta} > {output.fasta}
        minimap2 -ax splice -t {threads} --secondary=no -C5 {input.mmi} {output.fasta} > {output.sam} 2>/dev/null
        samtools view -@ {threads} -b {output.sam} > {output.tmp}
        samtools sort -@ {threads} {output.tmp} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# rule remove_internal_priming_reads:
#     input:
#         bam = rules.minimap2_align.output.bam,
#         fasta = config["genome_fasta"]
#     output:
#         bam = outdir + "/mapping/{sample}.filtered.bam",
#         bai = outdir + "/mapping/{sample}.filtered.bam.bai",
#         tsv = outdir + "/mapping/{sample}.filtered.tsv"
#     shell:
#         """
#         ./scripts/remove_internal_priming_reads.py {input.bam} {input.fasta} {output.bam} > {output.tsv}
#         samtools index {output.bam}
#         """