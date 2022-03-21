#!/usr/bin/env python
configfile: "config.yaml"
samples = config["samples"]
outdir = "results/gatk"

rule all:
    input:
        expand(outdir + "/bam/{sample}.bam", sample=samples),
        expand(outdir + "/haplotype/{sample}.gvcf.gz", sample=samples),

rule add_rg_tag:
    input:
        bam = "results/mapping/markdup/{sample}.bam"
    output:
        bam = outdir + "/bam/{sample}.bam",
        bai = outdir + "/bam/{sample}.bam.bai"
    params:
        header = "@RG\\tID:{sample}\\tLB:{sample}\\tPU:{sample}\\tSM:{sample}"
    shell:
        """
        samtools addreplacerg -O BAM -r '{params.header}' -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule call_haplotype: # HaplotypeCaller
    input:
        bam = rules.add_rg_tag.output.bam,
        bai = rules.add_rg_tag.output.bam + ".bai",
        fasta = config["genome"]
    output:
        gvcf = outdir + "/haplotype/{sample}.gvcf.gz"
    log:
        log = outdir + "/haplotype/{sample}.log"
    threads:
        4
    shell:
        """
        gatk --java-options -Xmx4G HaplotypeCaller -I {input.bam} -O {output.gvcf} -R {input.fasta} --emit-ref-confidence GVCF &> {log}
        """
