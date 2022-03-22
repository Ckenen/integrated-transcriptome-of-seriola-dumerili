#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/rmdup"
outdir = "results/expression"
refs = ["ncbi", "ensembl", "taco"]

# 流程： 计算FPKM、合并、计算相关性

rule all:
    input:
        expand(outdir + "/featureCounts/{ref}/{sample}.txt", ref=refs, sample=samples),
        expand(outdir + "/featureCounts/{ref}.feature_count.tsv", ref=refs),
        expand(outdir + "/featureCounts/{ref}.feature_count.corr.pdf", ref=refs),
        expand(outdir + "/stringtie/{ref}/{sample}", ref=refs, sample=samples),
        expand(outdir + "/stringtie/{ref}.gene_abund.tsv", ref=refs),
        expand(outdir + "/stringtie/{ref}.gene_abund.corr.pdf", ref=refs),
        expand(outdir + "/stringtie/{ref}.transcript_fpkm.tsv", ref=refs),
        expand(outdir + "/stringtie/{ref}.transcript_fpkm.corr.pdf", ref=refs),


def get_gtf(wildcards):
    ref = wildcards.ref
    if ref == "ncbi":
        return config["ncbi"]
    if ref == "ensembl":
        return config["ensembl"]
    if ref == "taco":
        return config["assembly"]
    assert False


rule uncompress_gtf:
    input:  
        gtf = lambda wildcards: get_gtf(wildcards)
    output:
        gtf = outdir + "/{ref}.gtf"
    shell:
        """
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output.gtf}
        """

# featureCounts

rule feature_counts:
    input:
        bam = indir + "/{sample}.bam",
        gtf = rules.uncompress_gtf.output.gtf
    output:
        txt = outdir + "/featureCounts/{ref}/{sample}.txt"
    log:
        log = outdir + "/featureCounts/{ref}/{sample}.log"
    threads:
        8
    shell:
        """
        featureCounts -T {threads} -s 2 -p -B -a {input.gtf} -o {output.txt} {input.bam} &> {log}
        """

rule merge_feature_counts:
    input:
        [outdir + "/featureCounts/{ref}/%s.txt" % s for s in samples]
    output:
        tsv = outdir + "/featureCounts/{ref}.feature_count.tsv"
    shell:
        """
        ../common/scripts/expression/merge_feature_counts.py {input} {output.tsv}
        """

# StringTie

rule stringtie_fpkm:
    input:
        bam = indir + "/{sample}.bam",
        gtf = rules.uncompress_gtf.output.gtf
    output:
        out = directory(outdir + "/stringtie/{ref}/{sample}")
    threads:
        8
    shell:
        """
        mkdir {output.out}
        stringtie {input.bam} --rf -e -A {output.out}/gene_abund.tab -C {output.out}/cov_refs.gtf \
            -B -p {threads} -G {input.gtf} -o {output.out}/transcripts.gtf
        """

rule merge_stringtie_gene_abund:
    input:
        [outdir + "/stringtie/{ref}/%s" % s for s in samples]
    output:
        tsv = outdir + "/stringtie/{ref}.gene_abund.tsv"
    shell:
        """
        ../common/scripts/expression/merge_stringtie_gene_abund.py {input} {output.tsv}
        """

rule merge_stringtie_transcript_fpkm:
    input:
        [outdir + "/stringtie/{ref}/%s" % s for s in samples]
    output:
        tsv = outdir + "/stringtie/{ref}.transcript_fpkm.tsv"
    shell:
        """
        ../common/scripts/expression/merge_stringtie_transcript_fpkm.py {input} {output.tsv}
        """

# common rules

rule plot_correlation:
    input:
        tsv = "{prefix}.tsv"
    output:
        pdf = "{prefix}.corr.pdf"
    shell:
        """
        ../common/scripts/expression/plot_correlation.py {input.tsv} {output.pdf}
        """
