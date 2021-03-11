#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "results/mapping/filtered"
outdir = "results/stringtie"


rule all:
    input:
        expand(outdir + "/round1/{sample}.gtf", sample=samples),
        outdir + "/merged/merged.gtf",
        outdir + "/merged/merged.sorted.gtf.gz",
        expand(outdir + "/round2/{sample}/gene_abund.tab", sample=samples),
        outdir + "/taco/merged",
        # outdir + "/taco/assembly.sorted.gtf.gz"


# StringTie 
# 1. de novo isoform;
# 2. merge gtf;
# 3. calculate expression (FPKM)

rule strintie_round1:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
        gtf = "data/genome/annotation.gtf"
    output:
        gtf = outdir + "/round1/{sample}.gtf",
    log:
        outdir + "/round1/{sample}.log"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf -G {input.gtf} -p {threads} -o {output.gtf} -l {wildcards.sample} &> {log}
        """

rule strintie_merge:
    input:
        gtfs = expand(rules.strintie_round1.output.gtf, sample=samples),
        gtf = "data/genome/annotation.gtf"
    output:
        txt = outdir + "/merged/merged.gtf.list",
        gtf = outdir + "/merged/merged.gtf"
    log:
        outdir + "/merged/merged.log"
    shell:
        """(
        for gtf in {input.gtfs}; do echo $gtf; done > {output.txt}
        stringtie --merge -G {input.gtf} -o {output.gtf} -l merge {output.txt} ) &> {log}
        """

rule stringtie_round2:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
        gtf = rules.strintie_merge.output.gtf
    output:
        gtf = outdir + "/round2/{sample}/transcripts.gtf", 
        tab = outdir + "/round2/{sample}/gene_abund.tab",
        ref = outdir + "/round2/{sample}/cov_refs.gtf",
        ct1 = outdir + "/round2/{sample}/i_data.ctab", # For Ballgown
        ct2 = outdir + "/round2/{sample}/e_data.ctab",
        ct3 = outdir + "/round2/{sample}/t_data.ctab",
        ct4 = outdir + "/round2/{sample}/i2t.ctab",
        ct5 = outdir + "/round2/{sample}/e2t.ctab",
    log:
        outdir + "/round2/{sample}.log"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf -e -A {output.tab} -C {output.ref} -B -p {threads} -G {input.gtf} -o {output.gtf} &> {log}
        """

# TACO

rule taco_merge:
    input:
        gtfs = expand(rules.strintie_round1.output.gtf, sample=samples)
    output:
        txt = outdir + "/taco/merged.gtf.txt",
        out = directory(outdir + "/taco/merged")
    log:
        outdir + "/taco/merged.log"
    params:
        gtf = outdir + "/taco/merged/assembly.gtf",
        bed = outdir + "/taco/merged/assembly.bed"
    threads:
        8
    shell:
        """
        for f in {input.gtfs}; do echo $f; done > {output.txt}
        taco_run -p {threads} -o {output.out} --gtf-expr-attr FPKM --filter-min-length 200 \
            --filter-min-expr 0.5 --isoform-frac 0.05 {output.txt} &> {log}
        """

# Common rules

rule gtf_index:
    input:
        gtf = "{prefix}.gtf"
    output:
        gtf = "{prefix}.sorted.gtf.gz",
        tbi = "{prefix}.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -header -i {input.gtf} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

