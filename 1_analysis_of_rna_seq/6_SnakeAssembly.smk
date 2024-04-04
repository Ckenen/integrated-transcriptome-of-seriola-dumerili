#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_FINAL
INDIR = "results/denovo_mapping/rmdup"
OUTDIR = "results/assembly"

rule all:
    input:
        expand(OUTDIR + "/gtfs/{sample}.gtf", sample=SAMPLES),
        OUTDIR + "/stringtie_merged_all_samples.gtf"

rule stringtie:
    input:
        bam = INDIR + "/{sample}.bam"
    output:
        gtf = OUTDIR + "/gtfs/{sample}.gtf"
    log:
        OUTDIR + "/gtfs/{sample}.log"
    conda:
        "stringtie"
    threads:
        4
    shell:
        """
        stringtie {input.bam} \
            --rf \
            --conservative \
            -p {threads} \
            -o {output.gtf} \
            -l {wildcards.sample} &> {log}
        """

rule merge_gtfs:
    input:
        gtfs = expand(rules.stringtie.output.gtf, sample=SAMPLES)
    output:
        gtf = OUTDIR + "/stringtie_merged_all_samples.gtf",
        gtf2 = OUTDIR + "/stringtie_merged_all_samples.gtf.gz"
    log:
        OUTDIR + "/stringtie_merged_all_samples.log"
    conda:
        "stringtie"
    shell:
        """
        stringtie --merge {input.gtfs} > {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2}
        """
