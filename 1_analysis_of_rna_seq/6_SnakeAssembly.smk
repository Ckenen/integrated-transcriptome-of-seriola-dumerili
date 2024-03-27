#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_FINAL
INDIR = "results/denovo_mapping/rmdup"
OUTDIR = "results/assembly"

rule all:
    input:
        expand(OUTDIR + "/stringtie/gtf/{sample}.gtf", sample=SAMPLES[:1]),
        OUTDIR + "/stringtie/merged_all_samples.gtf",
        OUTDIR + "/taco/stringtie.gtf"

# StringTie

rule stringtie:
    input:
        bam = INDIR + "/{sample}.bam"
    output:
        gtf = OUTDIR + "/stringtie/gtf/{sample}.gtf"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf --conservative -p {threads} -o {output.gtf} -l {wildcards.sample}
        """

rule stringtie_merge:
    input:
        gtfs = [OUTDIR + "/stringtie/gtf/%s.gtf" % s for s in SAMPLES]
    output:
        gtf = OUTDIR + "/stringtie/merged_all_samples.gtf",
        gtf_gz = OUTDIR + "/stringtie/merged_all_samples.sorted.gtf.gz"
    shell:
        """
        stringtie --merge {input.gtfs} > {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """

# TACO

rule stringtie_taco:
    input:
        gtfs = expand(rules.stringtie.output.gtf, sample=SAMPLES)
    output:
        txt = OUTDIR + "/taco/stringtie.filelist",
        out = directory(OUTDIR + "/taco/stringtie"),
        gtf = OUTDIR + "/taco/stringtie.gtf",
        gtf_gz = OUTDIR + "/taco/stringtie.gtf.gz"
    log:
        log = OUTDIR + "/taco/stringtie.log"
    threads:
        8
    shell:
        """
        for f in {input.gtfs}; do echo $f; done > {output.txt}
        taco_run -p {threads} -o {output.out} --filter-min-expr 1.0 {output.txt} &> {log}
        cp {output.out}/assembly.gtf {output.gtf}
        bedtools sort -header -i {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """
