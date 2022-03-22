#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "results/mapping/rmdup"
outdir = "results/assembly"

rule all:
    input:
        expand(outdir + "/stringtie/assemblied/{sample}.gtf", sample=samples),
        outdir + "/taco/stringtie.gtf.gz"

# StringTie

rule stringtie:
    input:
        bam = indir + "/{sample}.bam"
    output:
        gtf = outdir + "/stringtie/assemblied/{sample}.gtf"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf --conservative -p {threads} -o {output.gtf} -l {wildcards.sample}
        """

# TACO

rule stringtie_taco:
    input:
        gtfs = expand(rules.stringtie.output.gtf, sample=samples)
    output:
        txt = temp(outdir + "/taco/stringtie.filelist.txt"),
        out = directory(outdir + "/taco/stringtie"),
        gtf = outdir + "/taco/stringtie.gtf.gz",
        tbi = outdir + "/taco/stringtie.gtf.gz.tbi"
    log:
        log = outdir + "/taco/stringtie.log"
    threads:
        8
    shell:
        """
        for f in {input.gtfs}; do echo $f; done > {output.txt}
        taco_run -p {threads} -o {output.out} --filter-min-expr 1.0 {output.txt} &> {log}
        bedtools sort -header -i {output.out}/assembly.gtf | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """
