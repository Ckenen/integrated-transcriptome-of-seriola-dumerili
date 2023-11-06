#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
samples = final_samples
indir = "results/denovo_mapping/rmdup"
outdir = "results/assembly"

rule all:
    input:
        expand(outdir + "/stringtie/gtf/{sample}.gtf", sample=samples),
        outdir + "/stringtie/merged_all_samples.gtf",
        outdir + "/taco/stringtie.gtf"

# StringTie

rule stringtie:
    input:
        bam = indir + "/{sample}.bam"
    output:
        gtf = outdir + "/stringtie/gtf/{sample}.gtf"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf --conservative -p {threads} -o {output.gtf} -l {wildcards.sample}
        """

rule stringtie_merge:
    input:
        gtfs = [outdir + "/stringtie/gtf/%s.gtf" % s for s in samples]
    output:
        gtf = outdir + "/stringtie/merged_all_samples.gtf",
        gtf_gz = outdir + "/stringtie/merged_all_samples.sorted.gtf.gz"
    shell:
        """
        stringtie --merge {input.gtfs} > {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """

# TACO

rule stringtie_taco:
    input:
        gtfs = expand(rules.stringtie.output.gtf, sample=samples)
    output:
        txt = outdir + "/taco/stringtie.filelist",
        out = directory(outdir + "/taco/stringtie"),
        gtf = outdir + "/taco/stringtie.gtf",
        gtf_gz = outdir + "/taco/stringtie.gtf.gz"
    log:
        log = outdir + "/taco/stringtie.log"
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
