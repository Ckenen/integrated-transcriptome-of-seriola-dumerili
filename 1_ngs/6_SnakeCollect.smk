#!/usr/bin/env snakemake

names = [
    "stringtie_stringtie", 
    "stringtie_taco",
    "cufflink_cuffmerge",
    "cufflink_taco"
]

rule all:
    input:
        expand("results/collect/{name}.gtf", name=names),
        expand("results/collect/{name}.sorted.gtf.gz", name=names)


def get_input_gtf(name):
    if name == "stringtie_stringtie":
        return "results/stringtie/merged/merged.gtf"
    elif name == "stringtie_taco":
        return "results/stringtie/taco/merged/assembly.gtf"
    elif name == "cufflink_cuffmerge":
        return "results/cufflinks/merged.filtered.gtf"
    elif name == "cufflink_taco":
        return "results/cufflinks/taco/merged/assembly.gtf"


rule collect:
    input:
        gtf = lambda wildcards: get_input_gtf(wildcards.name)
    output:
        gtf = "results/collect/{name}.gtf"
    shell:
        """
        cp {input.gtf} {output.gtf}
        """

rule compress:
    input:
        gtf = "results/collect/{name}.gtf"
    output:
        gtf = "results/collect/{name}.sorted.gtf.gz",
        tbi = "results/collect/{name}.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -i {input.gtf} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """