#!/usr/bin/env snakemake

names = [
    "stringtie_stringtie", 
    "stringtie_taco",
    "cufflink_cuffmerge",
    "cufflink_taco"
]


rule all:
    input:
        expand("results/tama/{name}.bed", name=names),


rule gtf_to_bed:
    input:
        
    output:
        "results/tama/{name}.bed"
    shell:
        """
        """