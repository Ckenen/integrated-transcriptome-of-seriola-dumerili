#!/usr/bin/env python
configfile: "config.yaml"
names = config["builds"]
outdir = "results/compare"

pairs = []
for name1 in names:
    for name2 in names:
        if name1 != name2:
            pairs.append("%s_vs_%s" % (name1, name2))

rule all:
    input:
        expand(outdir + "/sqanti3/{pair}", pair=pairs),
        expand(outdir + "/gffcompare/{pair}.tracking", pair=pairs),

# SQANTI3

rule sqanti3:
    input:
        ref = lambda wildcards: config[wildcards.ref],
        gtf = lambda wildcards: config[wildcards.que],
        fasta = config["genome"]
    output:       
        directory(outdir + "/sqanti3/{ref}_vs_{que}")
    log:
        log = outdir + "/sqanti3/{ref}_vs_{que}.log"
    threads:
        8
    shell:
        """(
        mkdir {output}
        gzip -d -c {input.ref} | awk '$3!="gene"' > {output}/reference.gtf
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output}/query.gtf
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py -d {output} -t {threads} \
            --gtf {output}/query.gtf {output}/reference.gtf {input.fasta} 
        conda deactivate ) &> {log}
        """

# gffcompare

rule gffcompare:
    input:
        ref = lambda wildcards: config[wildcards.ref],
        gtf = lambda wildcards: config[wildcards.que]
    output:
        ref = outdir + "/gffcompare/{ref}_vs_{que}.ref.gtf",
        gtf = outdir + "/gffcompare/{ref}_vs_{que}.query.gtf",
        out1 = outdir + "/gffcompare/{ref}_vs_{que}.loci",
        out2 = outdir + "/gffcompare/{ref}_vs_{que}.tracking",
        out3 = outdir + "/gffcompare/{ref}_vs_{que}.annotated.gtf",
    log:
        log = outdir + "/gffcompare/{ref}_vs_{que}.log"
    params:
        prefix = outdir + "/gffcompare/{ref}_vs_{que}"
    shell:
        """
        gzip -d -c {input.ref} | awk '$3!="gene"' > {output.ref}
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output.gtf}
        gffcompare -r {output.ref} -R -o {params.prefix} {output.gtf} &> {log}
        """