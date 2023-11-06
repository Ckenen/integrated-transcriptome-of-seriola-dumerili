#!/usr/bin/env python
include: "0_SnakeCommon.smk"
outdir = "results/compare"
names = ["ncbi", "ensembl", "ngs", "tgs"]

pairs = []
for name1 in names:
    for name2 in names:
        if name1 != name2:
            pairs.append("%s_vs_%s" % (name1, name2))

rule all:
    input:
        expand(outdir + "/gtfs/{name}.gtf", name=names),
        expand(outdir + "/sqanti3/{pair}", pair=pairs),
        expand(outdir + "/gffcompare/{pair}.tracking", pair=pairs),

rule make_gtf:
    input:
        gtf = lambda wildcards: GTFS[wildcards.name]
    output:
        gtf = outdir + "/gtfs/{name}.gtf"
    shell:
        """
        zcat {input.gtf} > {output.gtf}
        """

# SQANTI3

rule sqanti3:
    input:
        ref = outdir + "/gtfs/{ref}.gtf",
        gtf = outdir + "/gtfs/{que}.gtf",
        fasta = GENOME_FASTA
    output:       
        directory(outdir + "/sqanti3/{ref}_vs_{que}")
    log:
        log = outdir + "/sqanti3/{ref}_vs_{que}.log"
    threads:
        8
    shell:
        """
        mkdir {output}
        awk '$3!="gene"' {input.ref} > {output}/reference.gtf
        awk '$3!="gene"' {input.gtf} > {output}/query.gtf
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py -d {output} -t {threads} --skipORF --gtf {output}/query.gtf {output}/reference.gtf {input.fasta} &> {log}
        conda deactivate
        """

# gffcompare

rule gffcompare:
    input:
        ref = outdir + "/gtfs/{ref}.gtf",
        gtf = outdir + "/gtfs/{que}.gtf"
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
        awk '$3!="gene"' {input.ref} > {output.ref}
        awk '$3!="gene"' {input.gtf} > {output.gtf}
        gffcompare -r {output.ref} -R -o {params.prefix} {output.gtf} &> {log}
        """