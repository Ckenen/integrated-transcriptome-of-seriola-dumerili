#!/usr/bin/env python
include: "0_SnakeCommon.smk"
outdir = "results/assembly"
refs = ["ncbi", "ensembl"]

rule all:
    input:
        outdir + "/asm.gtf",
        expand(outdir + "/sqanti3/{ref}_vs_asm", ref=refs),
        expand(outdir + "/gffcompare/{ref}_vs_asm.tracking", ref=refs),
        outdir + "/asm.final.sorted.gtf.gz",
        # outdir + "/fpkm.tsv",
        
rule merge_ngs_and_tgs:
    input:
        ngs = GTFS["ngs"],
        tgs = GTFS["tgs"],
        fpkm_ngs = "results/expression/stringtie/ngs.transcript_fpkm.tsv",
        fpkm_tgs = "results/expression/stringtie/tgs.transcript_fpkm.tsv",
        sqanti = "results/compare/sqanti3/tgs_vs_ngs/query_classification.txt"        
    output:
        gtf = outdir + "/asm.gtf",
        gtf_gz = outdir + "/asm.sorted.gtf.gz",
    log:
        outdir + "/asm.log"
    shell:
        """
        ./scripts/intergrate_ngs_tgs_isoforms.py {input} {output.gtf} &> {log}
        bedtools sort -i {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """

# SQANTI3

rule sqanti3:
    input:
        ref = lambda wildcards: GTFS[wildcards.ref],
        gtf = rules.merge_ngs_and_tgs.output.gtf_gz,
        fasta = GENOME_FASTA
    output:       
        out = directory(outdir + "/sqanti3/{ref}_vs_asm")
    log:
        outdir + "/sqanti3/{ref}_vs_asm.log"
    threads:
        8
    shell:
        """
        mkdir {output}
        gzip -d -c {input.ref} | awk '$3!="gene"' > {output}/reference.gtf
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output}/query.gtf
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py -d {output} -t {threads} --gtf {output}/query.gtf {output}/reference.gtf {input.fasta} &> {log}
        conda deactivate
        """

# gffcompare

rule gffcompare:
    input:
        ref = lambda wildcards: GTFS[wildcards.ref],
        gtf = rules.merge_ngs_and_tgs.output.gtf_gz
    output:
        ref = outdir + "/gffcompare/{ref}_vs_asm.ref.gtf",
        gtf = outdir + "/gffcompare/{ref}_vs_asm.query.gtf",
        out1 = outdir + "/gffcompare/{ref}_vs_asm.loci",
        out2 = outdir + "/gffcompare/{ref}_vs_asm.tracking",
        out3 = outdir + "/gffcompare/{ref}_vs_asm.annotated.gtf",
    log:
        outdir + "/gffcompare/{ref}_vs_asm.log"
    params:
        prefix = outdir + "/gffcompare/{ref}_vs_asm"
    shell:
        """
        gzip -d -c {input.ref} | awk '$3!="gene"' > {output.ref}
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output.gtf}
        gffcompare -r {output.ref} -R -o {params.prefix} {output.gtf} &> {log}
        """

# final gtf

rule final_gtf:
    input:
        gtf = outdir + "/sqanti3/ncbi_vs_asm"
    output:
        gtf = outdir + "/asm.final.sorted.gtf.gz"
    shell:
        """
        bedtools sort -i {input.gtf}/query_corrected.gtf.cds.gff | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# rule get_fpkm:
#     input:
#         ngs = "results/quality_control/expression.coincided/ngs.stringtie.taco.fpkm.tsv",
#         tgs = "results/quality_control/expression.coincided/tgs.cupcake.fpkm.tsv",
#         gtf = outdir + "/final.gtf.gz"
#     output:
#         tsv = outdir + "/fpkm.tsv"
#     shell:
#         """
#         ./scripts/get_fpkm.py {input} {output}
#         """

