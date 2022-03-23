#!/usr/bin/env python
# 合并二代和三代组装结果
configfile: "config.yaml"
outdir = "results/assembly"
refs = ["ncbi", "ensembl"]

rule all:
    input:
        outdir + "/ngs_tgs.merged.gtf.gz",
        expand(outdir + "/sqanti3/{ref}_vs_ngs_tgs.merged", ref=refs),
        expand(outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.tracking", ref=refs),
        outdir + "/final.gtf.gz"
        # outdir + "/fpkm.tsv",
        

# 合并NGS和TGS，过滤掉chrM的注释

rule merge_ngs_and_tgs:
    input:
        ngs = "../1_RNA-seq/results/assembly/taco/stringtie.gtf.gz",
        tgs = "../2_Iso-seq/results/assembly/cupcake.gtf.gz",
        fpkm_ngs = "../1_RNA-seq/results/expression/stringtie/taco.transcript_fpkm.tsv",
        fpkm_tgs = "../2_Iso-seq/results/expression/stringtie/cupcake.transcript_fpkm.tsv",
        sqanti = "results/compare/sqanti3/tgs_vs_ngs/query_classification.txt"        
    output:
        tmp = temp(outdir + "/ngs_tgs.merged.gtf"),
        gtf = outdir + "/ngs_tgs.merged.gtf.gz",
        tbi = outdir + "/ngs_tgs.merged.gtf.gz.tbi"
    log:
        outdir + "/ngs_tgs.merged.log"
    shell:
        """(
        ./scripts/intergrate_ngs_tgs_isoforms.py {input} {output.tmp}
        bedtools sort -i {output.tmp} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf} ) &> {log}
        """

# SQANTI3

rule sqanti3:
    input:
        ref = "../common/ncbi/serDum.{ref}.optimized.gtf.gz",
        gtf = rules.merge_ngs_and_tgs.output.gtf,
        fasta = config["GENOME"]
    output:       
        directory(outdir + "/sqanti3/{ref}_vs_ngs_tgs.merged")
    log:
        log = outdir + "/sqanti3/{ref}_vs_ngs_tgs.merged.log"
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
        ref = "../common/ncbi/serDum.{ref}.optimized.gtf.gz",
        gtf = rules.merge_ngs_and_tgs.output.gtf
    output:
        ref = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.ref.gtf",
        gtf = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.query.gtf",
        out1 = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.loci",
        out2 = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.tracking",
        out3 = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.annotated.gtf",
    log:
        log = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged.log"
    params:
        prefix = outdir + "/gffcompare/{ref}_vs_ngs_tgs.merged"
    shell:
        """
        gzip -d -c {input.ref} > {output.ref}
        gzip -d -c {input.gtf} > {output.gtf}
        gffcompare -r {output.ref} -R -o {params.prefix} {output.gtf} &> {log}
        """

# 将与NCBI比较之后的结果最为最终注释，并且添加chrM注释

rule final_gtf:
    input:
        gtf = outdir + "/sqanti3/ncbi_vs_ngs_tgs.merged"
    output:
        gtf = outdir + "/final.gtf.gz",
        tbi = outdir + "/final.gtf.gz.tbi"
    shell:
        """
        ( zcat ../common/ncbi/serDum.ncbi.optimized.gtf.gz | grep NC_016870.1
        bedtools sort -i {input.gtf}/query_corrected.gtf.cds.gff ) | bgzip -c > {output.gtf}
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

