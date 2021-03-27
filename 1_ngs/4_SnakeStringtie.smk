#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
indir = "results/mapping/filtered"
outdir = "results/stringtie"


rule all:
    input:
        # StringTie
        expand(outdir + "/round1/{sample}.gtf", sample=samples),
        outdir + "/merged/merged.gtf",
        outdir + "/merged/merged.sorted.gtf.gz",    # final output
        expand(outdir + "/round2/{sample}/gene_abund.tab", sample=samples),
        # TACO
        outdir + "/taco/merged",
        outdir + "/taco/assembly.sorted.gtf.gz",    # final output
        # outdir + "/taco/cpat",
        # outdir + "/taco/orf_finder/orf.fasta",
        # outdir + "/taco/blastp/blastp.tsv"


# StringTie 
# 1. de novo isoform;
# 2. merge gtf;
# 3. calculate expression (FPKM)

rule strintie_round1:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
        gtf = "data/genome/annotation.gtf"
    output:
        gtf = outdir + "/round1/{sample}.gtf",
    log:
        outdir + "/round1/{sample}.log"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf -G {input.gtf} -p {threads} \
            -o {output.gtf} -l {wildcards.sample} &> {log}
        """

# StringTie合并的结果里面有一些feature是没有strand信息的（极少部分）
# 这将导致下游分析报错，所以要过滤掉这些feature。
rule strintie_merge:
    input:
        gtfs = expand(rules.strintie_round1.output.gtf, sample=samples),
        gtf = "data/genome/annotation.gtf",
    output:
        txt = outdir + "/merged/merged.gtf.list",
        tmp = temp(outdir + "/merged/merged.tmp.gtf"),
        gtf = outdir + "/merged/merged.gtf"
    log:
        outdir + "/merged/merged.log"
    shell:
        """(
        for gtf in {input.gtfs}; do echo $gtf; done > {output.txt}
        stringtie --merge -G {input.gtf} -o {output.tmp} -l StringTie {output.txt} 
        cat {output.tmp} | awk '$7!="."' > {output.gtf} ) &> {log}
        """

rule stringtie_round2:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
        gtf = rules.strintie_merge.output.gtf
    output:
        gtf = outdir + "/round2/{sample}/transcripts.gtf", 
        tab = outdir + "/round2/{sample}/gene_abund.tab",
        ref = outdir + "/round2/{sample}/cov_refs.gtf",
        ct1 = outdir + "/round2/{sample}/i_data.ctab", # For Ballgown
        ct2 = outdir + "/round2/{sample}/e_data.ctab",
        ct3 = outdir + "/round2/{sample}/t_data.ctab",
        ct4 = outdir + "/round2/{sample}/i2t.ctab",
        ct5 = outdir + "/round2/{sample}/e2t.ctab",
    log:
        outdir + "/round2/{sample}.log"
    threads:
        8
    shell:
        """
        stringtie {input.bam} --rf -e -A {output.tab} -C {output.ref} -B -p {threads} -G {input.gtf} -o {output.gtf} &> {log}
        """

# TACO

rule taco_merge:
    input:
        gtfs = expand(rules.strintie_round1.output.gtf, sample=samples)
    output:
        txt = outdir + "/taco/merged.gtf.txt",
        out = directory(outdir + "/taco/merged")
    log:
        outdir + "/taco/merged.log"
    params:
        gtf = outdir + "/taco/merged/assembly.gtf",
        bed = outdir + "/taco/merged/assembly.bed"
    threads:
        8
    shell:
        """
        for f in {input.gtfs}; do echo $f; done > {output.txt}
        taco_run -p {threads} -o {output.out} --gtf-expr-attr FPKM --filter-min-length 200 \
            --filter-min-expr 0.5 --isoform-frac 0.05 {output.txt} &> {log}
        """

rule taco_gtf_index:
    input:
        outdir + "/taco/merged"
    output:
        gtf = outdir + "/taco/assembly.sorted.gtf.gz",
        tbi = outdir + "/taco/assembly.sorted.gtf.gz.tbi"
    params:
        gtf = outdir + "/taco/merged/assembly.gtf"
    shell:
        """
        bedtools sort -header -i {params.gtf} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule taco_cpat:
    input:
        fsa = "data/genome/genome.fasta",
        rda = "../ncbi/Sdu_1.0/GTS.logit.RData",
        tsv = "../ncbi/Sdu_1.0/GTS_Hexamer.tsv",
        tac = outdir + "/taco/merged"
    output:
        directory(outdir + "/taco/cpat")
    log:
        outdir + "/taco/cpat.log"
    params:
        bed = outdir + "/taco/merged/assembly.bed"
    shell:
        """
        mkdir {output}
        cpat.py -r {input.fsa} -g {params.bed} -d {input.rda} -x {input.tsv} -o {output}/out &> {log}
        """

rule taco_orf_finder:
    input:
        fsa = "data/genome/genome.fasta",
        tac = outdir + "/taco/merged"
    output:
        seq = outdir + "/taco/orf_finder/seq.fasta",
        orf = outdir + "/taco/orf_finder/orf.fasta"
    log:
        outdir + "/taco/orf_finder/orf.log"
    params:
        bed = outdir + "/taco/merged/assembly.bed"
    shell:
        """
        bedtools getfasta -s -split -name -fi {input.fsa} -bed {params.bed} | sed 's/([0-9.+-]*)//g' | sed 's/G[0-9]*|//g' > {output.seq}
        ORFfinder -in {output.seq} -strand plus -out {output.orf} &> {log}
        """

rule taco_blastp:
    input:
        orf = outdir + "/taco/orf_finder/orf.fasta",
    output:
        tsv = outdir + "/taco/blastp/blastp.tsv"
    log:
        outdir + "/taco/blastp/blastp.log"
    params:
        pdb = "../ncbi/swissprot/swissprot",
        # pdb = "../ncbi/protein/nr",
        max_target_seqs = 6,
        evalue = "1e-6"
    threads:
        80
    shell:
        """
        nice blastp -db {params.pdb} -outfmt 6 -evalue {params.evalue} -num_threads {threads} \
            -max_target_seqs {params.max_target_seqs} -query {input.orf} -out {output} &> {log}
        """

# Common rules

rule compress_gtf:
    input:
        gtf = "{prefix}.gtf"
    output:
        gtf = "{prefix}.sorted.gtf.gz",
        tbi = "{prefix}.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -header -i {input.gtf} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

