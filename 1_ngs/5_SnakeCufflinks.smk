#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/cufflinks"


rule all:
    input:
        expand(outdir + "/clout/{sample}", sample=samples),
        expand(outdir + "/postlink/{sample}.gtf", sample=samples),
        "results/cufflinks/merged",
        outdir + "/merged.filtered.gtf",
        outdir + "/taco/merged",
        outdir + "/taco/assembly.sorted.gtf.gz",


# Cufflinks：主要的输出结果是transcripts.gtf，里面有一些feature是没有strand信息的（极少数）
# 输出的feature的范围包括soft-clip？这将导致start或者end范围越界

rule cufflinks:
    input:
        bam = indir + "/{sample}.bam",
        gtf = "data/genome/annotation.gtf"
    output:
        directory(outdir + "/clout/{sample}")
    log:
        outdir + "/clout/{sample}.log"
    threads:
        threads
    shell:
        """
        cufflinks -p {threads} -o {output} -g {input.gtf} --library-type fr-firststrand \
            -L {wildcards.sample} {input.bam} &> {log}
        """

# 由于cufflinks的输出有一些问题，例如越界，无strand信息等，所以添加这一步，过滤掉这些问题
rule postlinks:
    input:
        gsz = "data/genome/genome.sizes",
        gtf = outdir + "/clout/{sample}"
    output:
        gtf1 = outdir + "/postlink/{sample}.gtf",
        gtf2 = outdir + "/postlink/{sample}.sorted.gtf.gz",
        tbi = outdir + "/postlink/{sample}.sorted.gtf.gz.tbi"
    params:
        gtf = outdir + "/clout/{sample}/transcripts.gtf"
    log:
        log = outdir + "/postlink/{sample}.log"
    shell:
        """(
        ./scripts/post_links.py {input.gsz} {params.gtf} > {output.gtf1}
        bedtools sort -i {output.gtf1} | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2} ) &> {log}
        """

rule cuffmerge:
    input:
        gtf = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta",
        lis = [outdir + "/postlink/{sample}.gtf".format(sample=sample) for sample in samples]
    output:
        txt = outdir + "/merged.gtf.txt",
        out = directory(outdir + "/merged")
    log:
        outdir + "/merged.log"
    threads:
        threads
    shell:
        """
        for f in {input.lis}; do echo $f; done > {output.txt}
        cuffmerge -g {input.gtf} -s {input.fsa} -o {output.out} -p {threads} {output.txt}  &> {log}
        """

rule post_cuffmerge:
    input:
        outdir + "/merged"
    output:
        out1 = outdir + "/merged.filtered.gtf",
        out2 = outdir + "/merged.filtered.sorted.gtf.gz",
        out3 = outdir + "/merged.filtered.sorted.gtf.gz.tbi"
    shell:
        """
        awk '$7!="."' {input}/merged.gtf > {output.out1}
        bedtools sort -i {output.out1} | bgzip -c > {output.out2}
        tabix -p gff {output.out2}
        """
        
# TACO

rule taco_merge:
    input:
        gtfs = [outdir + "/postlink/{sample}.gtf".format(sample=sample) for sample in samples]
    output:
        txt = outdir + "/taco/merged.gtf.txt",
        out = directory(outdir + "/taco/merged")
    log:
        outdir + "/taco/merged.log"
    params:
        gtf = outdir + "/taco/merged/assembly.gtf",
        bed = outdir + "/taco/merged/assembly.bed"
    threads:
        threads
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