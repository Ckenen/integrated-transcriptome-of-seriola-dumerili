#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
names = ["ngs", "tgs"]
outdir = "results/function"

DIAMOND_DB = "../common/ncbi_nr/db/nr.dmnd"

rule all:
    input:
        expand(outdir + "/isoforms/{name}.fa", name=names),
        expand(outdir + "/ORFfinder/{name}.fa", name=names),
        expand(outdir + "/diamond/{name}.nr.tsv", name=names),
        expand(outdir + "/eggnog/{name}.emapper.hits", name=names),

# Fetching isoform sequences from GTF annotation

rule get_isoform_sequences:
    input:
        gtf = lambda wildcards: GTFS[wildcards.name],
        fa = GENOME_FASTA
    output:
        fa = outdir + "/isoforms/{name}.fa"
    shell:
        """
        ./scripts/get_isoform_sequence_from_gtf.py {input.gtf} {input.fa} > {output.fa}
        """

# Predicting ORF regions and translated amino acids

rule ORFfinder:
    input:
        fa = rules.get_isoform_sequences.output.fa
    output:
        fa = outdir + "/ORFfinder/{name}.fa"
    log:
        outdir + "/ORFfinder/{name}.log"
    shell:
        """
        ORFfinder -in {input.fa} -strand plus -out {output.fa} &> {log}
        """

# Finding protein similarity

rule diamond:
    input:
        fa = rules.ORFfinder.output.fa,
        db = DIAMOND_DB
    output:
        tsv = outdir + "/diamond/{name}.nr.tsv"
    log:
        log = outdir + "/diamond/{name}.nr.log"
    threads:
        80
    shell:
        """
        diamond blastp -p {threads} --evalue 0.00001 --db {input.db} --query {input.fa} --strand plus --out {output.tsv} &> {log}
        """

rule EggNOG:
    input:
        fa = rules.ORFfinder.output.fa
    output:
        ann = outdir + "/eggnog/{name}.emapper.annotations",
        hit = outdir + "/eggnog/{name}.emapper.hits",
        ort = outdir + "/eggnog/{name}.emapper.seed_orthologs"
    log:
        outdir + "/eggnog/{name}.emapper.log"
    params:
        prefix = outdir + "/eggnog/{name}"
    threads:
        80
    shell:
        """
        set +u; source activate EggNOG
        ../common/software/EggNOG/eggnog-mapper-2.1.7/emapper.py -m diamond -i {input.fa} --output {params.prefix} -d euk --usemem --cpu {threads} &> {log}
        """
