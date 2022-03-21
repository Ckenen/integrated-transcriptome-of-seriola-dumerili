#!/usr/bin/env snakemake
configfile: "config.yaml"
import glob
names = [os.path.basename(p).strip(".fa") for p in glob.glob("results/function/interpro/faSplited/*.fa")]
outdir = "results/function"

rule all:
    input:
        outdir + "/isoforms.fa",
        outdir + "/ORFfinder/ORFs.fa",
        # outdir + "/GeneMarkST/ORFs.gff",
        outdir + "/diamond/nr.tsv",
        outdir + "/eggnog/out.emapper.hits",
        outdir + "/interpro/faSplited",
        # expand(outdir + "/interpro/pfam/{name}.tsv", name=names),
        # expand(outdir + "/interpro/PANTHER/{name}.tsv", name=names),
        outdir + "/interpro/pfam.tsv",
        outdir + "/interpro/PANTHER.tsv",

# Fetching isoform sequences from GTF annotation

rule get_isoform_sequences:
    input:
        gtf = config["assembly"],
        fa = "../common/ncbi/serDum.ncbi.fasta"
    output:
        fa = outdir + "/isoforms.fa"
    shell:
        """
        ../common/scripts/assembly/get_sequence_from_gtf.py \
            {input.gtf} {input.fa} > {output.fa}
        """

# Predicting ORF regions and translated amino acids

rule ORFfinder:
    input:
        fa = rules.get_isoform_sequences.output.fa
    output:
        fa = outdir + "/ORFfinder/ORFs.fa"
    log:
        outdir + "/ORFfinder/ORFs.log"
    shell:
        """
        ORFfinder -in {input.fa} -strand plus -out {output.fa} &> {log}
        """

# Deprecated

rule GeneMarkST:
    input:
        fa = rules.get_isoform_sequences.output.fa
    output:
        gff = outdir + "/GeneMarkST/ORFs.gff",
        fa = outdir + "/GeneMarkST/ORFs.gff.faa"
    log:
        outdir + "/GeneMarkST/ORFs.log"
    shell:  
        """
        ../common/software/GeneMarkST/gmst.pl \
            --strand direct --filter 0 --gcode 1 --output {output.gff} \
            --format GFF --fnn --faa {input.fa} &> {log}
        """

# Finding protein similarity

rule diamond:
    input:
        fa = rules.ORFfinder.output.fa,
        dmnd = "../common/ncbi_nr/db/nr.dmnd"
    output:
        tsv = outdir + "/diamond/nr.tsv"
    log:
        log = outdir + "/diamond/nr.log"
    threads:
        80
    shell:
        """
        diamond blastp -p {threads} --evalue 0.00001 --db {input.dmnd} \
            --query {input.fa} --strand plus --out {output.tsv} &> {log}
        """

rule EggNOG:
    input:
        fa = rules.ORFfinder.output.fa
    output:
        ann = outdir + "/eggnog/out.emapper.annotations",
        hit = outdir + "/eggnog/out.emapper.hits",
        ort = outdir + "/eggnog/out.emapper.seed_orthologs"
    log:
        outdir + "/eggnog/out.emapper.log"
    params:
        prefix = outdir + "/eggnog/out"
    threads:
        80
    shell:
        """
        set +u; source activate EggNOG
        ../common/software/EggNOG/eggnog-mapper-2.1.7/emapper.py \
            -m diamond -i {input.fa} --output {params.prefix} \
            -d euk --usemem --cpu {threads} &> {log}
        """

# InterProScan

rule split_fasta:
    input:
        fa = rules.ORFfinder.output.fa
    output:
        out = directory(outdir + "/interpro/faSplited")
    shell:
        """
        mkdir {output}
        faSplit about {input.fa} 1000000 {output}/
        """

rule InterProScan_Pfam:
    input:
        fas = rules.split_fasta.output.out
    output:
        tsv = outdir + "/interpro/pfam/{name}.tsv",
        gff = outdir + "/interpro/pfam/{name}.gff3"
    log:
        log = outdir + "/interpro/pfam/{name}.log"
    threads:
        4
    params:
        prefix = outdir + "/interpro/pfam/{name}"
    shell:
        """
        set +u; source activate InterProScan
        ~/software/interproscan-5.52-86.0/interproscan.sh \
            -i {input.fas}/{wildcards.name}.fa -f tsv,gff3 \
            -appl Pfam -goterms -pathways \
            --cpu {threads} -b {params.prefix} -dp &> {log}
        conda deactivate
        """

rule InterProScan_Pfam_merge:
    input:
        tsvs = expand(outdir + "/interpro/pfam/{name}.tsv", name=names)
    output:
        tsv = outdir + "/interpro/pfam.tsv"
    shell:
        """
        cat {input.tsvs} > {output.tsv}
        """

rule InterProScan_PANTHER:
    input:
        fas = rules.split_fasta.output.out
    output:
        tsv = outdir + "/interpro/PANTHER/{name}.tsv",
        gff = outdir + "/interpro/PANTHER/{name}.gff3"
    log:
        log = outdir + "/interpro/PANTHER/{name}.log"
    threads:
        4
    params:
        prefix = outdir + "/interpro/PANTHER/{name}"
    shell:
        """
        set +u; source activate InterProScan
        ~/software/interproscan-5.52-86.0/interproscan.sh \
            -i {input.fas}/{wildcards.name}.fa -f tsv,gff3 \
            -appl PANTHER -goterms -pathways \
            --cpu {threads} -b {params.prefix} -dp &> {log}
        conda deactivate
        """

rule InterProScan_PANTHER_merge:
    input:
        tsvs = expand(outdir + "/interpro/PANTHER/{name}.tsv", name=names)
    output:
        tsv = outdir + "/interpro/PANTHER.tsv"
    shell:
        """
        cat {input.tsvs} > {output.tsv}
        """
