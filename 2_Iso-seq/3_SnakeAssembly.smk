#!/usr/bin/env snakemake
outdir = "results/assembly"
tgs = ["Ad_Fe", "Ad_Ma", "Ju_Mi"]

rule all:
    input:
        # change read_name and merge
        expand(outdir + "/change_read_name/{name}.bam", name=tgs),
        outdir + "/merged_bam/merged.bam",
        expand(outdir + "/change_read_name_in_fasta/{name}.hq.fasta.gz", name=tgs),
        outdir + "/merged_fasta/merged.hq.fasta.gz",
        expand(outdir + "/change_read_name_in_csv/{name}.cluster_report.csv", name=tgs),
        outdir + "/merged_csv/merged.cluster_report.csv",
        # cupcake
        outdir + "/cupcake/collapsed/all.collapsed.gff",
        outdir + "/cupcake/collapsed/all.collapsed.abundance.txt",
        outdir + "/cupcake/collapsed/all.collapsed.filtered.gff",
        outdir + "/cupcake/collapsed/all.filtered.gtf.gz", # final output

# Prepare

rule change_read_name:
    input:
        bam = "results/mapping/clean/{name}.bam"
    output:
        bam = outdir + "/change_read_name/{name}.bam",
        bai = outdir + "/change_read_name/{name}.bam.bai"
    shell:
        """
        ./scripts/change_read_name.py {input.bam} {wildcards.name} {output.bam}
        samtools index {output.bam}
        """

rule merge_bam:
    input:
        bam1 = outdir + "/change_read_name/Ad_Fe.bam",
        bam2 = outdir + "/change_read_name/Ad_Ma.bam",
        bam3 = outdir + "/change_read_name/Ju_Mi.bam"
    output:
        bam = outdir + "/merged_bam/merged.bam",
        bai = outdir + "/merged_bam/merged.bam.bai"
    shell:
        """
        samtools merge {output.bam} {input.bam1} {input.bam2} {input.bam3}
        samtools index {output.bam}
        """

rule change_read_name_in_fasta:
    input:
        fasta = "results/isoseq/polished/{name}.hq.fasta.gz"
    output:
        fasta = outdir + "/change_read_name_in_fasta/{name}.hq.fasta.gz"
    shell:
        """
        gzip -d -c {input.fasta} | sed 's/^>/>{wildcards.name}./g' | gzip -c > {output.fasta}
        """

rule merge_fasta:
    input:
        fasta1 = outdir + "/change_read_name_in_fasta/Ad_Fe.hq.fasta.gz",
        fasta2 = outdir + "/change_read_name_in_fasta/Ad_Ma.hq.fasta.gz",
        fasta3 = outdir + "/change_read_name_in_fasta/Ju_Mi.hq.fasta.gz"
    output:
        fasta = outdir + "/merged_fasta/merged.hq.fasta.gz"
    shell:
        """
        zcat {input} | gzip -c > {output.fasta}
        """

rule change_read_name_in_csv:
    input:
        csv = "results/isoseq/polished/{name}.cluster_report.csv"
    output:
        csv = outdir + "/change_read_name_in_csv/{name}.cluster_report.csv"
    shell:
        """
        sed 's/^transcript/{wildcards.name}.transcript/g' {input.csv} > {output.csv}
        """

rule merge_csv:
    input:
        csv1 = outdir + "/change_read_name_in_csv/Ad_Fe.cluster_report.csv",
        csv2 = outdir + "/change_read_name_in_csv/Ad_Ma.cluster_report.csv",
        csv3 = outdir + "/change_read_name_in_csv/Ju_Mi.cluster_report.csv"
    output:
        csv = outdir + "/merged_csv/merged.cluster_report.csv"
    shell:
        """
        head -n 1 {input.csv1} > {output.csv}
        cat {input} | grep -v cluster_id >> {output.csv}
        """

# Cupcake

rule collapse:
    input:
        fasta = rules.merge_fasta.output.fasta, # 序列名称里面包含有FLNC的数量
        bam = rules.merge_bam.output.bam
    output:
        fasta = temp(outdir + "/cupcake/collapsed/all.hq.fasta"),
        gff = outdir + "/cupcake/collapsed/all.collapsed.gff",
        gtf = outdir + "/cupcake/collapsed/all.gtf.gz",
        tbi = outdir + "/cupcake/collapsed/all.gtf.gz.tbi"
    params:
        prefix = outdir + "/cupcake/collapsed/all",
    log:
        log = outdir + "/cupcake/collapsed/collapse.log"
    shell:
        """(
        gzip -d -c {input.fasta} > {output.fasta}
        collapse_isoforms_by_sam.py --input {output.fasta} --bam {input.bam} --dun-merge-5-shorter -o {params.prefix}
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf} ) &> {log}
        """

rule get_abundance_post_collapse: # 获取每个isoform上的FLNC的数量
    input:
        csv = rules.merge_csv.output.csv,
        gff = rules.collapse.output.gff
    output:
        txt = outdir + "/cupcake/collapsed/all.collapsed.abundance.txt",
    params:
        prefix = outdir + "/cupcake/collapsed/all.collapsed",
    log:
        outdir + "/cupcake/collapsed/all.collapsed.abundance.log"
    shell:
        """
        get_abundance_post_collapse.py {params.prefix} {input.csv} &> {log}
        """

rule filter_away_subset: # Filter away 5' degraded isoforms
    input:
        gff = rules.collapse.output.gff,
        txt = rules.get_abundance_post_collapse.output.txt
    output:
        gff = outdir + "/cupcake/collapsed/all.collapsed.filtered.gff",
        gtf = outdir + "/cupcake/collapsed/all.filtered.gtf.gz",
        tbi = outdir + "/cupcake/collapsed/all.filtered.gtf.gz.tbi",
    params:
        prefix = outdir + "/cupcake/collapsed/all.collapsed",
    shell:
        """
        filter_away_subset.py {params.prefix}
        bedtools sort -i {output.gff} | grep -v 'NC_016870.1' | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """
