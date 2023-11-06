#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
outdir = "results/mapping"

rule all:
    input:
        outdir + "/minimap2.splice.mmi",
        expand(outdir + "/minimap2/{sample}.mp{p}.bam", sample=samples, p=min_passes_list),
        expand(outdir + "/filtered/{sample}.mp{p}.bam", sample=samples, p=min_passes_list),
        expand(outdir + "/stat_clip/{sample}.mp{p}.bam", sample=samples, p=min_passes_list),
        expand(outdir + "/stat_internal_priming/{sample}.mp{p}.bam", sample=samples, p=min_passes_list),
        expand(outdir + "/filtered_internal_priming/{sample}.mp{p}.bam", sample=samples, p=min_passes_list),

# Minimap2

rule build_mm2_index:
    input:
        fa = GENOME_FASTA
    output:
        mmi = outdir + "/minimap2.splice.mmi"
    threads:
        threads
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fa}
        """

rule minimap2:
    input:
        mmi = rules.build_mm2_index.output.mmi,
        fq = "results/smrt/min_passes_{p}/polished/{sample}.hq.fastq.gz"
    output:
        bam = outdir + "/minimap2/{sample}.mp{p}.bam"
    log:
        outdir + "/minimap2/{sample}.mp{p}.log"
    params:
        rg = "@RG\\tID:{sample}\\tLB:{sample}\\tSM:{sample}"
    threads:
        threads
    shell:
        """(
        minimap2 -a -x splice:hq -u f -Y --MD -R '{params.rg}' -t {threads} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u \
            | samtools sort -@ {threads} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = outdir + "/filtered/{sample}.mp{p}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -F 2308 -q 30 --expr 'rname =~ "^NW_"' -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/stat_clip/{sample}.mp{p}.bam",
        txt = outdir + "/stat_clip/{sample}.mp{p}.tsv"
    log:
        outdir + "/stat_clip/{sample}.mp{p}.log"
    threads:
        4
    shell:
        """
        ./scripts/stat_clip.py --max-clip-3 5 -s {output.txt} -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# internal priming

rule stat_internal_priming:
    input:
        bam = rules.stat_clip.output.bam,
        fa = GENOME_FASTA
    output:
        bam = outdir + "/stat_internal_priming/{sample}.mp{p}.bam"
    threads:
        4
    shell:
        """
        ./scripts/stat_internal_priming.py {input.bam} {input.fa} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule filter_internal_priming:
    input:
        bam = rules.stat_internal_priming.output.bam,
    output:
        bam = outdir + "/filtered_internal_priming/{sample}.mp{p}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --tag 'XP:0' -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """
