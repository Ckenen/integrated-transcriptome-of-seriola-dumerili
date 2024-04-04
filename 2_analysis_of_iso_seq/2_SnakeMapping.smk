#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = "results/mapping"

rule all:
    input:
        OUTDIR + "/minimap2.splice.mmi",
        expand(OUTDIR + "/minimap2/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/stat_clip/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/stat_internal_priming/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/filtered_internal_priming/{sample}.bam", sample=SAMPLES),

# Minimap2

rule minimap2_index:
    input:
        fa = GENOME_FASTA
    output:
        mmi = OUTDIR + "/minimap2.splice.mmi"
    log:
        OUTDIR + "/minimap2.splice.log"
    conda:
        "minimap2"
    threads:
        THREADS
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        mmi = rules.minimap2_index.output.mmi,
        fq = "results/smrt/polished/{sample}.hq.fastq.gz"
    output:
        bam = OUTDIR + "/minimap2/{sample}.bam"
    log:
        OUTDIR + "/minimap2/{sample}.log"
    conda:
        "minimap2"
    params:
        rg = "@RG\\tID:{sample}\\tLB:{sample}\\tSM:{sample}"
    threads:
        THREADS
    shell:
        """(
        minimap2 -a -x splice:hq -u f -Y --MD -R '{params.rg}' \
            -t {threads} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u \
            | samtools sort -@ {threads} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/filtered/{sample}.bam"
    log:
        OUTDIR + "/filtered/{sample}.log"
    conda:
        "minimap2"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -F 2308 -q 30 \
            --expr 'rname =~ "^NW_"' \
            -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

# rule stat_clip:
#     input:
#         bam = rules.filter_bam.output.bam
#     output:
#         bam = OUTDIR + "/stat_clip/{sample}.bam",
#         txt = OUTDIR + "/stat_clip/{sample}.tsv"
#     log:
#         OUTDIR + "/stat_clip/{sample}.log"
#     shell:
#         """(
#         ./scripts/stat_clip.py \
#             --max-clip-3 5 \
#             -s {output.txt} \
#             -o {output.bam} {input.bam}
#         samtools index {output.bam} ) &> {log}
#         """

# rule stat_internal_priming:
#     input:
#         bam = rules.stat_clip.output.bam,
#         fa = GENOME_FASTA
#     output:
#         bam = OUTDIR + "/stat_internal_priming/{sample}.bam"
#     log:
#         OUTDIR + "/stat_internal_priming/{sample}.log"
#     shell:
#         """(
#         ./scripts/stat_internal_priming.py \
#             {input.bam} {input.fa} {output.bam}
#         samtools index {output.bam} ) &> {log}
#         """

# rule filter_internal_priming:
#     input:
#         bam = rules.stat_internal_priming.output.bam,
#     output:
#         bam = OUTDIR + "/filtered_internal_priming/{sample}.bam"
#     log:
#         OUTDIR + "/filtered_internal_priming/{sample}.log"
#     threads:
#         4
#     shell:
#         """(
#         samtools view -@ {threads} \
#             --tag 'XP:0' \
#             -o {output.bam} {input.bam}
#         samtools index -@ {threads} {output.bam} ) &> {log}
#         """
