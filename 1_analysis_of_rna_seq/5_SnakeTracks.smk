#!/usr/bin/env python
include: "0_SnakeCommon.smk"
SAMPLES = SAMPLES_FINAL
INDIR = "results/denovo_mapping/rmdup"
OUTDIR = "results/tracks"

rule all:
    input:
        expand(OUTDIR + "/bed/{sample}.bed", sample=SAMPLES),
        expand(OUTDIR + "/library_size/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.bg", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.+.bg", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.-.bg", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.norm.bg", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.norm.+.bg", sample=SAMPLES),
        expand(OUTDIR + "/bg/{sample}.norm.-.bg", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.bw", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.+.bw", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.-.bw", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.norm.bw", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.norm.+.bw", sample=SAMPLES),
        expand(OUTDIR + "/bw/{sample}.norm.-.bw", sample=SAMPLES)

# BED

rule make_fragment_bed:
    input:  
        bam = INDIR + "/{sample}.bam"
    output:
        bed = OUTDIR + "/bed/{sample}.bed"
    shell:
        """
        ./scripts/make_fragment_bed.py {input.bam} {output.bed}
        """

# Library size

rule get_library_size:
    input:
        bed = rules.make_fragment_bed.output.bed
    output:
        txt = OUTDIR + "/library_size/{sample}.txt"
    shell:
        """
        wc -l {input.bed} > {output.txt}
        """

# BEDGRAPH

rule bed_to_bedgraph:
    input:
        bed = rules.make_fragment_bed.output.bed,
        gsize = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.bg",
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -i {input.bed} -g {input.gsize} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}

        """

rule bed_to_bedgraph_pos:
    input:
        bed = OUTDIR + "/bed/{sample}.bed",
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.+.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.+.bg"
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -strand + -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_neg:
    input:
        bed = OUTDIR + "/bed/{sample}.bed",
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.-.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.-.bg"
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -strand - -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_norm:
    input:
        bed = rules.make_fragment_bed.output.bed,
        txt = rules.get_library_size.output.txt,
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.norm.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.norm.bg"
    shell:
        """
        mkdir {output.tmp}
        scale=`cat {input.txt} | awk '{{print 1000000/$0}}'`
        bedtools genomecov -bg -split -scale $scale -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_norm_pos:
    input:
        bed = rules.make_fragment_bed.output.bed,
        txt = rules.get_library_size.output.txt,
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.norm.+.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.norm.+.bg"
    shell:
        """
        mkdir {output.tmp}
        scale=`cat {input.txt} | awk '{{print 1000000/$0}}'`
        bedtools genomecov -bg -split -scale $scale -strand + -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_norm_neg:
    input:
        bed = rules.make_fragment_bed.output.bed,
        txt = rules.get_library_size.output.txt,
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(OUTDIR + "/bg/{sample}.norm.-.bg.SORT_TMP")),
        bg = OUTDIR + "/bg/{sample}.norm.-.bg"
    shell:
        """
        mkdir {output.tmp}
        scale=`cat {input.txt} | awk '{{print 1000000/$0}}'`
        bedtools genomecov -bg -split -scale $scale -strand - -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """


# BIGWIG

rule bedgraph_to_bigwig:
    input:
        bg = OUTDIR + "/bg/{prefix}.bg",
        txt = GENOME_SIZES
    output:
        bw = OUTDIR + "/bw/{prefix}.bw"
    shell:
        """
        bedGraphToBigWig {input.bg} {input.txt} {output.bw}
        """

