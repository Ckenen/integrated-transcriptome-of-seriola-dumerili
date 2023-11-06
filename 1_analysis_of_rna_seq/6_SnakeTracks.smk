#!/usr/bin/env python
include: "0_SnakeCommon.smk"
samples = final_samples
indir = "results/denovo_mapping/rmdup"
outdir = "results/tracks"

rule all:
    input:
        # expand(outdir + "/bed/{sample}.bed", sample=samples),
        expand(outdir + "/library_size/{sample}.txt", sample=samples),
        # expand(outdir + "/bg/{sample}.bg", sample=samples),
        # expand(outdir + "/bg/{sample}.+.bg", sample=samples),
        # expand(outdir + "/bg/{sample}.-.bg", sample=samples),
        # expand(outdir + "/bg/{sample}.norm.bg", sample=samples),
        # expand(outdir + "/bg/{sample}.norm.+.bg", sample=samples),
        # expand(outdir + "/bg/{sample}.norm.-.bg", sample=samples),
        expand(outdir + "/bw/{sample}.bw", sample=samples),
        expand(outdir + "/bw/{sample}.+.bw", sample=samples),
        expand(outdir + "/bw/{sample}.-.bw", sample=samples),
        expand(outdir + "/bw/{sample}.norm.bw", sample=samples),
        expand(outdir + "/bw/{sample}.norm.+.bw", sample=samples),
        expand(outdir + "/bw/{sample}.norm.-.bw", sample=samples)

# BED

rule make_fragment_bed:
    input:  
        bam = indir + "/{sample}.bam"
    output:
        bed = outdir + "/bed/{sample}.bed"
    shell:
        """
        ./scripts/make_fragment_bed.py {input.bam} {output.bed}
        """

# Library size

rule get_library_size:
    input:
        bed = rules.make_fragment_bed.output.bed
    output:
        txt = outdir + "/library_size/{sample}.txt"
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
        tmp = temp(directory(outdir + "/bg/{sample}.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.bg",
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -i {input.bed} -g {input.gsize} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}

        """

rule bed_to_bedgraph_pos:
    input:
        bed = outdir + "/bed/{sample}.bed",
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(outdir + "/bg/{sample}.+.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.+.bg"
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -strand + -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_neg:
    input:
        bed = outdir + "/bed/{sample}.bed",
        genome = GENOME_SIZES
    output:
        tmp = temp(directory(outdir + "/bg/{sample}.-.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.-.bg"
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
        tmp = temp(directory(outdir + "/bg/{sample}.norm.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.norm.bg"
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
        tmp = temp(directory(outdir + "/bg/{sample}.norm.+.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.norm.+.bg"
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
        tmp = temp(directory(outdir + "/bg/{sample}.norm.-.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.norm.-.bg"
    shell:
        """
        mkdir {output.tmp}
        scale=`cat {input.txt} | awk '{{print 1000000/$0}}'`
        bedtools genomecov -bg -split -scale $scale -strand - -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """


# BIGWIG

rule bedgraph_to_bigwig:
    input:
        bg = outdir + "/bg/{prefix}.bg",
        txt = GENOME_SIZES
    output:
        bw = outdir + "/bw/{prefix}.bw"
    shell:
        """
        bedGraphToBigWig {input.bg} {input.txt} {output.bw}
        """

