#!/usr/bin/env python
configfile: "config.yaml"
samples = config["samples"]
indir = "results/mapping/rmdup"
outdir = "results/tracks"

rule all:
    input:
        # expand(outdir + "/bed/{sample}.bed", sample=samples),
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

rule bam_to_bed:
    input:  
        fam = indir + "/{sample}.fam"
    output:
        bed = outdir + "/bed/{sample}.bed",
        txt = outdir + "/bed/{sample}.txt"
    shell:
        """
        ../common/scripts/mapping/fam_to_bed.py {input.fam} {output.bed} > {output.txt}
        """

# BEDGRAPH

rule bed_to_bedgraph:
    input:
        bed = outdir + "/bed/{sample}.bed",
        genome = config["genome_sizes"]
    output:
        tmp = temp(directory(outdir + "/bg/{sample}.bg.SORT_TMP")),
        bg = outdir + "/bg/{sample}.bg"
    shell:
        """
        mkdir {output.tmp}
        bedtools genomecov -bg -split -i {input.bed} -g {input.genome} | sort -T {output.tmp} -k1,1 -k2,2n > {output.bg}
        """

rule bed_to_bedgraph_pos:
    input:
        bed = outdir + "/bed/{sample}.bed",
        genome = config["genome_sizes"]
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
        genome = config["genome_sizes"]
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
        bed = outdir + "/bed/{sample}.bed",
        txt = outdir + "/bed/{sample}.txt",
        genome = config["genome_sizes"]
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
        bed = outdir + "/bed/{sample}.bed",
        txt = outdir + "/bed/{sample}.txt",
        genome = config["genome_sizes"]
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
        bed = outdir + "/bed/{sample}.bed",
        txt = outdir + "/bed/{sample}.txt",
        genome = config["genome_sizes"]
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
        genome = config["genome_sizes"]
    output:
        bw = outdir + "/bw/{prefix}.bw"
    shell:
        """
        bedGraphToBigWig {input.bg} {input.genome} {output.bw}
        """

