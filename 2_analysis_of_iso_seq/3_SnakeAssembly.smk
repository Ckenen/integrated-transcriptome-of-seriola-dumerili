#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/filtered"
OUTDIR = "results/assembly"

rule all:
    input:
        expand(OUTDIR + "/collapsed/{sample}", sample=SAMPLES),
        OUTDIR + "/tama_merged_all_samples",
        OUTDIR + "/tama_merged_all_samples.bed",
        OUTDIR + "/tama_merged_all_samples.gtf",
        OUTDIR + "/tama_merged_all_samples.filtered.bed",
        OUTDIR + "/tama_merged_all_samples.filtered.gtf",
        
rule tama_collapse:
    input:
        bam = INDIR + "/{sample}.bam",
        fa = GENOME_FASTA
    output:
        sam = temp(OUTDIR + "/collapsed/{sample}.sam"),
        out = directory(OUTDIR + "/collapsed/{sample}")
    log:
        OUTDIR + "/collapsed/{sample}.log"
    conda:
        "TAMA"
    params:
        prefix = OUTDIR + "/collapsed/{sample}/{sample}"
    shell:
        """
        samtools view -h {input.bam} > {output.sam}
        mkdir {output.out}
        python {TAMA_ROOT}/tama_collapse.py \
            -s {output.sam} \
            -f {input.fa} \
            -x no_cap \
            -p {params.prefix} &> {log}
        """

rule tama_merge:
    input:
        beddir1 = OUTDIR + "/collapsed/Ad_Ma",
        beddir2 = OUTDIR + "/collapsed/Ad_Fe",
        beddir3 = OUTDIR + "/collapsed/Ju_Mi"
    output:
        txt = OUTDIR + "/tama_merged_all_samples.filelist",
        out = directory(OUTDIR + "/tama_merged_all_samples")
    log:
        OUTDIR + "/tama_merged_all_samples.log"
    conda:
        "TAMA"
    params:
        bed1 = OUTDIR + "/collapsed/Ad_Ma/Ad_Ma.bed",
        bed2 = OUTDIR + "/collapsed/Ad_Fe/Ad_Fe.bed",
        bed3 = OUTDIR + "/collapsed/Ju_Mi/Ju_Mi.bed"
    shell:
        """
        echo -e "{params.bed1}\\tno_cap\\t1,1,1\\tAdMa" >> {output.txt}
        echo -e "{params.bed2}\\tno_cap\\t1,1,1\\tAdFe" >> {output.txt}
        echo -e "{params.bed3}\\tno_cap\\t1,1,1\\tJuMi" >> {output.txt}
        mkdir {output.out}
        python {TAMA_ROOT}/tama_merge.py \
            -f {output.txt} -p {output.out}/all_samples &> {log}
        """

rule filter_zero_length_exon_isoform:
    input:
        beddir = rules.tama_merge.output.out
    output:
        bed = OUTDIR + "/tama_merged_all_samples.bed",
        bed2 = OUTDIR + "/tama_merged_all_samples.bed.gz"
    shell:
        """
        cat {input.beddir}/all_samples.bed \
            | ./scripts/filter_zero_length_exon_isoform.py > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2}
        """

rule bed_to_gtf:
    input:
        bed = rules.filter_zero_length_exon_isoform.output.bed
    output:
        gtf = OUTDIR + "/tama_merged_all_samples.gtf",
        gtf2 = OUTDIR + "/tama_merged_all_samples.gtf.gz"
    conda:
        "TAMA"
    shell:
        """
        python {TAMA_ROOT}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
            {input.bed} {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2}
        """

rule filter_internal_primer_isoform:
    input:
        bed = rules.filter_zero_length_exon_isoform.output.bed,
        fasta = GENOME_FASTA
    output:
        bed = OUTDIR + "/tama_merged_all_samples.filtered.bed",
        bed2 = OUTDIR + "/tama_merged_all_samples.filtered.bed.gz"
    shell:
        """
        ./scripts/filter_internal_primer_isoform.py {input} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2}
        """

rule tama_convert_bed_gtf_ensembl_no_cds:
    input:
        bed = rules.filter_internal_primer_isoform.output.bed
    output:
        gtf = OUTDIR + "/tama_merged_all_samples.filtered.gtf",
        gtf2 = OUTDIR + "/tama_merged_all_samples.filtered.gtf.gz"
    conda:
        "TAMA"
    shell:
        """
        python {TAMA_ROOT}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
            {input.bed} {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2}
        """