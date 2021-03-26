#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = config["samples"]
threads: 20


rule all:
    input:
        expand("results/tama/collapsed/{sample}.bed", sample=samples),
        "results/tama/merged/tama_merged.bed",
        "results/tama/merged/tama_merged.sorted.gtf.gz",
        "results/tama/read_support/all_read_support.txt",
        "results/tama/remove_polya/filelist.txt",
        "results/tama/remove_polya/remove_polya.sorted.gtf.gz",
        

# Collapse
# 这一步需要很多运行内存，不要同时运行多个任务

rule tama_collapse:
    input:
        bam = "results/isoseq/aligned/{sample}.bam",
        bai = "results/isoseq/aligned/{sample}.bam.bai",
        fsa = "data/genome/genome.fasta"
    output:
        out1 = "results/tama/collapsed/{sample}.bed",
        out2 = "results/tama/collapsed/{sample}_local_density_error.txt",
        out3 = "results/tama/collapsed/{sample}_polya.txt",
        out4 = "results/tama/collapsed/{sample}_read.txt",
        out5 = "results/tama/collapsed/{sample}_strand_check.txt",
        out6 = "results/tama/collapsed/{sample}_trans_read.bed",
        out7 = "results/tama/collapsed/{sample}_trans_report.txt",
        out8 = "results/tama/collapsed/{sample}_varcov.txt",
        out9 = "results/tama/collapsed/{sample}_variants.txt",
    log:
        log = "results/tama/collapsed/{sample}.log"
    params:
        prefix = "results/tama/collapsed/{sample}"
    shell:
        """
        set +u; source activate py27
        python ~/software/tama/tama_collapse.py \
            -b BAM -s {input.bam} -f {input.fsa} -p {params.prefix} -x no_cap &> {log}
        conda deactivate
        """

# Merge

rule tama_merge:
    input:
        bed1 = "results/tama/collapsed/Adult_Female.bed",
        bed2 = "results/tama/collapsed/Adult_Male.bed",
        bed3 = "results/tama/collapsed/Juvenile.bed",
    output:
        out1 = "results/tama/merged/tama_merged.bed",
        out2 = "results/tama/merged/tama_merged_gene_report.txt",
        out3 = "results/tama/merged/tama_merged_merge.txt",
        out4 = "results/tama/merged/tama_merged_trans_report.txt",
        out5 = "results/tama/merged/tama_merged.filelist.txt"
    log:
        log = "results/tama/merged/tama_merged.log"
    params:
        prefix = "results/tama/merged/tama_merged"
    shell:
        """
        cat > {output.out5} << EOF
{input.bed1}	no_cap	1,1,1	AdFe
{input.bed2}	no_cap	1,1,1	AdMa
{input.bed3}	no_cap	1,1,1	JuMi
EOF
        set +u; source activate py27
        python ~/software/tama/tama_merge.py \
            -f {output.out5} -p {params.prefix} &> {log}
        conda deactivate
        """

rule compress_merged_gtf:
    input:
        bed = "results/tama/merged/tama_merged.bed"
    output:
        gtf1 = "results/tama/merged/tama_merged.gtf",
        gtf2 = "results/tama/merged/tama_merged.sorted.gtf.gz",
        tbi = "results/tama/merged/tama_merged.sorted.gtf.gz.tbi"
    shell:
        """
        set +u; source activate py27
        python ~/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
            {input.bed} {output.gtf1}
        conda deactivate
        bedtools sort -i {output.gtf1} | awk '$3!="gene"' | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2}
        """


# Support read

rule tama_read_support_levels:
    input:
        in1 = "results/tama/collapsed/Adult_Female_trans_read.bed",
        in2 = "results/tama/collapsed/Adult_Male_trans_read.bed",
        in3 = "results/tama/collapsed/Juvenile_trans_read.bed",
        in4 = "results/tama/merged/tama_merged_merge.txt",
    output:
        out1 = "results/tama/read_support/all_read_support.filelist.txt",
        out2 = "results/tama/read_support/all_read_support.txt"
    log:
        log = "results/tama/read_support/read_support.log"
    params:
        prefix = "results/tama/read_support/all"
    shell:
        """
        cat > {output.out1} << EOF
AdFe	{input.in1}	trans_read
AdMa	{input.in2}	trans_read
JuMi	{input.in3}	trans_read
EOF
        set +u; source activate py27
        python ~/software/tama/tama_go/read_support/tama_read_support_levels.py \
			-f {output.out1} -m {input.in4} -o {params.prefix} &> {log}
        conda deactivate
        """

# Filter polyA

rule tama_remove_polya_models_levels:
    input:
        in1 = "results/tama/collapsed/Adult_Female_polya.txt",
        in2 = "results/tama/collapsed/Adult_Male_polya.txt",
        in3 = "results/tama/collapsed/Juvenile_polya.txt",
        in4 = "results/tama/merged/tama_merged.bed",
        in5 = "results/tama/read_support/all_read_support.txt"
    output:
        out1 = "results/tama/remove_polya/filelist.txt",
        out2 = "results/tama/remove_polya/remove_polya_trash_polya.bed",
        out3 = "results/tama/remove_polya/remove_polya_polya_report.txt",
        out4 = "results/tama/remove_polya/remove_polya_polya_support.txt",
        out5 = "results/tama/remove_polya/remove_polya.bed",
    log:
        log = "results/tama/remove_polya/remove_polya.log"
    params:
        prefix = "results/tama/remove_polya/remove_polya"
    shell:
        """
        cat > {output.out1} << EOF
AdFe	{input.in1}
AdMa	{input.in2}
JuMi	{input.in3}
EOF
        set +u; source activate py27
        python ~/software/tama/tama_go/filter_transcript_models/tama_remove_polya_models_levels.py \
            -b {input.in4} -f {output.out1} -r {input.in5} -o {params.prefix} > {log}
        conda deactivate
        """

rule compress_filtered_gtf:
    input:
        bed = "results/tama/remove_polya/remove_polya.bed"
    output:
        gtf1 = "results/tama/remove_polya/remove_polya.gtf",
        gtf2 = "results/tama/remove_polya/remove_polya.sorted.gtf.gz",
        tbi = "results/tama/remove_polya/remove_polya.sorted.gtf.gz.tbi"
    shell:
        """
        set +u; source activate py27
        python ~/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
            {input.bed} {output.gtf1}
        conda deactivate
        bedtools sort -i {output.gtf1} | awk '$3!="gene"' | bgzip -c > {output.gtf2}
        tabix -p gff {output.gtf2}
        """
