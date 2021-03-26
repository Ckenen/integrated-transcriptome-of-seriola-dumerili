#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = ["Adult_Male", "Adult_Female", "Juvenile"]
threads = 20

rule all:
    input:
        # Mapping
        expand("results/cupcake/mapping/{sample}.hq.fasta", sample=samples),
        "results/cupcake/genome.mmi",
        expand("results/cupcake/mapping/{sample}.sam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.sam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.bam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.bam.bai", sample=samples),
        # Collapse
        expand("results/cupcake/collapsed/{sample}/out.collapsed.gff", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz", sample=samples),
        # Chain
        "results/cupcake/chain/all_samples.chained.gff",
        "results/cupcake/chain/all_samples.chained.sorted.gtf.gz",
        # Fusion
        expand("results/cupcake/fusion/{sample}/out.abundance.txt", sample=samples),
        expand("results/cupcake/fusion.sqanti3/{sample}/out_classification.txt", sample=samples),
        expand("results/cupcake/fusion/{sample}/out.annotated.txt", sample=samples),


# HQ FASTA

rule get_hq_fasta:
    input:
        fsa = "results/isoseq/polished/{sample}.hq.fasta.gz"
    output:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta"
    shell:
        """
        gzip -d -c {input.fsa} > {output.fsa}
        """

# 使用Minimap2进行比对

rule minimap2_index:
    input:
        fsa = "data/genome/genome.fasta"
    output:
        mmi = "results/cupcake/genome.mmi"
    threads:
        threads
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fsa}
        """

rule minimap2_align:
    input:
        mmi = "results/cupcake/genome.mmi",
        fsa = "results/cupcake/mapping/{sample}.hq.fasta"
    output:
        sam = "results/cupcake/mapping/{sample}.sam"
    log:
        log = "results/cupcake/mapping/{sample}.log"
    threads:
        threads
    shell:
        """
        minimap2 -ax splice -t {threads} --secondary=no -C5 {input.mmi} {input.fsa} > {output.sam} 2> {log}
        """

rule sort_sam:
    input:
        sam = "results/cupcake/mapping/{sample}.sam"
    output:
        sam = "results/cupcake/mapping/{sample}.sorted.sam"
    shell:
        """
        cat {input.sam} | grep '^@' > {output.sam}
        cat {input.sam} | grep -v '^@' | sort -k3,3 -k4,4n >> {output.sam}
        """

rule sam_to_bam:
    input:
        sam = "results/cupcake/mapping/{sample}.sorted.sam"
    output:
        tmp = temp("results/cupcake/mapping/{sample}.tmp.bam"),
        bam = "results/cupcake/mapping/{sample}.sorted.bam"
    shell:
        """
        samtools view -b {input.sam} > {output.tmp}
        samtools sort {output.tmp} > {output.bam}
        """

rule bam_index:
    input:
        bam = "results/cupcake/mapping/{sample}.sorted.bam"
    output:
        bai = "results/cupcake/mapping/{sample}.sorted.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Cupcake

rule collapse:
    input:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta", # 序列名称里面包含有FLNC的数量
        sam = "results/cupcake/mapping/{sample}.sorted.sam"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.gff.unfuzzy",
        out3 = "results/cupcake/collapsed/{sample}/out.collapsed.group.txt",
        out4 = "results/cupcake/collapsed/{sample}/out.collapsed.group.txt.unfuzzy",
        out5 = "results/cupcake/collapsed/{sample}/out.ignored_ids.txt",
        out6 = "results/cupcake/collapsed/{sample}/out.collapsed.rep.fa",
    params:
        out = "results/cupcake/collapsed/{sample}/out"
    log:
        log = "results/cupcake/collapsed/{sample}/out.log"
    shell:
        """
        set +u; source activate SQANTI3.env
        collapse_isoforms_by_sam.py --input {input.fsa} -s {input.sam} \
            --dun-merge-5-shorter -o {params.out} &> {log}
        conda deactivate
        """

rule compress_collapsed_gff:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff"
    output:
        gtf = "results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz",
        tbi = "results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz.tbi",
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule get_abundance_post_collapse: # 获取每个isoform上的FLNC的数量
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        csv = "results/isoseq/polished/{sample}.cluster_report.csv"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.read_stat.txt"
    params:
        prefix = "results/cupcake/collapsed/{sample}/out.collapsed"
    shell:
        """
        set +u; source activate SQANTI3.env
        get_abundance_post_collapse.py {params.prefix} {input.csv} &> /dev/null
        conda deactivate
        """

# Filter away 5' degraded isoforms
rule filter_away_subset:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        abd = "results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.rep.fa",
        out3 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.abundance.txt"
    params:
        prefix = "results/cupcake/collapsed/{sample}/out.collapsed"
    shell:
        """
        set +u; source activate SQANTI3.env
        filter_away_subset.py {params.prefix}
        conda deactivate
        """
    
rule compress_filtered_gff:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff"
    output:
        gtf = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz",
        tbi = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# 使用cupcake合并需要collapsed的其他结果文件，适用范围较小
rule chain:
    input:
        in1 = "results/cupcake/collapsed/Adult_Female/out.collapsed.filtered.gff",
        in2 = "results/cupcake/collapsed/Adult_Male/out.collapsed.filtered.gff",
        in3 = "results/cupcake/collapsed/Juvenile/out.collapsed.filtered.gff",
    output:
        cfg = "results/cupcake/chain/config.txt",
        out1 = "results/cupcake/chain/all_samples.chained.gff",
        out2 = "results/cupcake/chain/all_samples.chained_ids.txt",
        out3 = "results/cupcake/chain/all_samples.chained_count.txt",
    log:
        log = "results/cupcake/chain/chain.log"
    threads:
        8
    shell:
        """
        cat > {output.cfg} << EOF
SAMPLE=Adult_Female;results/cupcake/collapsed/Adult_Female
SAMPLE=Adult_Male;results/cupcake/collapsed/Adult_Male
SAMPLE=Juvenile;results/cupcake/collapsed/Juvenile
GROUP_FILENAME=out.collapsed.group.txt
GFF_FILENAME=out.collapsed.filtered.gff
COUNT_FILENAME=out.collapsed.filtered.abundance.txt
FASTA_FILENAME=out.collapsed.filtered.rep.fa
EOF
        set +u; source activate SQANTI3.env
        chain_samples.py --cpus {threads} {output.cfg} count_fl &> {log}
        conda deactivate
        mv all_samples.chained.gff all_samples.chained_ids.txt all_samples.chained_count.txt `dirname {output.cfg}`
        """

rule compress_chain_gff:
    input:
        gff = "results/cupcake/chain/all_samples.chained.gff"
    output:
        gtf = "results/cupcake/chain/all_samples.chained.sorted.gtf.gz",
        tbi = "results/cupcake/chain/all_samples.chained.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# Fusion gene

rule fusion_finder:
    input:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta",
        sam = "results/cupcake/mapping/{sample}.sorted.sam",
        csv = "results/isoseq/polished/{sample}.cluster_report.csv"
    output:
        out1 = "results/cupcake/fusion/{sample}/out.abundance.txt",
        out2 = "results/cupcake/fusion/{sample}/out.abundance.txt.bak",
        out3 = "results/cupcake/fusion/{sample}/out.gff",
        out4 = "results/cupcake/fusion/{sample}/out.group.txt",
        out5 = "results/cupcake/fusion/{sample}/out.read_stat.txt",
        out6 = "results/cupcake/fusion/{sample}/out.rep.fa",
    params:
        prefix = "results/cupcake/fusion/{sample}/out"
    log:
        log = "results/cupcake/fusion/{sample}/out.fusion.log"
    shell:
        """
        set +u; source activate SQANTI3.env
        fusion_finder.py --input {input.fsa} -s {input.sam} --cluster_report {input.csv} \
            -o {params.prefix} --min_locus_coverage_bp 500 -d 1000000 &> {log}
        conda deactivate
        mv {output.out1} {output.out2}
        grep -v '#' {output.out2} > {output.out1}
        """

rule sqanti3_qc_for_fusion:
    input:
        gff = "results/cupcake/fusion/{sample}/out.gff",
        gtf = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta",
    output:
        out1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
        out2 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.fasta",
        out3 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.genePred",
        out4 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf",
        out5 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf.cds.gff",
        out6 = "results/cupcake/fusion.sqanti3/{sample}/out_junctions.txt",
        out7 = "results/cupcake/fusion.sqanti3/{sample}/out.params.txt",
        out8 = "results/cupcake/fusion.sqanti3/{sample}/out_sqanti_report.pdf",
        out9 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred"
    log:
        log = "results/cupcake/fusion.sqanti3/{sample}.log"
    params:
        prefix = "results/cupcake/fusion.sqanti3/{sample}"
    shell:
        """
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py --gtf {input.gff} {input.gtf} {input.fsa} \
            --is_fusion -d {params.prefix} &> {log}
        conda deactivate
        """

rule fusion_collate_info:
    input:
        in1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
        in2 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred",
        fsa = "data/genome/genome.fasta"
    output:
        out1 = "results/cupcake/fusion/{sample}/out.annotated.txt",
        out2 = "results/cupcake/fusion/{sample}/out.annotated_ignored.txt"
    params:
        prefix = "results/cupcake/fusion/{sample}/out"
    shell:
        """
        set +u; source activate SQANTI3.env
        fusion_collate_info.py {params.prefix} {input.in1} {input.in2} --genome {input.fsa}
        conda deactivate
        """

