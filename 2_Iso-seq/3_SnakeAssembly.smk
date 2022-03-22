#!/usr/bin/env snakemake
configfile: "config.yaml"
import yaml
samples = config["samples"]
with open("../1_ngs/config.yaml") as f:
    rnaseq_list = yaml.load(f, Loader=yaml.BaseLoader)["samples"]
threads = 20
outdir = "results/assembly"
# 流程：比对、注释、合并、过滤内部启动isoforms

rule all:
    input:
        expand(outdir + "/cupcake/collapsed/{sample}/out.collapsed.gff", sample=samples),
        expand(outdir + "/cupcake/collapsed/{sample}/out.collapsed.abundance.txt", sample=samples),
        expand(outdir + "/cupcake/collapsed/{sample}/out.collapsed.filtered.gff", sample=samples),
        outdir + "/cupcake/chain/all_samples.chained.gff",
        # expand(outdir + "/cupcake.collapsed/{sample}.filtered.gtf.gz", sample=samples),
        # outdir + "/cupcake.chain/all_samples.chained.gtf.gz",
        outdir + "/cupcake/chain/all_samples.chained.corrected.gtf.gz",
        outdir + "/cupcake.gtf.gz",
        # expand(outdir + "cupcake.fusion/{sample}/out.gff", sample=samples),
        # expand(outdir + "cupcake.fusion.sqanti3/{sample}/out_classification.txt", sample=samples),
        # expand(outdir + "cupcake.fusion/{sample}/out.annotated.txt", sample=samples),

# Cupcake

rule collapse:
    input:
        fasta = "results/isoseq/polished/{sample}.hq.fasta.gz", # 序列名称里面包含有FLNC的数量
        bam = "results/mapping/clean/{sample}.bam"
    output:
        fasta = temp(outdir + "/cupcake/collapsed/{sample}/out.hq.fasta"),
        gff = outdir + "/cupcake/collapsed/{sample}/out.collapsed.gff",
        # out2 = outdir + "/cupcake.collapsed/{sample}/out.collapsed.gff.unfuzzy",
        # out3 = outdir + "/cupcake.collapsed/{sample}/out.collapsed.group.txt",
        # out4 = outdir + "/cupcake.collapsed/{sample}/out.collapsed.group.txt.unfuzzy",
        # out5 = outdir + "/cupcake.collapsed/{sample}/out.ignored_ids.txt",
        # out6 = outdir + "/cupcake.collapsed/{sample}/out.collapsed.rep.fa",
        gtf = outdir + "/cupcake/collapsed/{sample}.gtf.gz",
        tbi = outdir + "/cupcake/collapsed/{sample}.gtf.gz.tbi"
    params:
        prefix = outdir + "/cupcake/collapsed/{sample}/out",
        script = "../software/cDNA_Cupcake-master/cupcake/tofu/collapse_isoforms_by_sam.py"
    log:
        log = outdir + "/cupcake/collapsed/{sample}/collapse.log"
    shell:
        """
        gzip -d -c {input.fasta} > {output.fasta}
        python {params.script} --max_3_diff 5 --input {output.fasta} --bam {input.bam} --dun-merge-5-shorter -o {params.prefix} &> {log}
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule get_abundance_post_collapse: # 获取每个isoform上的FLNC的数量
    input:
        csv = "results/isoseq/polished/{sample}.cluster_report.csv",
        gff = rules.collapse.output.gff
    output:
        txt = outdir + "/cupcake/collapsed/{sample}/out.collapsed.abundance.txt",
        # out2 = outdir + "/cupcake.collapsed/{sample}/out.collapsed.read_stat.txt"       
    params:
        prefix = outdir + "/cupcake/collapsed/{sample}/out.collapsed",
        script = "../software/cDNA_Cupcake-master/cupcake/tofu/get_abundance_post_collapse.py"
    shell:
        """
        python {params.script} {params.prefix} {input.csv} &> /dev/null
        """

# Filter away 5' degraded isoforms
rule filter_away_subset:
    input:
        gff = rules.collapse.output.gff,
        txt = rules.get_abundance_post_collapse.output.txt
    output:
        gff = outdir + "/cupcake/collapsed/{sample}/out.collapsed.filtered.gff",
        # out2 = outdir + "/cupcake/collapsed/{sample}/out.collapsed.filtered.rep.fa",
        # out3 = outdir + "/cupcake/collapsed/{sample}/out.collapsed.filtered.abundance.txt",
        gtf = outdir + "/cupcake/collapsed/{sample}.filtered.gtf.gz",
        tbi = outdir + "/cupcake/collapsed/{sample}.filtered.gtf.gz.tbi",
    params:
        prefix = outdir + "/cupcake/collapsed/{sample}/out.collapsed",
        script = "../software/cDNA_Cupcake-master/cupcake/tofu/filter_away_subset.py"
    shell:
        """
        python {params.script} {params.prefix}
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# 使用cupcake合并需要collapsed的其他结果文件，适用范围较小

rule chain:
    input:
        in1 = outdir + "/cupcake/collapsed/Ad_Fe/out.collapsed.filtered.gff",
        in2 = outdir + "/cupcake/collapsed/Ad_Ma/out.collapsed.filtered.gff",
        in3 = outdir + "/cupcake/collapsed/Ju_Mi/out.collapsed.filtered.gff",
    output:
        cfg = temp(outdir + "/cupcake/chain/config.txt"),
        gff = outdir + "/cupcake/chain/all_samples.chained.gff",
        # out2 = outdir + "/cupcake/chain/all_samples.chained_ids.txt",
        # out3 = outdir + "/cupcake/chain/all_samples.chained_count.txt",
        gtf = outdir + "/cupcake/chain/all_samples.chained.gtf.gz",
        tbi = outdir + "/cupcake/chain/all_samples.chained.gtf.gz.tbi"
    log:
        log = outdir + "/cupcake/chain/chain.log"
    params:
        script = "../software/cDNA_Cupcake-master/cupcake/tofu/counting/chain_samples.py"
    threads:
        8
    shell:
        """
        cat > {output.cfg} << EOF
SAMPLE=Ad_Fe;results/assembly/cupcake/collapsed/Ad_Fe
SAMPLE=Ad_Ma;results/assembly/cupcake/collapsed/Ad_Ma
SAMPLE=Ju_Mi;results/assembly/cupcake/collapsed/Ju_Mi
GROUP_FILENAME=out.collapsed.group.txt
GFF_FILENAME=out.collapsed.filtered.gff
COUNT_FILENAME=out.collapsed.filtered.abundance.txt
FASTA_FILENAME=out.collapsed.filtered.rep.fa
EOF
        python {params.script} --cpus {threads} {output.cfg} count_fl &> {log}
        mv all_samples.chained.gff all_samples.chained_ids.txt all_samples.chained_count.txt `dirname {output.cfg}`
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule correct_gene_id_of_cupcake_chain_gtf:
    input:
        gff = rules.chain.output.gff
    output:
        gff = outdir + "/cupcake/chain/all_samples.chained.corrected.gff",
        gtf = outdir + "/cupcake/chain/all_samples.chained.corrected.gtf.gz",
        tbi = outdir + "/cupcake/chain/all_samples.chained.corrected.gtf.gz.tbi"
    shell:
        """
        ./scripts/correct_gene_id_of_cupcake_chain_gtf.py {input.gff} {output.gff}
        bedtools sort -i {output.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# Filtered (线粒体上的isoforms和长度小于200的isoforms也会被过滤掉)

rule filter_internal_primer_isoform:
    input:
        gtf = rules.correct_gene_id_of_cupcake_chain_gtf.output.gtf,
        fasta = config["genome_fasta"],
    output:
        gtf = outdir + "/cupcake.gtf.gz",
        tbi = outdir + "/cupcake.gtf.gz.tbi"
    shell:
        """
        ./scripts/filter_internal_primer_isoform.py {input.gtf} {input.fasta} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# Fusion gene

# rule fusion_finder:
#     input:
#         fsa = "results/cupcake/mapping/{sample}.hq.fasta",
#         sam = "results/cupcake/mapping/{sample}.sorted.sam",
#         csv = "results/isoseq/polished/{sample}.cluster_report.csv"
#     output:
#         out1 = "results/cupcake/fusion/{sample}/out.abundance.txt",
#         out2 = "results/cupcake/fusion/{sample}/out.abundance.txt.bak",
#         out3 = "results/cupcake/fusion/{sample}/out.gff",
#         out4 = "results/cupcake/fusion/{sample}/out.group.txt",
#         out5 = "results/cupcake/fusion/{sample}/out.read_stat.txt",
#         out6 = "results/cupcake/fusion/{sample}/out.rep.fa",
#     params:
#         prefix = "results/cupcake/fusion/{sample}/out"
#     log:
#         log = "results/cupcake/fusion/{sample}/out.fusion.log"
#     shell:
#         """
#         set +u; source activate SQANTI3.env
#         fusion_finder.py --input {input.fsa} -s {input.sam} --cluster_report {input.csv} \
#             -o {params.prefix} --min_locus_coverage_bp 500 -d 1000000 &> {log}
#         conda deactivate
#         mv {output.out1} {output.out2}
#         grep -v '#' {output.out2} > {output.out1}
#         """

# rule sqanti3_qc_for_fusion:
#     input:
#         gff = "results/cupcake/fusion/{sample}/out.gff",
#         gtf = "data/genome/annotation.gtf",
#         fsa = "data/genome/genome.fasta",
#     output:
#         out1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
#         out2 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.fasta",
#         out3 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.genePred",
#         out4 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf",
#         out5 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf.cds.gff",
#         out6 = "results/cupcake/fusion.sqanti3/{sample}/out_junctions.txt",
#         out7 = "results/cupcake/fusion.sqanti3/{sample}/out.params.txt",
#         out8 = "results/cupcake/fusion.sqanti3/{sample}/out_sqanti_report.pdf",
#         out9 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred"
#     log:
#         log = "results/cupcake/fusion.sqanti3/{sample}.log"
#     params:
#         prefix = "results/cupcake/fusion.sqanti3/{sample}"
#     shell:
#         """
#         set +u; source activate SQANTI3.env
#         ~/software/SQANTI3/sqanti3_qc.py --gtf {input.gff} {input.gtf} {input.fsa} \
#             --is_fusion -d {params.prefix} &> {log}
#         conda deactivate
#         """

# rule fusion_collate_info:
#     input:
#         in1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
#         in2 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred",
#         fsa = "data/genome/genome.fasta"
#     output:
#         out1 = "results/cupcake/fusion/{sample}/out.annotated.txt",
#         out2 = "results/cupcake/fusion/{sample}/out.annotated_ignored.txt"
#     params:
#         prefix = "results/cupcake/fusion/{sample}/out"
#     shell:
#         """
#         set +u; source activate SQANTI3.env
#         fusion_collate_info.py {params.prefix} {input.in1} {input.in2} --genome {input.fsa}
#         conda deactivate
#         """

# Common rules

rule bam_stats:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input.bam} > {output.txt}
        """
