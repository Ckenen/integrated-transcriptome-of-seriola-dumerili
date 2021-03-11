#!/usr/bin/env python
configfile: "config.yaml"
samples = config["samples"]
threads = 80
outdir = "results/cupcake"


rule all:
    input:
        outdir + "/genome.mmi",
        expand(outdir + "/mapping/{sample}.hq.fasta", sample=samples),
        expand(outdir + "/mapping/{sample}.sam", sample=samples),
        expand(outdir + "/mapping/{sample}.sorted.bam", sample=samples),
        expand(outdir + "/mapping/{sample}.sorted.sam", sample=samples),
        expand(outdir + "/collapsed/{sample}/out.collapsed.gff", sample=samples),
        expand(outdir + "/collapsed/{sample}/out.collapsed.abundance.txt", sample=samples),
        # expand(outdir + "/collapsed/{sample}/out.collapsed.filtered.gff", sample=samples),
        # outdir + "/chained",
        expand(outdir + "/fusion/{sample}/out.gff", sample=samples),
        expand(outdir + "/fusion.sqanti3/{sample}", sample=samples),
        expand("results/cupcake/fusion/{sample}/out.annotated.txt", sample=samples),

# Minimap2

rule minimap2_index:
    input:
        fa = "data/genome/genome.fasta"
    output:
        mmi = outdir + "/genome.mmi"
    threads: 
        threads
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fa}
        """

rule unzip_hq_fasta:
    input:
        fa = "results/isoseq/polished/{sample}.hq.fasta.gz"
    output:
        fa = outdir + "/mapping/{sample}.hq.fasta"
    threads: 
        4
    shell:
        """
        pigz -p {threads} -d -c {input.fa} > {output.fa}
        """

rule minimap2_mapping:
    input:
        mmi = rules.minimap2_index.output.mmi,
        fa = rules.unzip_hq_fasta.output.fa
    output:
        sam = outdir + "/mapping/{sample}.sam"
    log:
        outdir + "/mapping/{sample}.log"
    threads:
        threads
    shell:
        """
        minimap2 -ax splice -t {threads} --secondary=no -C5 {input.mmi} {input.fa} > {output.sam} 2> {log}
        """

rule sort_sam:
    input:
        sam = rules.minimap2_mapping.output.sam
    output:
        sam = outdir + "/mapping/{sample}.sorted.sam"
    shell:
        """
        sort -k 3,3 -k 4,4n {input.sam} > {output.sam}
        """

rule sam_to_bam:
    input:
        sam = rules.minimap2_mapping.output.sam
    output:
        tmp = temp(outdir + "/mapping/{sample}.tmp.bam"),
        bam = outdir + "/mapping/{sample}.sorted.bam",
        bai = outdir + "/mapping/{sample}.sorted.bam.bai"
    shell:
        """
        samtools view -b {input.sam} > {output.tmp}
        bamtools sort -in {output.tmp} -out {output.bam}
        bamtools index -in {output.bam}
        """

# Collapse redundant isoforms (has genome)

rule collapse_isoforms_by_sam:
    input:
        fa = rules.unzip_hq_fasta.output.fa,
        sam = rules.sort_sam.output.sam
    output:
        gff1 = outdir + "/collapsed/{sample}/out.collapsed.gff",
        gff2 = outdir + "/collapsed/{sample}/out.collapsed.gff.unfuzzy",
        txt1 = outdir + "/collapsed/{sample}/out.collapsed.group.txt",
        txt2 = outdir + "/collapsed/{sample}/out.collapsed.group.txt.unfuzzy",
        fasa = outdir + "/collapsed/{sample}/out.collapsed.rep.fa",
        txt3 = outdir + "/collapsed/{sample}/out.ignored_ids.txt"
    log:
        outdir + "/collapsed/{sample}.log"
    params:
        prefix = outdir + "/collapsed/{sample}/out"
    shell:
        """
        set +u
        source activate SQANTI3.env
        collapse_isoforms_by_sam.py --input {input.fa}  -s {input.sam} --dun-merge-5-shorter -o {params.prefix} &> {log}
        conda deactivate
        """

# Obtain associated count information

rule get_abundance_post_collapse:
    input:
        gff = outdir + "/collapsed/{sample}/out.collapsed.gff",
        csv = "results/isoseq/polished/{sample}.cluster_report.csv"
    output:
        txt1 = outdir + "/collapsed/{sample}/out.collapsed.abundance.txt",
        txt2 = outdir + "/collapsed/{sample}/out.collapsed.read_stat.txt"
    params:
        prefix = outdir + "/collapsed/{sample}/out.collapsed"
    shell:
        """
        set +u
        source activate SQANTI3.env
        get_abundance_post_collapse.py {params.prefix} {input.csv} > /dev/null 2>&1
        conda deactivate
        """

# Filter collapse results by minimum FL count support
# Does not filter

# Filter away 5' degraded isoforms

rule filter_away_subset:
    input:
        gff = outdir + "/collapsed/{sample}/out.collapsed.gff",
        txt = outdir + "/collapsed/{sample}/out.collapsed.abundance.txt"
    output:
        gff = outdir + "/collapsed/{sample}/out.collapsed.filtered.gff",
        fsa = outdir + "/collapsed/{sample}/out.collapsed.filtered.rep.fa",
        txt = outdir + "/collapsed/{sample}/out.collapsed.filtered.abundance.txt"
    params:
        prefix = outdir + "/collapsed/{sample}/out.collapsed"
    shell:
        """
        set +u
        source activate SQANTI3.env
        filter_away_subset.py {params.prefix}
        conda deactivate
        """

# Chain samples together

rule chain_samples:
    input:
        gff1 = outdir + "/collapsed/Adult_Female/out.collapsed.filtered.gff",
        gff2 = outdir + "/collapsed/Adult_Male/out.collapsed.filtered.gff",
        gff3 = outdir + "/collapsed/Juvenile/out.collapsed.filtered.gff",
    output:
        cfg = "cupcake.txt",
        out = directory(outdir + "/chained")
    log:
        outdir + "/chained.log"
    params:
        out1 = "all_samples.chained_ids.txt",
        out2 = "all_samples.chained_count.txt",
        out3 = "all_samples.chained.gff",
        out4 = "all_samples.chained.rep.fa"
    shell:
        """( 
        echo "SAMPLE=Adult_Female;results/cupcake/collapsed/Adult_Female"
        echo "SAMPLE=Adult_Male;results/cupcake/collapsed/Adult_Male"
        echo "SAMPLE=Juvenile;results/cupcake/collapsed/Juvenile"
        echo "GROUP_FILENAME=out.collapsed.group.txt" 
        echo "GFF_FILENAME=out.collapsed.filtered.gff"
        echo "COUNT_FILENAME=out.collapsed.filtered.abundance.txt"
        echo "FASTA_FILENAME=out.collapsed.filtered.rep.fa" ) > {output.cfg}
        set +u
        source activate SQANTI3.env
        chain_samples.py {output.cfg} count_fl &> {log}
        conda deactivate
        mkdir -p {output.out}
        mv {params.out1} {params.out2} {params.out3} {params.out4} {output.out}
        """

# rule chained_sqanti3:
#     input:
#     output:
#     log:
#     shell:
#         """

#         """

# Finding fusion genes

rule fusion_finder:
    input:
        fsa = "results/cupcake/mapping/Adult_Female.hq.fasta",
        sam = "results/cupcake/mapping/Adult_Female.sorted.sam",
        csv = "results/isoseq/polished/Adult_Female.cluster_report.csv"
    output:
        gff = "results/cupcake/fusion/{sample}/out.gff",
        fsa = "results/cupcake/fusion/{sample}/out.rep.fa",
        ab1 = "results/cupcake/fusion/{sample}/out.abundance.txt",
        ab2 = "results/cupcake/fusion/{sample}/out.abundance.txt.bak",
        sta = "results/cupcake/fusion/{sample}/out.read_stat.txt",
    log:
        "results/cupcake/fusion/{sample}/out.fusion.log"
    params:
        prefix = "results/cupcake/fusion/{sample}/out"
    shell:
        """
        set +u
        source activate SQANTI3.env
        fusion_finder.py --input {input.fsa} -s {input.sam} --cluster_report {input.csv} \
            -o {params.prefix} --min_locus_coverage_bp 500 -d 1000000 &> {log}
        mv {output.ab1} {output.ab2}
        grep -v '#' {output.ab2} > {output.ab1}
        conda deactivate
        """

# SQANTI3

rule fusion_sqanti3:
    input:
        gff = outdir + "/fusion/{sample}/out.gff",
        gtf = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta"
    output:
        directory(outdir + "/fusion.sqanti3/{sample}")
    log:
        outdir + "/fusion.sqanti3/{sample}.log"
    shell:
        """
        set +u
        source activate SQANTI3.env
        mkdir -p {output}
        ~/software/SQANTI3/sqanti3_qc.py --gtf {input.gff} {input.gtf} {input.fsa} --is_fusion -d {output} &> {log}
        conda deactivate
        """

rule fusion_collate_info:
    input:
        gno = "data/genome/genome.fasta",
        gff = "results/cupcake/fusion/{sample}/out.gff",
        fsa = "results/cupcake/fusion/{sample}/out.rep.fa",
        txt = "results/cupcake/fusion/{sample}/out.abundance.txt",
        sqc = "results/cupcake/fusion.sqanti3/{sample}"
    output:
        txt = "results/cupcake/fusion/{sample}/out.annotated.txt"
    params:
        prefix = "results/cupcake/fusion/{sample}/out",
        sqc = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
        gpr = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred"
    shell:
        """
        set +u
        source activate SQANTI3.env
        fusion_collate_info.py {params.prefix} {params.sqc} {params.gpr} --genome {input.gno}
        conda deactivate
        """
        
# Common 

rule pbindex:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.pbi"
    shell:
        """
        pbindex {input}
        """

rule bam_stats:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """
