#!/usr/bin/env snakemake
configfile: "config.yaml"
samples = config["samples"]
# GENOME_FASTA = "../common/ncbi/serDum.ncbi.fasta"
threads = 20
outdir = "results"

# Reference:
# Data pre-processing for variant discovery (https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
# Germline short variant discovery (SNPs + Indels) (https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

rule all:
    input:
        # expand(outdir + "/sra/{sample}.sra", sample=samples),
        # expand(outdir + "/fastq/{sample}_{r}.fastq.gz", sample=samples, r=["1", "2"]),
        # expand(outdir + "/fastq/{sample}_{r}_fastqc.html", sample=samples, r=["1", "2"]),
        # expand(outdir + "/trimmomatic/{sample}_1.fastq.gz", sample=samples),
        outdir + "/bwa/index",
        expand(outdir + "/bwa/mapped/{sample}.bam", sample=samples),
        expand(outdir + "/bwa/mapped/{sample}.bam.bai", sample=samples),
        expand(outdir + "/bwa/mapped/{sample}.stats", sample=samples),
        # expand("results/bwa/mapped/{sample}.rc.txt", sample=samples),
        expand(outdir + "/bwa/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/bwa/filtered/{sample}.bam.bai", sample=samples),
        expand(outdir + "/bwa/filtered/{sample}.stats", sample=samples),
        expand(outdir + "/gatk/markdup/{sample}.bam", sample=samples),
        expand(outdir + "/gatk/markdup/{sample}.bam.bai", sample=samples),
        expand(outdir + "/gatk/markdup/{sample}.stats", sample=samples),
        # expand(outdir + "/gatk/base_recalibrator/{sample}.txt", sample=samples),
        # expand(outdir + "/rmdup/{sample}.bam", sample=samples),
        # expand(outdir + "/rmdup/{sample}.bam.bai", sample=samples),
        # expand(outdir + "/rmdup/{sample}.stats", sample=samples),
        expand(outdir + "/gatk/haplotype/{sample}.gvcf.gz", sample=samples),
        outdir + "/gatk/combined/all.gvcf.gz",

        # expand("results/haplotype/splited/{sample}", sample=samples),
        # expand("results/haplotype/single/{sample}.gvcf.gz", sample=samples),
        # "results/haplotype/merged/all.gvcf.gz",
        # "results/haplotype/genotype.vcf.gz",
        # "results/haplotype/genotype.snp.vcf.gz",
        # "results/haplotype/genotype.snp.marked.vcf.gz",
        # "results/haplotype/genotype.snp.marked.clean.vcf.gz",

# Prepare 

rule prefetch:
    output:
        sra = outdir + "/sra/{sample}.sra"
    log:
        outdir + "/sra/{sample}.log"
    shell:
        """
        prefetch --max-size 50000000 -o {output.sra} {wildcards.sample} &> {log}
        """

rule fastq_dump:
    input:
        sra = rules.prefetch.output.sra
    output:
        tmp1 = temp(outdir + "/fastq/{sample}_1.fastq"),
        tmp2 = temp(outdir + "/fastq/{sample}_2.fastq"),
        fq1 = outdir + "/fastq/{sample}_1.fastq.gz",
        fq2 = outdir + "/fastq/{sample}_2.fastq.gz"
    log:
        outdir + "/fastq/{sample}.log"
    params:
        odir = outdir + "/fastq"
    shell:
        """(
        fastq-dump --split-3 --outdir {params.odir} {input.sra}
        gzip -c {output.tmp1} > {output.fq1}
        gzip -c {output.tmp2} > {output.fq2} ) &> {log}
        """

rule fastqc:
    input:
        fastq = "{prefix}.fastq.gz"
    output:
        html = "{prefix}_fastqc.html",
        zf = "{prefix}_fastqc.zip"
    log:
        "{prefix}_fastqc.log"
    shell:
        """
        fastqc --outdir `dirname {output.html}` {input.fastq} &> {log}
        """

rule trimmomatic:
    input:
        read1 = rules.fastq_dump.output.fq1,
        read2 = rules.fastq_dump.output.fq2,
        fasta = "/home/chenzonggui/software/Trimmomatic-main/adapters/TruSeq2-PE.fa"
    output:
        pair1 = outdir + "/trimmomatic/{sample}_1.fastq.gz",
        unpair1 = outdir + "/trimmomatic/{sample}_1.unpair.fastq.gz",
        pair2 = outdir + "/trimmomatic/{sample}_2.fastq.gz",
        unpair2 = outdir + "/trimmomatic/{sample}_2.unpair.fastq.gz"
    log:
        log = outdir + "/trimmomatic/{sample}.log"
    threads:
        threads
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 {input.read1} {input.read2} \
            {output.pair1} {output.unpair1} {output.pair2} {output.unpair2} \
            ILLUMINACLIP:{input.fasta}:2:30:10:1:TRUE LEADING:20 \
            TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36 &> {log}
        """

# Mapping 

rule bwa_index:
    input:
        fasta = config["genome"]
    output:
        directory("results/bwa/index")
    log:
        "results/bwa/index.log"
    shell:
        """
        mkdir {output}
        bwa index -a is -p {output}/ref {input.fasta} &> {log}
        """

rule bwa_align:
    input:
        index = rules.bwa_index.output,
        read1 = rules.fastq_dump.output.fq1,
        read2 = rules.fastq_dump.output.fq2
    output:
        tmp = outdir + temp("/bwa/mapped/{sample}.tmp.bam"),
        bam = outdir + "/bwa/mapped/{sample}.bam"
    log:
        log = outdir + "/bwa/mapped/{sample}.log"
    params:
        header = "@RG\\tID:{sample}\\tLB:{sample}\\tPL:Illumina\\tPU:{sample}\\tSM:{sample}"
    threads:
        threads
    shell:
        """(
        bwa mem -t {threads} -R "{params.header}" {input.index}/ref {input.read1} {input.read2} \
            | samtools view -b - > {output.tmp}
        samtools sort -@ {threads} -o {output.bam} {output.tmp} ) &> {log}
        """

rule filter_alignments:
    input:
        bam = outdir + "/bwa/mapped/{sample}.bam"
    output:
        bam = outdir + "/bwa/filtered/{sample}.bam"
    shell:
        """
        ./scripts/filter_alignments.py {input.bam} {output.bam}
        """

# GATK

rule mark_duplicates: # MarkDuplicates
    input:
        bam = rules.filter_alignments.output.bam,
        bai = rules.filter_alignments.output.bam + ".bai"
    output:
        bam = outdir + "/gatk/markdup/{sample}.bam",
        txt = outdir + "/gatk/markdup/{sample}_metrics.txt"
    log:
        outdir + "/gatk/markdup/{sample}.log"
    threads:
        threads
    shell:
        """
        gatk MarkDuplicates --TAG_DUPLICATE_SET_MEMBERS true --TAGGING_POLICY All \
            -I {input.bam} -O {output.bam} -M {output.txt} &> {log}
        """

# rule BaseRecalibrator:
#     input:
#         bam = rules.mark_duplicates.output.bam,
#         fasta = config["genome"]
#     output:
#         tbl = outdir + "/gatk/base_recalibrator/{sample}.txt"
#     log:
#         outdir + "/gatk/base_recalibrator/{sample}.log"
#     shell:
#         """
#         gatk BaseRecalibrator -I {input.bam} -R {input.fasta} -O {output.tbl} &> {log}
#         """

# gatk CreateSequenceDictionary -R ../common/ncbi/serDum.ncbi.fasta -O ../common/ncbi/serDum.ncbi.dict

rule call_haplotype: # HaplotypeCaller
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bam + ".bai",
        fasta = config["genome"]
    output:
        gvcf = outdir + "/gatk/haplotype/{sample}.gvcf.gz"
    log:
        log = outdir + "/gatk/haplotype/{sample}.log"
    threads:
        4
    shell:
        """
        gatk --java-options -Xmx4G HaplotypeCaller -I {input.bam} -O {output.gvcf} -R {input.fasta} --emit-ref-confidence GVCF &> {log}
        """

rule combine_gvcf: # CombineGVCFs
    input:
        fasta = config["genome"],
        gvcf_list = [rules.call_haplotype.output.gvcf.format(sample=s) for s in samples]
    output:
        gvcf = outdir + "/gatk/combined/all.gvcf.gz"
    log:
        outdir + "/gatk/combined/all.log"
    params:
        gvcf_list = " ".join("-V results/gatk/haplotype/%s.gvcf.gz" % s for s in samples)
    shell:
        """
        gatk CombineGVCFs -R {input.fasta} -O {output.gvcf} {params.gvcf_list} &> {log}
        """

rule genotyping_gvcf: # GenotypeGVCFs
    input:
        fasta = config["genome"],
        gvcf = rules.combine_gvcf.output.gvcf
    output:
        vcf = outdir + "/gatk/combined/genotyped.vcf.gz"
    log:
        log = outdir + "/gatk/combined/genotyped.log"
    shell:
        """
        gatk GenotypeGVCFs -R {input.fasta} -V {input.gvcf} -O {output.vcf} &> {log}
        """

rule fetch_snp: # SelectVariants
    input:
        vcf = rules.genotyping_gvcf.output.vcf
    output:
        vcf = outdir + "/gatk/combined/genotyped.snp.vcf.gz"
    log:
        outdir + "/gatk/combined/genotyped.snp.log"
    shell:
        """
        gatk SelectVariants -V {input.vcf} -O {output.vcf} --select-type-to-include SNP &> {log}
        """

rule mark_snp: # VariantFiltration
    input:
        vcf = rules.fetch_snp.output.vcf
    output:
        vcf = outdir + "/gatk/combined/genotyped.snp.marked.vcf.gz"
    log:
        outdir + "/gatk/combined/genotyped.snp.marked.log"
    shell:
        """
        gatk VariantFiltration -O {output.vcf} -V {input.vcf} \
            --filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 || SOR > 4.0' \
            --filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 \
            --missing-values-evaluate-as-failing &> {log} 
        """
    
rule filter_snp:
    input:
        vcf = rules.mark_snp.output.vcf
    output:
        vcf = outdir + "/gatk/genotyped.snp.marked.filtered.vcf.gz",
        tbi = outdir + "/gatk/genotyped.snp.marked.filtered.vcf.gz.tbi"
    shell:
        """
        bcftools filter -i 'FILTER=="PASS"' {input.vcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Common rules

rule bam_index:
    input:
        bam = "{prefix}.bam"
    output:
        bai = "{prefix}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule bam_stats:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input.bam} > {output.txt}
        """

rule get_record_count:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.rc.txt"
    shell:
        """
        samtools view {input.bam} | awk '{{print $1}}' | sort \
            | uniq -c | awk -v OFS='\\t' '{{print $2,$1}}' > {output.txt}
        """
