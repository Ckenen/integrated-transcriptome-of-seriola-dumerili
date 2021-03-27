#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = config["threads"]
indir = "data/datasets"
outdir = "results/mapping"

rule all:
    input:
        outdir + "/star_index",
        expand(outdir + "/mapped.1st/{sample}", sample=samples),
        outdir + "/all.SJ.out.tab",
        expand(outdir + "/mapped.2nd/{sample}", sample=samples),
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/filtered/{sample}.bam.bai", sample=samples),
        expand(outdir + "/filtered/{sample}.infer.txt", sample=samples),
        expand(outdir + "/filtered/{sample}.stats", sample=samples),
        # expand(outdir + "/filtered/{sample}.fam", sample=samples),
        # expand(outdir + "/filtered/{sample}.bw", sample=samples),
        # expand(outdir + "/uniq/{sample}.bam", sample=samples),
        # expand(outdir + "/uniq/{sample}.bam.bai", sample=samples),
        # expand(outdir + "/uniq/{sample}.stats", sample=samples),
        # expand(outdir + "/uniq/{sample}.fam", sample=samples),
        # expand(outdir + "/uniq/{sample}.bw", sample=samples),


# Mapping to genome

rule star_index:
    input:
        fas = "data/genome/genome.fasta",
        gtf = "data/genome/annotation.gtf"
    output:
        directory(outdir + "/star_index")
    log:
        outdir + "/star_index.log"
    threads:
        threads
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeSAindexNbases 10 --genomeDir {output} \
            --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fas} &> {log}
        """

rule star_mapping_1st:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz",
        idx = rules.star_index.output
    output:
        directory(outdir + "/mapped.1st/{sample}")
    log:
        outdir + "/mapped.1st/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.1st/{sample}",
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads --outSAMtype None  --outFilterMismatchNoverLmax 0.2 --outFilterMismatchNmax 20 \
            --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule merge_splice_junction:
    input:
        tabs = expand(outdir + "/mapped.1st/{sample}", sample=samples)
    output:
        outdir + "/all.SJ.out.tab"
    shell:
        """
        for tab in {input.tabs}; do cat $tab/SJ.out.tab; done | awk '$7>=5' | sort -k1,1 -k2,2n -k3,3n > {output}
        """

rule star_mapping_2nd:
    input:
        fq1 = indir + "/{sample}_R1.fastq.gz",
        fq2 = indir + "/{sample}_R2.fastq.gz", 
        idx = outdir + "/star_index",
        spl = rules.merge_splice_junction.output
    output:
        directory(outdir + "/mapped.2nd/{sample}")
    log:
        outdir + "/mapped.2nd/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.2nd/{sample}"
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 5 \
            --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS --outFilterMismatchNoverLmax 0.1 \
            --outFilterMismatchNmax 10 --limitBAMsortRAM 10000000000 \
            --quantMode TranscriptomeSAM --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} --sjdbFileChrStartEnd {input.spl} --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

# Filtered

rule filter_alignment:
    input:
        bam = rules.star_mapping_2nd.output,
    output:
        bam = outdir + "/filtered/{sample}.bam"
    shell:
        """
        bamtools filter -in {input.bam}/Aligned.sortedByCoord.out.bam -out {output} -tag NH:1 -isProperPair true -isPrimaryAlignment true
        """

rule infer_experiment:
    input:
        bam = rules.filter_alignment.output,
        bed = "data/genome/transcript.bed"
    output:
        txt = outdir + "/filtered/{sample}.infer.txt"
    log:
        outdir + "/filtered/{sample}.infer.log"
    shell:
        """
        set +u
        source activate py27
        infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        conda deactivate 
        """

# # Uniq

rule remove_duplicate:
    input:
        bam = rules.filter_alignment.output.bam,
        bai = rules.filter_alignment.output.bam + ".bai"
    output:
        bam = outdir + "/uniq/{sample}.bam",
        txt = outdir + "/uniq/{sample}_metrics.txt"
    log:
        outdir + "/uniq/{sample}.log"
    threads:
        10
    shell:
        """
        set +u
        source activate picard
        picard MarkDuplicates REMOVE_DUPLICATES=true I={input.bam} O={output.bam} M={output.txt} &> {log}
        conda deactivate
        """

# Common rules

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        bamtools index -in {input}
        """

rule bam_stat:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """

rule bam_to_fam:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai"
    output:
        fam = "{prefix}.fam"
    threads:
        10
    shell:
        """
        famtools.py build {input.bam} {output.fam}
        """

rule fam_to_bigwig:
    input:
        fam = "{prefix}.fam"
    output:
        bw = "{prefix}.bw",
    shell:
        """
        famtools.py convert -f bigwig -s 2 -S {input.fam} {output.bw}
        """

# rule fam_to_bigwig:
#     input:
#         fam = "{prefix}.fam"
#     output:
#         bw1 = "{prefix}.bw",
#         bw2 = "{prefix}.+.bw",
#         bw3 = "{prefix}.-.bw",
#         bw4 = "{prefix}.normalized.bw",
#         bw5 = "{prefix}.normalized.+.bw",
#         bw6 = "{prefix}.normalized.-.bw",
#     params:
#         prefix = "{prefix}"
#     shell:
#         """
#         ./scripts/fam_to_bigwig.sh --rf {input} {params.prefix}
#         """
