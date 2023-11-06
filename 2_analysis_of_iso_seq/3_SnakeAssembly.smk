#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/filtered"
outdir = "results/assembly"

rule all:
    input:
        expand(outdir + "/stringtie/gtf/{sample}.mp{p}.gtf", sample=samples, p=min_passes_list),
        expand(outdir + "/stringtie/merged/all_samples.mp{p}.gtf", p=min_passes_list),
        expand(outdir + "/tama/collapsed/{sample}.mp{p}", sample=samples, p=min_passes_list),
        expand(outdir + "/tama/merged/all_samples.mp{p}", p=min_passes_list),
        expand(outdir + "/tama/merged_bed/all_samples.mp{p}.bed", p=min_passes_list),
        expand(outdir + "/tama/merged_gtf/all_samples.mp{p}.gtf", p=min_passes_list),
        expand(outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.bed" , p=min_passes_list),
        expand(outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.gtf" , p=min_passes_list),
        

rule stringtie_collapse:
    input:
        bam = indir + "/{sample}.mp{p}.bam"
    output:
        gtf = outdir + "/stringtie/gtf/{sample}.mp{p}.gtf",
        gtf_gz = outdir + "/stringtie/gtf/{sample}.mp{p}.sorted.gtf.gz"
    threads:
        8
    shell:
        """
        stringtie {input.bam} -p {threads} --fr -L | bedtools sort | awk '$7!="."' > {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        """

rule stringtie_merge:
    input:
        gtfs = [outdir + "/stringtie/gtf/%s.mp{p}.gtf" % s for s in samples]
    output:
        gtf = outdir + "/stringtie/merged/all_samples.mp{p}.gtf",
        gtf_gz = outdir + "/stringtie/merged/all_samples.mp{p}.sorted.gtf.gz"
    log:
        outdir + "/stringtie/merged/all_samples.mp{p}.log"
    shell:
        """
        stringtie --merge {input.gtfs} > {output.gtf} 
        grep -v '#' {output.gtf} | sort -k1,1 -k4,4n | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """

rule tama_collapse:
    input:
        bam = indir + "/{sample}.mp{p}.bam",
        fa = GENOME_FASTA
    output:
        sam = temp(outdir + "/tama/collapsed/{sample}.mp{p}.sam"),
        out = directory(outdir + "/tama/collapsed/{sample}.mp{p}")
    log:
        outdir + "/tama/collapsed/{sample}.mp{p}.log"
    params:
        prefix = outdir + "/tama/collapsed/{sample}.mp{p}/{sample}"
    shell:
        """
        samtools view -h {input.bam} > {output.sam}
        mkdir {output.out}
        set +u; source activate TAMA
        python ~/software/tama/tama_collapse.py -s {output.sam} -f {input.fa} -x no_cap -p {params.prefix} &> {log}
        """

rule tama_merge:
    input:
        beddir1 = outdir + "/tama/collapsed/Ad_Ma.mp{p}",
        beddir2 = outdir + "/tama/collapsed/Ad_Fe.mp{p}",
        beddir3 = outdir + "/tama/collapsed/Ju_Mi.mp{p}"
    output:
        txt = outdir + "/tama/merged/all_samples.mp{p}.filelist",
        out = directory(outdir + "/tama/merged/all_samples.mp{p}")
    log:
        outdir + "/tama/merged/all_samples.mp{p}.log"
    params:
        bed1 = outdir + "/tama/collapsed/Ad_Ma.mp{p}/Ad_Ma.bed",
        bed2 = outdir + "/tama/collapsed/Ad_Fe.mp{p}/Ad_Fe.bed",
        bed3 = outdir + "/tama/collapsed/Ju_Mi.mp{p}/Ju_Mi.bed",
        prefix = outdir + "/tama/merged/all_samples.mp{p}/all_samples"
    shell:
        """
        echo -e "{params.bed1}\\tno_cap\\t1,1,1\\tAdMa" >> {output.txt}
        echo -e "{params.bed2}\\tno_cap\\t1,1,1\\tAdFe" >> {output.txt}
        echo -e "{params.bed3}\\tno_cap\\t1,1,1\\tJuMi" >> {output.txt}
        set +u; source activate TAMA; mkdir {output.out}
        python ~/software/tama/tama_merge.py -f {output.txt} -p {params.prefix} &> {log}
        """

rule filter_zero_length_exon_isoform:
    input:
        beddir = rules.tama_merge.output.out
    output:
        bed = outdir + "/tama/merged_bed/all_samples.mp{p}.bed",
        bed_gz = outdir + "/tama/merged_bed/all_samples.mp{p}.sorted.bed.gz"
    shell:
        """
        cat {input.beddir}/all_samples.bed | ./scripts/filter_zero_length_exon_isoform.py > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """

rule bed_to_gtf:
    input:
        bed = rules.filter_zero_length_exon_isoform.output.bed
    output:
        gtf = outdir + "/tama/merged_gtf/all_samples.mp{p}.gtf",
        gtf_gz = outdir + "/tama/merged_gtf/all_samples.mp{p}.sorted.gtf.gz"
    shell:
        """
        set +u; source activate TAMA
        python ~/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {input.bed} {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """

rule filter_internal_primer_isoform:
    input:
        bed = rules.filter_zero_length_exon_isoform.output.bed,
        fasta = GENOME_FASTA
    output:
        bed = outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.bed",
        bed_gz = outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.sorted.bed.gz",
        gtf = outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.gtf",
        gtf_gz = outdir + "/tama/filtered_internal_primer/all_samples.mp{p}.sorted.gtf.gz"
    shell:
        """
        ./scripts/filter_internal_primer_isoform.py {input.bed} {input.fasta} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        
        set +u; source activate TAMA
        python ~/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {output.bed} {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """