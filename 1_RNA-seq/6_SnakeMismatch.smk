#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]

indir = "results/mapping/rmdup"
outdir = "results/mismatch"

rule all:
    input:
        # sample
        expand(outdir + "/single/reference_base/{sample}.tsv", sample=samples),
        expand(outdir + "/single/base_coverage/{sample}.{s}.bw", sample=samples, s=["+", "-"]),
        expand(outdir + "/single/mismatch_events/{sample}.tsv.gz", sample=samples),
        # remove SNPs
        expand(outdir + "/noSNPs/mismatch_events/{sample}.tsv.gz", sample=samples),
        expand(outdir + "/noSNPs/mismatch_ratio/{sample}", sample=samples),
        # candidate
        expand(outdir + "/candidate/mismatch_sites/{sample}.tsv.gz", sample=samples),
        outdir + "/candidate/mismatch_sites.union.tsv.gz",
        expand(outdir + "/candidate/base_coverage/{sample}.tsv.gz", sample=samples),
        expand(outdir + "/candidate/mismatch_events/{sample}.tsv.gz", sample=samples),
        # Merge all info
        outdir + "/merged/all.tsv.gz",

# Sample

rule get_reference_base:
    input:
        fam = indir + "/{sample}.fam",
        fasta = config["genome_fasta"]
    output:
        tsv = outdir + "/single/reference_base/{sample}.tsv"
    log:
        outdir + "/single/reference_base/{sample}.log"
    shell:
        """
        ./scripts/mismatch/get_reference_base.py {input.fam} {input.fasta} {output.tsv} &> {log}
        """

rule get_base_coverage:
    input:
        fam = indir + "/{sample}.fam",
        txt = config["genome_size"]
    output:
        bw1 = outdir + "/single/base_coverage/{sample}.+.bw",
        bw2 = outdir + "/single/base_coverage/{sample}.-.bw"
    log:
        outdir + "/single/base_coverage/{sample}.log"
    params:
        prefix = outdir + "/single/base_coverage/{sample}"
    shell:
        """(
        ./scripts/mismatch/get_base_coverage.py {input.fam} {params.prefix}.+.bedGraph {params.prefix}.-.bedGraph
        bedGraphToBigWig {params.prefix}.+.bedGraph {input.txt} {params.prefix}.+.bw
        bedGraphToBigWig {params.prefix}.-.bedGraph {input.txt} {params.prefix}.-.bw
        rm {params.prefix}.+.bedGraph {params.prefix}.-.bedGraph ) &> {log}
        """

rule get_mismatch_events:
    input:
        fam = indir + "/{sample}.fam"
    output:
        tmp = temp(outdir + "/single/mismatch_events/{sample}.tmp.tsv"),
        tsv = outdir + "/single/mismatch_events/{sample}.tsv.gz"
    log:
        outdir + "/single/mismatch_events/{sample}.log"
    shell:
        """(
        ./scripts/mismatch/get_mismatch_events.py {input.fam} {output.tmp} 
        sort -k1,1 -k2,2n -k3,3 -T ./ {output.tmp} | gzip -c > {output.tsv} ) &> {log}
        """

# Remove SNPs

rule remove_mismatch_events:
    input:
        tsv = rules.get_mismatch_events.output.tsv,
        bed = "../8_wgs/results/mismatch/confidence/mismatch_sites.conf.bed.gz"
    output:
        tsv = outdir + "/noSNPs/mismatch_events/{sample}.tsv.gz"
    shell:
        """
        ./scripts/mismatch/filter_mismatch_events.py {input.tsv} {input.bed} | gzip -c > {output.tsv}
        """

rule calculate_mismatch_ratio:
    input:
        tsv1 = rules.get_reference_base.output.tsv,
        tsv2 = rules.remove_mismatch_events.output.tsv
    output:
        directory(outdir + "/noSNPs/mismatch_ratio/{sample}")
    params:
        prefix = outdir + "/noSNPs/mismatch_ratio/{sample}/{sample}"
    shell:
        """
        mkdir -p {output}
        ./scripts/mismatch/calculate_mismatch_ratio.py {input.tsv1} {input.tsv2} {params.prefix}
        """

# Candidate

rule get_candidate_mismatch_sites:
    input:
        tsv = rules.get_mismatch_events.output.tsv
    output:
        tsv = outdir + "/candidate/mismatch_sites/{sample}.tsv.gz"
    shell:
        """
        ./scripts/mismatch/get_candidate_mismatch_sites.py {input.tsv} {output.tsv}
        """

rule merge_candidate_mismatch_sites:
    input:
        expand(rules.get_candidate_mismatch_sites.output.tsv, sample=samples)
    output:
        tsv = outdir + "/candidate/mismatch_sites.union.tsv.gz"
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2n -k3,3 -u | gzip -c > {output.tsv}
        """

rule get_candidate_base_coverage:
    input:
        tsv = rules.merge_candidate_mismatch_sites.output.tsv,
        bw1 = rules.get_base_coverage.output.bw1,
        bw2 = rules.get_base_coverage.output.bw2
    output:
        tsv = outdir + "/candidate/base_coverage/{sample}.tsv.gz"
    shell:
        """
        ./scripts/mismatch/get_candidate_base_coverage.py {input.tsv} {input.bw1} {input.bw2} {output.tsv}
        """

rule get_candidate_mismatch_events:
    input:
        tsv1 = rules.merge_candidate_mismatch_sites.output.tsv,
        tsv2 = rules.get_mismatch_events.output.tsv
    output:
        tsv = outdir + "/candidate/mismatch_events/{sample}.tsv.gz"
    shell:
        """
        ./scripts/mismatch/get_candidate_mismatch_events.py {input.tsv1} {input.tsv2} {output.tsv}
        """

# Merge

rule merge_all_info:
    input:
        tsv = rules.merge_candidate_mismatch_sites.output.tsv,
        covs = [outdir + "/candidate/base_coverage/%s.tsv.gz" % s for s in samples],
        events = [outdir + "/candidate/mismatch_events/%s.tsv.gz" % s for s in samples]
    output:
        tsv = outdir + "/merged/all.tsv.gz"
    shell:
        """
        ./scripts/mismatch/merge_all_info.py {input.tsv} {input.covs} {input.events} {output.tsv}
        """