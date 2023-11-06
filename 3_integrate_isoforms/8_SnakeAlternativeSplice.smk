#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
samples = ngs_samples
groups = list(set([s[:-2] for s in samples]))
groups.sort()
# groups = groups[:3]

# pairs1 = []
# for i in range(len(samples)):
#     for j in range(i + 1, len(samples)):
#         pairs1.append("%s_vs_%s" % (samples[i], samples[j]))
# print(len(pairs1))

pairs2 = []
for i in range(len(groups)):
    for j in range(i + 1, len(groups)):
        pairs2.append("%s_vs_%s" % (groups[i], groups[j]))
print(len(pairs2))

indir = "../1_analysis_of_rna_seq/results/mapping/rmdup"
outdir = "results/rmats"

rule all:
    input:
        outdir + "/ref.gtf",
        expand(outdir + "/prep/{sample}.tmp", sample=samples),
        # expand(outdir + "/pairs1/{pair}", pair=pairs1),
        expand(outdir + "/pairs2/{pair}", pair=pairs2)

rule make_ref_gtf:
    input:
        gtf = GTFS["asm"]
    output:
        gtf = outdir + "/ref.gtf"
    shell:
        """
        gzip -d -c {input.gtf} | awk '$3!="gene"' > {output.gtf}
        """

rule rmat_prep: # v4.2.0
    input:
        bam = indir + "/{sample}.bam",
        gtf = rules.make_ref_gtf.output.gtf
    output:
        txt = outdir + "/prep/{sample}.txt",
        out = directory(outdir + "/prep/{sample}"),
        tmp = directory(outdir + "/prep/{sample}.tmp")
    log:
        outdir + "/prep/{sample}.log"
    threads:
        4
    shell:
        """
        echo {input.bam} > {output.txt}
        set +u; source activate rMATS
        rmats.py --b1 {output.txt} --gtf {input.gtf} -t paired --readLength 100 --variable-read-length \
            --libType fr-firststrand --nthread {threads} --od {output.out} --tmp {output.tmp} --task prep &> {log}
        conda deactivate
        """

# rule rmats_1:
#     input:
#         bam1 = indir + "/{sample1}.bam",
#         bam2 = indir + "/{sample2}.bam",
#         prep1 = outdir + "/prep/{sample1}.tmp",
#         prep2 = outdir + "/prep/{sample2}.tmp",
#         gtf = outdir + "/ref.gtf"
#     output:
#         txt1 = temp(outdir + "/pairs1/{sample1}_vs_{sample2}.batch1.txt"),
#         txt2 = temp(outdir + "/pairs1/{sample1}_vs_{sample2}.batch2.txt"),
#         tmp = temp(directory(outdir + "/pairs1/{sample1}_vs_{sample2}.tmp")),
#         out = directory(outdir + "/pairs1/{sample1}_vs_{sample2}")
#     log:
#         log = outdir + "/pairs1/{sample1}_vs_{sample2}.log"
#     threads:
#         4
#     shell:
#         """
#         echo {input.bam1} > {output.txt1}
#         echo {input.bam2} > {output.txt2}
#         set +u; source activate rMATS420
#         mkdir {output.tmp}
#         python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.sample1}_ {output.tmp} {input.prep1}/*.rmats
#         python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.sample2}_ {output.tmp} {input.prep2}/*.rmats
#         rmats.py --b1 {output.txt1} --b2 {output.txt2} --gtf {input.gtf} -t paired --readLength 100 --variable-read-length \
#             --libType fr-firststrand --nthread {threads} --od {output.out} --tmp  {output.tmp} --task post &> {log}
#         conda deactivate
#         """

def get_rmats_2_input_list(wildcards):
    array1 = list(filter(lambda item: item.startswith(wildcards.group1), samples))
    array2 = list(filter(lambda item: item.startswith(wildcards.group2), samples))
    assert len(array1) == 2
    assert len(array2) == 2
    bam1 = indir + "/%s.bam" % array1[0]
    bam2 = indir + "/%s.bam" % array1[1]
    bam3 = indir + "/%s.bam" % array2[0]
    bam4 = indir + "/%s.bam" % array2[1]
    prep1 = outdir + "/prep/%s.tmp" % array1[0]
    prep2 = outdir + "/prep/%s.tmp" % array1[1]
    prep3 = outdir + "/prep/%s.tmp" % array2[0]
    prep4 = outdir + "/prep/%s.tmp" % array2[1]
    return [bam1, bam2, bam3, bam4, prep1, prep2, prep3, prep4]

rule rmats_2:
    input:
        paths = lambda wildcards: get_rmats_2_input_list(wildcards),
        gtf = outdir + "/ref.gtf"
    output:
        txt1 = outdir + "/pairs2/{group1}_vs_{group2}.batch1.txt",
        txt2 = outdir + "/pairs2/{group1}_vs_{group2}.batch2.txt",
        tmp = directory(outdir + "/pairs2/{group1}_vs_{group2}.tmp"),
        out = directory(outdir + "/pairs2/{group1}_vs_{group2}")
    log:
        log = outdir + "/pairs2/{group1}_vs_{group2}.log"
    threads:
        4
    shell:
        """
        echo {input.paths[0]},{input.paths[1]} > {output.txt1}
        echo {input.paths[2]},{input.paths[3]} > {output.txt2}
        set +u; source activate rMATS
        mkdir {output.tmp}
        python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.group1}_1_ {output.tmp} {input.paths[4]}/*.rmats
        python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.group1}_2_ {output.tmp} {input.paths[5]}/*.rmats
        python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.group2}_1_ {output.tmp} {input.paths[6]}/*.rmats
        python ~/software/rmats_turbo_v4_1_1/cp_with_prefix.py {wildcards.group2}_2_ {output.tmp} {input.paths[7]}/*.rmats
        rmats.py --b1 {output.txt1} --b2 {output.txt2} --gtf {input.gtf} -t paired --readLength 100 --variable-read-length \
            --libType fr-firststrand --nthread {threads} --od {output.out} --tmp  {output.tmp} --task post &> {log}
        conda deactivate
        """
