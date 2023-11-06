#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
samples = ngs_samples
groups = list(sorted(set([s[:-2] for s in samples])))
# print(groups)
pairs = []
for g1 in groups:
    for g2 in groups:
        if g1 == g2:
            continue
        pairs.append("%s/%s" % (g1, g2))
outdir = "results/expression/deseq2"

rule all:
    input:
        expand(outdir + "/prep/{pair}.tsv", pair=pairs),
        expand(outdir + "/stat/{pair}.tsv", pair=pairs),

rule prep:
    input:
        tsv = "results/expression/featureCounts/asm.feature_count.tsv"
    output:
        tsv = outdir + "/prep/{group1}/{group2}.tsv"
    run:
        import pandas as pd
        fpkms = pd.read_csv(input.tsv, sep="\t", index_col=0)
        cols1 = list(filter(lambda item: item.startswith(wildcards.group1), fpkms.columns))
        cols2 = list(filter(lambda item: item.startswith(wildcards.group2), fpkms.columns))
        fpkms[cols1 + cols2].to_csv(output.tsv, sep="\t")

rule deseq2:
    input:
        tsv = outdir + "/prep/{group1}/{group2}.tsv"
    output:
        tsv = outdir + "/stat/{group1}/{group2}.tsv"
    log:
        outdir + "/stat/{group1}/{group2}.log"
    shell:
        """
        set +u; source activate DESeq2
        Rscript ./scripts/deseq2.R {input.tsv} {output.tsv} &> {log}
        """
