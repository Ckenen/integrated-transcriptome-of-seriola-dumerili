#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
import glob
names = ["ngs", "tgs"]
outdir = "results/function"
name_batch_list = []
for name in names:
    for path in glob.glob(outdir + "/interpro/faSplited/%s/*.fa" % name):
        name_batch_list.append("%s/%s" % (name, path.split("/")[-1][:-3]))
# print(name_batch_list)

rule all:
    input:
        expand(outdir + "/interpro/faSplited/{name}", name=names),
        expand(outdir + "/interpro/pfam/{name_batch}.tsv", name_batch=name_batch_list),
        expand(outdir + "/interpro/PANTHER/{name_batch}.tsv", name_batch=name_batch_list),
        # outdir + "/interpro/pfam.tsv",
        # outdir + "/interpro/PANTHER.tsv",


rule split_fasta:
    input:
        fa = outdir + "/ORFfinder/{name}.fa"
    output:
        out = directory(outdir + "/interpro/faSplited/{name}")
    shell:
        """
        mkdir {output}
        faSplit about {input.fa} 1000000 {output.out}/
        """

rule InterProScan_Pfam:
    input:
        fas = rules.split_fasta.output.out
    output:
        tsv = outdir + "/interpro/pfam/{name}/{batch}.tsv",
        gff = outdir + "/interpro/pfam/{name}/{batch}.gff3"
    log:
        log = outdir + "/interpro/pfam/{name}/{batch}.log"
    threads:
        4
    params:
        prefix = outdir + "/interpro/pfam/{name}/{batch}"
    shell:
        """
        set +u; source activate InterProScan
        ~/software/interproscan-5.52-86.0/interproscan.sh -i {input.fas}/{wildcards.batch}.fa -f tsv,gff3 -appl Pfam -goterms -pathways --cpu {threads} -b {params.prefix} -dp &> {log}
        conda deactivate
        """

# rule InterProScan_Pfam_merge:
#     input:
#         tsvs = expand(outdir + "/interpro/pfam/{name}.tsv", name=names)
#     output:
#         tsv = outdir + "/interpro/pfam.tsv"
#     shell:
#         """
#         cat {input.tsvs} > {output.tsv}
#         """

rule InterProScan_PANTHER:
    input:
        fas = rules.split_fasta.output.out
    output:
        tsv = outdir + "/interpro/PANTHER/{name}/{batch}.tsv",
        gff = outdir + "/interpro/PANTHER/{name}/{batch}.gff3"
    log:
        log = outdir + "/interpro/PANTHER/{name}/{batch}.log"
    threads:
        4
    params:
        prefix = outdir + "/interpro/PANTHER/{name}/{batch}"
    shell:
        """
        set +u; source activate InterProScan
        ~/software/interproscan-5.52-86.0/interproscan.sh -i {input.fas}/{wildcards.batch}.fa -f tsv,gff3 -appl PANTHER -goterms -pathways --cpu {threads} -b {params.prefix} -dp &> {log}
        conda deactivate
        """

# rule InterProScan_PANTHER_merge:
#     input:
#         tsvs = expand(outdir + "/interpro/PANTHER/{name}.tsv", name=names)
#     output:
#         tsv = outdir + "/interpro/PANTHER.tsv"
#     shell:
#         """
#         cat {input.tsvs} > {output.tsv}
#         """
