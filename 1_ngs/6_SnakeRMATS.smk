#!/usr/bin/env snakemake

groups = [
    'Ad_Br_Fe', 'Ad_Br_Ma', 'Ad_Ey_Fe', 'Ad_Ey_Ma',
    'Ad_Gi_Fe', 'Ad_Gi_Ma', 'Ad_Go_Fe', 'Ad_Go_Ma',
    'Ad_He_Fe', 'Ad_He_Ma', 'Ad_In_Fe', 'Ad_In_Ma',
    'Ad_Ki_Fe', 'Ad_Ki_Ma', 'Ad_Li_Fe', 'Ad_Li_Ma',
    'Ad_Mu_Fe', 'Ad_Mu_Ma', 'Ad_Pi_Fe', 'Ad_Pi_Ma',
    'Ad_Sp_Fe', 'Ad_Sp_Ma', 'Ad_St_Fe', 'Ad_St_Ma',
    'Ju_Br_Mi', 'Ju_Ey_Mi', 'Ju_Gi_Mi', 'Ju_Go_Mi',
    'Ju_He_Mi', 'Ju_In_Mi', 'Ju_Ki_Mi', 'Ju_Li_Mi',
    'Ju_Mu_Mi'
]

pairs = []
for i in range(len(groups)):
    for j in range(i + 1, len(groups)):
        pairs.append("%s_vs_%s" % (groups[i], groups[j]))


rule all:
    input:
        expand("results/rMATS/pairs/{pair}", pair=pairs)


rule rmats:
    input:
        path1 = "results/rMATS/prep/{group1}/paths.txt",
        path2 = "results/rMATS/prep/{group2}/paths.txt",
        gtf = "results/stringtie/taco/merged/assembly.sorted.gtf"
    output:
        out = directory("results/rMATS/pairs/{group1}_vs_{group2}"),
        tmp = directory("results/rMATS/pairs/{group1}_vs_{group2}.tmp")
    log:
        "results/rMATS/pairs/{group1}_vs_{group2}.log"
    threads:
        4
    shell:
        """
        set +u
        source activate rMATS
        rmats.py --b1 {input.path1} --b2 {input.path2} --gtf {input.gtf} -t paired --readLength 100 --variable-read-length \
            --nthread {threads} --od {output.out} --tmp  {output.tmp} &> {log}
        conda deactivate
        """
