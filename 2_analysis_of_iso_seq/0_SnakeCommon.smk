#!/bin/sh/env snakemake
samples = [
    "Ad_Ma",
    "Ad_Fe",
    "Ju_Mi"
]

min_passes_list = [1, 2, 4]
threads = 20

GENOME_FASTA = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fa"
