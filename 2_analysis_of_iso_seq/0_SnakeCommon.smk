#!/bin/sh/env runsnakemake
SAMPLES = ["Ad_Ma", "Ad_Fe", "Ju_Mi"]
# MIN_PASSES_LIST = [1, 2, 4][2:]
MIN_PASSES = 4
THREADS = 20
GENOME_FASTA = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fa"
TAMA_ROOT = "/home/chenzonggui/software/tama"
