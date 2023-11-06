#!/usr/bin/env snakemake
import numpy as np
import pandas as pd
dat = pd.read_excel("../1_analysis_of_rna_seq/RNAseq.xls")
ngs_samples = list(dat[~np.isnan(dat["Replicate"])]["Sample"])
print("Samples: %d" % len(ngs_samples))

GENOME_FASTA = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fa"

GTFS = {
    "ncbi": "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.gtf.gz",
    "ensembl": "../common/ensembl_Sdu_1.0/Seriola_dumerili.Sdu_1.0.103.converted.clean.sorted.gtf.gz",
    "ngs": "../1_analysis_of_rna_seq/results/assembly/stringtie/merged_all_samples.sorted.gtf.gz",
    "tgs": "../2_analysis_of_iso_seq/results/assembly/tama/filtered_internal_primer/all_samples.mp4.sorted.gtf.gz",
    "asm": "results/assembly/asm.final.sorted.gtf.gz",
}