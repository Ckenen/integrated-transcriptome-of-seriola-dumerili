#!/usr/bin/env snakemake
import numpy as np
import pandas as pd
dat = pd.read_excel("RNAseq.xls")

samples = list(dat["Sample"])
print("Samples: %d" % len(samples))

tmp = dat[~np.isnan(dat["Replicate"])] # final samples
final_samples = list(tmp["Sample"])
print("Final samples: %d" % len(final_samples))

GENOME_FASTA = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fa"
GENOME_SIZES = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.sizes"
ANNOTATION_GTF = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.clean.gtf"
ANNOTATION_GTF_GZ = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.gtf.gz"
ANNOTATION_BED = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.bed"
ANNOTATION_BED_GZ = "../common/ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.bed.gz"
