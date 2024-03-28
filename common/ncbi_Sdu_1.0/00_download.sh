#!/bin/sh

# Download FASTA

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/260/705/GCF_002260705.1_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fna.gz
gzip -dc GCF_002260705.1_Sdu_1.0_genomic.fna.gz > GCF_002260705.1_Sdu_1.0_genomic.fa
samtools faidx GCF_002260705.1_Sdu_1.0_genomic.fa
faSize -detailed GCF_002260705.1_Sdu_1.0_genomic.fa > GCF_002260705.1_Sdu_1.0_genomic.sizes

# Download GTF

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/260/705/GCF_002260705.1_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.gtf.gz
gzip -dc GCF_002260705.1_Sdu_1.0_genomic.gtf.gz | grep -v '#' | grep -v NC_016870.1 | awk '$3=="gene"||$3=="exon"||$3=="CDS"' | grep -v unknown_transcript_1 > GCF_002260705.1_Sdu_1.0_genomic.clean.gtf
sort -k1,1 -k4,4n GCF_002260705.1_Sdu_1.0_genomic.clean.gtf | bgzip -c > GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.gtf.gz
tabix -p gff GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.gtf.gz

# Others

./ncbi_to_ddbj.py GCF_002260705.1_Sdu_1.0_genomic.fna.gz > GCF_002260705.1_Sdu_1.0_genomic.ncbi_id_to_ddbj_id.txt
./gtf_to_bed.py GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.gtf.gz > GCF_002260705.1_Sdu_1.0_genomic.clean.sorted.bed