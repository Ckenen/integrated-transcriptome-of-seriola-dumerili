#!/bin/sh

# FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/260/705/GCF_002260705.1_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.fna.gz
gzip -c -d GCF_002260705.1_Sdu_1.0_genomic.fna.gz > GCF_002260705.1_Sdu_1.0_genomic.fasta
samtools faidx GCF_002260705.1_Sdu_1.0_genomic.fasta
faSize -detailed GCF_002260705.1_Sdu_1.0_genomic.fasta > GCF_002260705.1_Sdu_1.0_genomic.sizes


# GFF
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/260/705/GCF_002260705.1_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.gff.gz
gzip -d -c GCF_002260705.1_Sdu_1.0_genomic.gff.gz > GCF_002260705.1_Sdu_1.0_genomic.gff
bedtools sort -i GCF_002260705.1_Sdu_1.0_genomic.gff > GCF_002260705.1_Sdu_1.0_genomic.sorted.gff
bgzip -c GCF_002260705.1_Sdu_1.0_genomic.sorted.gff > GCF_002260705.1_Sdu_1.0_genomic.sorted.gff.gz
tabix -p gff GCF_002260705.1_Sdu_1.0_genomic.sorted.gff.gz

# GTF
./gffToGTF.py GCF_002260705.1_Sdu_1.0_genomic.sorted.gff GCF_002260705.1_Sdu_1.0_genomic.temp.gtf
bedtools sort -i GCF_002260705.1_Sdu_1.0_genomic.temp.gtf > GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf
rm GCF_002260705.1_Sdu_1.0_genomic.temp.gtf
bgzip -c GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf > GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf.gz
tabix -p gff GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf.gz

# BED
./gtfToBed.py GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf > GCF_002260705.1_Sdu_1.0_genomic.bed
bedtools sort -i GCF_002260705.1_Sdu_1.0_genomic.bed > GCF_002260705.1_Sdu_1.0_genomic.sorted.bed
bgzip -c GCF_002260705.1_Sdu_1.0_genomic.sorted.bed > GCF_002260705.1_Sdu_1.0_genomic.sorted.bed.gz
tabix -p bed GCF_002260705.1_Sdu_1.0_genomic.sorted.bed.gz

# TSV
./gtfToTSV.py GCF_002260705.1_Sdu_1.0_genomic.sorted.gtf GCF_002260705.1_Sdu_1.0_genomic.sorted.bed GCF_002260705.1_Sdu_1.0_genomic.tsv

# Product
./gffInfo.py GCF_002260705.1_Sdu_1.0_genomic.sorted.gff GCF_002260705.1_Sdu_1.0_genomic.info.tsv

# lastz_32 genome.fasta[multiple] human.sox9.fasta --step=1 --ambiguous=iupac --nochain --strand=both --format=axt > map_human_sox9.axt
# axtChain -faT -faQ -linearGap=loose map_human_sox9.axt genome.fasta human.sox9.fasta map_human_sox9.chain

# Self-alignment
lastz GCF_002260705.1_Sdu_1.0_genomic.fasta[multiple] --self --nomirror --format=axt | gzip -c > lastz.self.align.axt.gz


# CPAT
bedtools getfasta -split -s -name+ -fi GCF_002260705.1_Sdu_1.0_genomic.fasta -bed cpat.protein_coding.cds.bed > cpat.protein_coding.cds.fasta
bedtools getfasta -split -s -name+ -fi GCF_002260705.1_Sdu_1.0_genomic.fasta -bed cpat.lncRNA.bed > cpat.lncRNA.fasta
make_hexamer_tab.py -c cpat.protein_coding.cds.fasta -n cpat.lncRNA.fasta > GTS_Hexamer.tsv
make_logitModel.py  -x GTS_Hexamer.tsv -c cpat.protein_coding.cds.fasta -n cpat.lncRNA.fasta -o GTS
