#!/bin/sh

# bowtie2 -p 80 --local -x bowtie2_index/ref -1 Cleandata/GTS-GN-11-2A/GTS-GN-11-2A_1.fq.gz -2 Cleandata/GTS-GN-11-2A/GTS-GN-11-2A_2.fq.gz -S GTS-GN-11-2A.sam > GTS-GN-11-2A.log
# bowtie2 -p 80 --local -x bowtie2_index/ref -1 Cleandata/GTS-GN-13-2A/GTS-GN-13-2A_1.fq.gz -2 Cleandata/GTS-GN-13-2A/GTS-GN-13-2A_2.fq.gz -S GTS-GN-13-2A.sam > GTS-GN-13-2A.log
# bowtie2 -p 80 --local -x bowtie2_index/ref -1 Cleandata/GTS-GN-14-2A/GTS-GN-14-2A_1.fq.gz -2 Cleandata/GTS-GN-14-2A/GTS-GN-14-2A_2.fq.gz -S GTS-GN-14-2A.sam > GTS-GN-14-2A.log

for smk in 2_SnakeMapping.smk  3_SnakeExpression.smk  3_SnakeExpressionUniq.smk; do
    snakemake -s $smk -j 80 --ri -k
done

cpat.py -r data/genome/genome.fasta -g results/taco/stringtie/assembly.bed -d ../ncbi/Sdu_1.0/GTS.logit.RData -x ../ncbi/Sdu_1.0/GTS_Hexamer.tsv  -o cpat_out



ORFfinder -in input.fasta -strand plus -out orf_finder_out.orf.fasta

head -n 1000 merged/assembly.bed > orf_finder.test.bed
bedtools getfasta -s -split -name -fi ../../../data/genome/genome.fasta -bed orf_finder.test.bed | sed 's/([0-9.+-]*)//g' | sed 's/G[0-9]*|//g' > orf_finder.test.fasta
ORFfinder -in orf_finder.test.fasta -strand plus -out orf_finder_out.fasta

blastp -num_threads 8 