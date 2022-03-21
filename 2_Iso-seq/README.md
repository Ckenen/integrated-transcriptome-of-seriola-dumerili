# Subreads结构（15%的错误率，构建CCS之后检查出来的）

AAGCAGTGGTATCAACGCAGAGTACATGGGG[mRNA]AAAAAAAAAA[8xBARCODE]GTACTCTGCGTTGATACCACTGCTT
AAGCAGTGGTATCAACGCAGAGTAC[8*BARCODE]TTTTTTTTTT[mRNA]CCCCATGTACTCTGCGTTGATACCACTGCTT

# 使用lima程序执行完之后，reads的结构中不包含有adapter，但包含有barcode，在polyA的尾部，8个碱基

[mRNA]AAAAAAAAAA[8xBARCODE]
[8*BARCODE]TTTTTTTTTT[mRNA]

# 

export PATH="${HOME}/software/smrtlink/smrtcmds/bin:${PATH}"

results/isoseq/collapsed/Adult_Female.gtf是isoseq3 collapse输出的文件，直接将gff后缀改为gtf后缀，不能排序，否则会报错。
~/software/SQANTI3/sqanti3_qc.py --gtf results/isoseq/collapsed/Adult_Female.gtf  data/genome/annotation.gtf data/genome/genome.fasta --fl_count results/isoseq/collapsed/Adult_Female.abundance.txt -t 30 -n 10



minimap2 -ax splice -t 30 -uf --secondary=no -C5 hg38.fa hq_isoforms.fastq > hq_isoforms.fastq.sam
sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
collapse_isoforms_by_sam.py --input hq_isoforms.fastq --fq -s hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o test



collapse_isoforms_by_sam.py --input results/isoseq/polished/Adult_Female.hq.fastq --fq -s results/isoseq/aligned/Adult_Female.bam --dun-merge-5-shorter -o test
使用pbmm2进行比对的结果出问题，报错了。


minimap2 -ax splice -t 40 --secondary=no -C5 data/genome/genome.fasta results/isoseq/polished/Adult_Female.hq.fastq > hq_isoforms.fastq.sam
sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
collapse_isoforms_by_sam.py --input results/isoseq/polished/Adult_Female.hq.fastq --fq -s hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o cupcake


filter_away_subset.py test.collapsed


fusion_finder.py --input results/cupcake/mapping/Adult_Female.hq.fastq --fq -s results/cupcake/mapping/Adult_Female.sorted.bam.sam --cluster_report results/isoseq/polished/Adult_Female.cluster_report.csv -o output.fusion --min_locus_coverage_bp 500 -d 1000000


~/software/SQANTI3/sqanti3_qc.py --gtf results/cupcake/fusion/Adult_Female/out.gff data/genome/annotation.gtf data/genome/genome.fasta --is_fusion -o results/cupcake/sqanti3/A


fusion_collate_info.py output.fusion output.fusion_classification.txt \
      refAnnotation_output.fusion.genePred \
      --genome hg38.fa

(py27) chenzonggui@luolab:~/4_cwork/2_tgs/results/tama$ python ~/software/tama/tama_go/filter_transcript_models/tama_remove_polya_models_levels.py -b merged/tama_merged.bed -f filelist.txt -o output -r R


# BLASTP 2.9.0+
# Query: lcl|ORF1_TU2:393:509 unnamed protein product
# Database: ../ncbi/swissprot/swissprot
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 7 hits found

cat annotation.gtf | awk '$3=="transcript"' | head | sed -e 's/^.*gene_id "\(\S*\)";.*/\1/g'

cat annotation.gtf | awk '$3=="transcript"' | sed -e 's/^.*gene_id "\(\S*\)";.*/\1/g' | sort | uniq -c | awk '{print $1}' | sort | uniq -c | sort -k2,2n | awk '{print $2"\t"$1}'