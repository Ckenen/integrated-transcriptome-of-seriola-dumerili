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