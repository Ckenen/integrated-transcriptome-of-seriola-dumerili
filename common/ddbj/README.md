# Paper

Whole Genome Sequencing of Greater Amberjack (Seriola dumerili) for SNP Identification on Aligned Scaffolds and Genome Structural Variation Analysis Using Parallel Resequencing

# Download from DDBJ database

BDQW.gz is DNA sequences (34655 scaffolds and contigs)

IACO.gz is transcript sequences (45109 transcripts)

# 格式转换

./genbank_to_fasta.py BDQW.gz > serDum.ddbj.genome.fasta
./genbank_to_fasta.py IACO.gz > serDum.ddbj.transcripts.fasta