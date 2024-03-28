#!/bin/sh

wget https://ftp.ensembl.org/pub/release-103/gtf/seriola_dumerili/Seriola_dumerili.Sdu_1.0.103.gtf.gz
./convert_id.py Seriola_dumerili.Sdu_1.0.103.gtf.gz ../ncbi_Sdu_1.0/GCF_002260705.1_Sdu_1.0_genomic.ncbi_id_to_ddbj_id.txt > Seriola_dumerili.Sdu_1.0.103.converted.gtf
cat Seriola_dumerili.Sdu_1.0.103.converted.gtf | grep -v '#' | grep -v NC_016870.1 | awk '$3=="gene"||$3=="transcript"||$3=="exon"||$3=="CDS"' > Seriola_dumerili.Sdu_1.0.103.converted.clean.gtf
sort -k1,1 -k4,4n Seriola_dumerili.Sdu_1.0.103.converted.clean.gtf | bgzip -c > Seriola_dumerili.Sdu_1.0.103.converted.clean.sorted.gtf.gz
tabix -p gff Seriola_dumerili.Sdu_1.0.103.converted.clean.sorted.gtf.gz

