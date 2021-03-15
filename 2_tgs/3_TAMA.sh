#!/bin/bash

# Adult_Female 
# Adult_Male 
# Juvenile

source activate py27

# Step 1: Collapse
# for name in Adult_Female Adult_Male Juvenile; do
#     python ./scripts/tama/tama_collapse.py \
# 		-b BAM -s results/isoseq/aligned/${name}.bam 
# 		-f data/genome/genome.fasta \
# 		-p results/tama/collapsed/${name} \
# 		-x no_cap &> results/tama/collapsed/${name}.log
# done

# Step 2: Merge
# mkdir -p results/tama/merged
# cat > results/tama/merged/filelist.txt << EOF
# results/tama/collapsed/Adult_Female.bed	no_cap	1,1,1	AdFe
# results/tama/collapsed/Adult_Male.bed	no_cap	1,1,1	AdMa
# results/tama/collapsed/Juvenile.bed	no_cap	1,1,1	JuMi
# EOF
# python ./scripts/tama/tama_merge.py \
# 	-f results/tama/merged/filelist.txt \
# 	-p results/tama/merged/tama_merged > results/tama/merged/tama_merged.log


# Step 3: Read support
# mkdir -p results/tama/read_support
# cat > results/tama/read_support/filelist.txt << EOF
# AdFe	results/tama/collapsed/Adult_Female_trans_read.bed	trans_read
# AdMa	results/tama/collapsed/Adult_Male_trans_read.bed	trans_read
# JuMi	results/tama/collapsed/Juvenile_trans_read.bed	trans_read
# EOF
# python ./scripts/tama/tama_go/read_support/tama_read_support_levels.py \
# 	-f results/tama/read_support/filelist.txt \
# 	-m results/tama/merged/tama_merged_merge.txt \
# 	-o results/tama/read_support/all > results/tama/read_support/all.log

# Step 4: Remove PolyA
# mkdir -p results/tama/remove_polya
# cat > results/tama/remove_polya/filelist.txt << EOF
# AdFe	results/tama/collapsed/Adult_Female_polya.txt
# AdMa	results/tama/collapsed/Adult_Male_polya.txt
# JuMi	results/tama/collapsed/Juvenile_polya.txt
# EOF
# python ./scripts/tama/tama_go/filter_transcript_models/tama_remove_polya_models_levels.py \
# 	-b results/tama/merged/tama_merged.bed \
# 	-f results/tama/remove_polya/filelist.txt \
# 	-r results/tama/read_support/all_read_support.txt \
# 	-o results/tama/remove_polya/remove_polya > results/tama/remove_polya/remove_polya.log

# End

conda deactivate
