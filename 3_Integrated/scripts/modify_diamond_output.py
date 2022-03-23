#!/usr/bin/env python
import sys
import pandas as pd

# 输出格式：
# qseqid query (e.g., unknown gene) sequence id
# sseqid subject (e.g., reference genome) sequence id
# pident percentage of identical matches
# length alignment length (sequence overlap)
# mismatch number of mismatches
# gapopen number of gap openings
# qstart start of alignment in query
# qend end of alignment in query
# sstart start of alignment in subject
# send end of alignment in subject
# evalue expect value
# bitscore bit score

infile, outfile = sys.argv[1:]
dat = pd.read_csv(infile, sep="\t", header=None)
dat.columns = ["ID", "Subject", "PercIdent", "ORF_Length", "Mismatch", "GapOpen", 
            "Query_Start", "Query_End", "Subject_Start", "Subject_End", 
            "Evalue", "BitScore"]
dat["Query_Start"] = dat["Query_Start"] - 1
dat["Subject_Start"] = dat["Subject_Start"] - 1
dat.to_csv(outfile, sep="\t", index=False)
