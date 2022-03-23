#!/usr/bin/env python
import sys
import pandas as pd
from Bio import SeqIO

# diamond modified output
infile1 = sys.argv[1]
# orffinder orf
infile2 = sys.argv[2]
# output
outfile = sys.argv[3]

dat = pd.read_csv(infile1, sep="\t", index_col=0)
array = set(dat["qseqid"])
with open(outfile, "w+") as fw:
    for record in SeqIO.parse(infile2, "fasta"):
        if record.id in array:
            SeqIO.write(record, fw, "fasta")