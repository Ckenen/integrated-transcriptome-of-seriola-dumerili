#!/usr/bin/env python
import sys
import pandas as pd

infile1, infile2, outfile = sys.argv[1:]
diamond = pd.read_csv(infile1, sep="\t")
subject_set = set(diamond["Subject"])
header = None
rows = []
with open(infile2) as f:
    for i, line in enumerate(f):
        line = line.strip("\n")
        row = line.split("\t")
        if i == 0:
            header = row
        else:
            if row[1] in subject_set:
                rows.append(row)
dat = pd.DataFrame(rows, columns=header)
dat = dat[dat.columns[1:]]
dat = diamond.merge(dat, left_on="Subject", right_on="Name")
dat = dat.sort_values(by="ID")
dat.to_csv(outfile, sep="\t", index=False)