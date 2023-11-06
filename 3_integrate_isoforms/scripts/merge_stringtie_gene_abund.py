#!/usr/bin/env python
import sys
import os
import pandas as pd

array = []
for infile in sys.argv[1:-1]:
    dat = pd.read_csv(infile + "/gene_abund.tab", sep="\t", index_col=0)
    dat = dat[~dat.index.duplicated()]
    series = dat["FPKM"]
    series.name = os.path.basename(infile)
    array.append(series)
dat = pd.concat(array, axis=1)
dat.index.name = "GeneID"
dat.to_csv(sys.argv[-1], sep="\t")