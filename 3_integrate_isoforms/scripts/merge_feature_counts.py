#!/usr/bin/env python
import sys
import os
import pandas as pd

array = []
for infile in sys.argv[1:-1]:
    basename = os.path.basename(infile)
    assert basename.endswith(".txt")
    name = basename[:-4]
    dat = pd.read_csv(infile, sep="\t", index_col=0, comment="#")
    series = dat[dat.columns[-1]]
    series.name = name
    array.append(series)
dat = pd.concat(array, axis=1)
dat.to_csv(sys.argv[-1], sep="\t")
