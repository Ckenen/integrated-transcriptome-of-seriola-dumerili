#!/usr/bin/env python
import sys
import pandas as pd

infile, group1, group2, outfile = sys.argv[1:]
counts = pd.read_csv(infile, sep="\t", index_col=0)
cs1 = list(filter(lambda item: item.startswith(group1), counts.columns))
cs2 = list(filter(lambda item: item.startswith(group2), counts.columns))
assert len(cs1) == 2 and len(cs2) == 2
counts[cs1 + cs2].to_csv(outfile, sep="\t")