#!/usr/bin/env python
import sys
import os
import re
import pandas as pd

array = []
for infile in sys.argv[1:-1]:
    sample = os.path.basename(infile)
    names = []
    values = []
    with open(infile + "/transcripts.gtf") as f:
        for line in f:
            ret = re.search("FPKM \"[0-9\.]+\";", line)
            if ret:
                fpkm = float(line[ret.start() + 6:ret.end() - 2])
                ret = re.search("transcript_id \"[\S]+\";", line)
                tid = line[ret.start() + 15:ret.end() - 2]
                names.append(tid)
                values.append(fpkm)
    series = pd.Series(values)
    series.index = names
    series.name = sample
    array.append(series)
dat = pd.concat(array, axis=1)
dat.index.name = "TranscriptID"
dat.to_csv(sys.argv[-1], sep="\t")