#!/usr/bin/env python
import sys
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder

with GtfFile(sys.argv[1]) as f:
    transcripts = list(GtfTranscriptBuilder(f))
dat = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
id_set = set(dat[dat.max(axis=1) < 0.5].index)

records1 = []
records2 = []
for transcript in transcripts:
    if transcript.chrom == "NC_016870.1":
        continue
    if transcript.name in id_set:
        records = records1
    else:
        records = records2
    for items in transcript.records.values():
        for item in items:
            records.append(item)

with open(sys.argv[3], "w+") as fw:
    for record in sorted(records1):
        fw.write(record.format() + "\n")
with open(sys.argv[4], "w+") as fw:
    for record in sorted(records2):
        fw.write(record.format() + "\n")