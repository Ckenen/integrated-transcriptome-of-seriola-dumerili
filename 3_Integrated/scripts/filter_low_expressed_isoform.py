#!/usr/bin/env python
import sys
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder

with GtfFile(sys.argv[1]) as f:
    transcripts = list(GtfTranscriptBuilder(f))
dat = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
id_set = set(dat[dat.max(axis=1) >= 0.5].index)
records = []
for transcript in transcripts:
    if transcript.name not in id_set:
        continue
    for items in transcript.records.values():
        for item in items:
            records.append(item)
for record in sorted(records):
    print(record.format())