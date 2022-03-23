#!/usr/bin/env python
import sys
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder
infile1, infile2, infile3, outfile = sys.argv[1:]
dat1 = pd.read_csv(infile1, sep="\t", index_col=0)
dat2 = pd.read_csv(infile2, sep="\t", index_col=0)
dat = pd.concat([dat1, dat2], axis=0)
with GtfFile(infile3) as f:
    transcripts = list(GtfTranscriptBuilder(f))
transcript_name_list = set([t.name for t in transcripts])
dat = dat[[index in transcript_name_list for index in dat.index]]
dat.to_csv(outfile, sep="\t")