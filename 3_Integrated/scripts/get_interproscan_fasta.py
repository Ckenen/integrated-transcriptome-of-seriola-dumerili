#!/usr/bin/env python
import sys
import pandas as pd
from Bio.Seq import Seq
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, FastaFile

# cpat output
cpat = sys.argv[1]
# diamond output
diamond = sys.argv[2]
orffinder = sys.argv[3]
gtf = sys.argv[4]
fasta = sys.argv[5]

with GtfFile(gtf) as f:
    transcripts = list(GtfTranscriptBuilder(f))
with FastaFile(fasta) as fasta:
    for transcript in transcripts:
        transcript.seq = fasta.fetch(obj=transcript)
transcript_dict = {t.name: t for t in transcripts}

dat1 = pd.read_csv(cpat, sep="\t", index_col=0)
dat2 = pd.read_csv(diamond, sep="\t", index_col=0)
dat3 = pd.read_csv(orffinder, sep="\t", index_col=0)
dat1 = dat1[dat1["Coding_Prob"] >= 0.5]
dat2 = dat2[(dat2["Evalue"] < 1e-6) & (dat2["PercIdent"] > 80)]
flags = []
stop_codon_list = ["TAA", "TGA", "TAG"]
for start_codon, stop_codon in dat3[["Start_Codon", "Stop_Codon"]].values:
    flags.append(start_codon == "ATG" and stop_codon in stop_codon_list)
set1 = set(dat1.index)
set2 = set(dat3[flags].index) & set(dat2.index)
set3 = set1 | set2
for v in list(set3):
    name, start, end = v.split(":")
    start = int(start)
    end = int(end)
    transcript = transcript_dict[name]
    cds = transcript.seq[start:end].upper()
    s = Seq(cds[:-3])
    s = str(s.translate())
    print(">%s" % v)
    print(s)