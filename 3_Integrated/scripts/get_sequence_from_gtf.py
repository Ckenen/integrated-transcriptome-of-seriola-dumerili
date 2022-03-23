#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, FastaFile

with GtfFile(sys.argv[1]) as f, FastaFile(sys.argv[2]) as fasta:
    for transcript in GtfTranscriptBuilder(f):
        seq = fasta.fetch(obj=transcript).upper()
        print(">%s" % transcript.name)
        print(seq)

