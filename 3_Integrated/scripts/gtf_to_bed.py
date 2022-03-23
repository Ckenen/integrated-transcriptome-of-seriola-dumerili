#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder

with GtfFile(sys.argv[1]) as f:
    for transcript in sorted(GtfTranscriptBuilder(f)):
        print(transcript.format())