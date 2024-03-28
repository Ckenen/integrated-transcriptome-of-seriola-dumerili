#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder


def main():
    with GtfFile(sys.argv[1]) as f:
        records = [x for x in f]
        transcripts = list(GtfTranscriptBuilder(records))
        for t in sorted(transcripts):
            print(t.format("BED"))


if __name__ == "__main__":
    main()