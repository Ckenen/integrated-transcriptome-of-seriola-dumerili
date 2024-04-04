#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import FamFileRandom

infile, outfile = sys.argv[1:]
with FamFileRandom(infile) as f, open(outfile, "w+") as fw:
    for i, fragment in enumerate(f):
        fragment.strand = fragment.mate2.strand
        fw.write(fragment.format("BED") + "\n")
