#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO

with gzip.open(sys.argv[1], "rt") as f:
    for record in SeqIO.parse(f, "genbank"):
        print(">" + record.name)
        seq = str(record.seq)
        for i in range(0, len(seq), 60):
            j = min(i + 60, len(seq))
            assert i < j
            print(seq[i:j])
