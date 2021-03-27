#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO


def main():
    path, outdir = sys.argv[1:]

    m = 0  # file number
    n = 0  # record number of one file
    mark = None
    fw = None  # file handle
    with gzip.open(path, "rt") as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            v1, v2, v3 = record.name.split("/")
            if mark is None or (v2 != mark and n >= 1000000):
                if fw:
                    fw.close()
                fw = open(outdir + "/%s.fasta" % m, "w+")
                m += 1
                n = 0
            mark = v2
            SeqIO.write(record, fw, "fasta")
            n += 1
    if fw:
        fw.close()


if __name__ == "__main__":
    main()
