#!/usr/bin/env python
import sys
import gzip
import re
from Bio import SeqIO


def main():
    with gzip.open(sys.argv[1], "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == "NC_016870.1":
                name = "MT"
            else:
                ret = re.search("BDQW[0-9]{8}\.1", record.description)
                assert ret
                name = record.description[ret.start():ret.end()]
            print(record.id, name, sep="\t")


if __name__ == "__main__":
    main()