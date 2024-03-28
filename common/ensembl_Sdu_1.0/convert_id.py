#!/usr/bin/env python
import sys
import gzip


def main():
    gtffile, mapperfile = sys.argv[1:]

    mapper = dict()
    for line in open(mapperfile):
        ncbi_id, ensembl_id = line.strip("\n").split("\t")
        mapper[ensembl_id] = ncbi_id

    with gzip.open(gtffile, "rt") as f:
        for line in f:
            if line.startswith("#"):
                print(line, end="")
            else:
                row = line.strip("\n").split("\t")
                row[0] = mapper[row[0]]
                print("\t".join(row))


if __name__ == "__main__":
    main()