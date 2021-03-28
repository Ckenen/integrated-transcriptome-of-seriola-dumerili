#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile


def main():
    path1, path2 = sys.argv[1:]

    chrom_length_dict = dict()
    with open(path1) as f:
        for line in f:
            chrom, length = line.split("\t")
            length = int(length)
            chrom_length_dict[chrom] = length

    with GtfFile(path2) as f:
        records = [record for record in f]

    blacklist = set()
    for record in records:
        if record.end > chrom_length_dict[record.chrom]:
            blacklist.add(record.attributes["transcript_id"])
        if record.strand == ".":
            blacklist.add(record.attributes["transcript_id"])

    # sys.stderr.write(str(blacklist) + "\n")
    for record in records:
        line = record.format("GTF")
        if record.attributes["transcript_id"] in blacklist:
            sys.stderr.write(line + "\n")
        else:
            sys.stdout.write(line + "\n")


if __name__ == "__main__":
    main()
