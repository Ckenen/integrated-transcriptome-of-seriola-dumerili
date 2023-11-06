#!/usr/bin/env python
import sys
import re
from pyBioInfo.IO.File import GtfFile


def main():
    path1, path2 = sys.argv[1:]
    with GtfFile(path1) as f, open(path2, "w+") as fw:
        for record in f:
            transcript_id = record.attributes["transcript_id"]
            assert re.match("^PB\.[0-9]+\.[0-9]+$", transcript_id)
            gene_id = ".".join(transcript_id.split(".")[:2])
            record.attributes["gene_id"] = gene_id
            fw.write(record.format() + "\n")


if __name__ == "__main__":
    main()
