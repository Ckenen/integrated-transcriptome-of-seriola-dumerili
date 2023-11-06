#!/usr/bin/env python
import sys
import pysam

def main():
    infile, name, outfile = sys.argv[1:]
    with pysam.AlignmentFile(infile) as f:
        header = f.header.as_dict()
        header["RG"] = [{"ID": name, "SM": name, "LB": name}]
        with pysam.AlignmentFile(outfile, "wb", header=header) as fw:
            for segment in f:
                segment.query_name = "%s.%s" % (name, segment.query_name)
                segment.set_tag("RG", name)
                fw.write(segment)


if __name__ == "__main__":
    main()