#!/usr/bin/env python
import sys
import pysam


def main():
    in_bam, out_bam = sys.argv[1:]

    with pysam.AlignmentFile(in_bam) as f, pysam.AlignmentFile(out_bam, "wb", f) as fw:
        for segment in f:
            if segment.is_unmapped:
                continue
            if not segment.is_proper_pair:
                continue
            if segment.is_secondary:
                continue
            tags = [tag[0] for tag in segment.get_tags()]
            if "XA" in tags:
                continue
            if "SA" in tags:
                continue
            fw.write(segment)


if __name__ == "__main__":
    main()
