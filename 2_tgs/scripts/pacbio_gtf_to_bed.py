#!/usr/bin/env python
import sys
import gzip

def main():
    path = sys.argv[1]

    # Load records
    records = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            values = line.strip("\n").split("\t")
            assert len(values) == 9
            attributes = dict()
            for item in values[-1].split(";"):
                item = item.strip()
                if item == "":
                    continue
                key, value = item.split()
                value = value[1:-1]
                attributes[key] = value
            values[-1] = attributes
            records.append(values)

    # Cluster exons
    exons = dict()
    for record in records:
        if record[2] == "exon":
            tid = record[-1]["transcript_id"]
            if tid not in exons:
                exons[tid] = list()
            exons[tid].append(record)

    # Generate BED format
    for tid, items in exons.items():
        chrom = None
        strand = None
        blocks = []
        for exon in items:
            if chrom is not None:
                assert exon[0] == chrom
            chrom = exon[0]
            if strand is not None:
                assert exon[6] == strand
            strand = exon[6]
            blocks.append([int(exon[3]) - 1, int(exon[4])])
        blocks = list(sorted(blocks, key=lambda item: item[0]))
        start = blocks[0][0]
        end = blocks[-1][1]
        sizes = ",".join(map(str, [by - bx for bx, by in blocks])) + ","
        offsets = ",".join(map(str, [bx - start for bx, by in blocks])) + ","
        vs = [chrom, start, end, tid, ".", strand, start, start, "0,0,0", len(blocks), sizes, offsets]
        s = "\t".join(map(str, vs))
        print(s)


if __name__ == "__main__":
    main()