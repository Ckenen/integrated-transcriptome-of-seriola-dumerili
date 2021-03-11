#!/usr/bin/env python
import sys


def main():
    path = sys.argv[1]

    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            columns = line.strip("\n").split("\t")
            attributes = dict()
            for item in columns[8].split(";"):
                item = item.strip()
                if item != "":
                    i = item.find(" ")
                    assert i != -1
                    key = item[:i]
                    value = item[i + 1:]
                    value = value.strip("\"")
                    attributes[key] = value
            columns[-1] = attributes
            records.append(columns)
    
    transcripts = dict()
    for record in records:
        if record[2] == "exon":
            parent = record[-1]["transcript_id"]
            if parent not in transcripts.keys():
                transcripts[parent] = [[], []]
            transcripts[parent][0].append(record)
        if record[2] == "CDS":
            parent = record[-1]["transcript_id"]
            if parent not in transcripts.keys():
                transcripts[parent] = [[], []]
            transcripts[parent][1].append(record)


    for tid, (exons, cdss) in transcripts.items():
        chrom = None
        strand = None
        blocks1 = []
        blocks2 = []
        for exon in exons:
            if chrom is not None:
                assert exon[0] == chrom
            chrom = exon[0]
            if strand is not None:
                assert exon[6] == strand
            strand = exon[6]
            blocks1.append([int(exon[3]) - 1, int(exon[4])])
        for cds in cdss:
            if chrom is not None:
                assert cds[0] == chrom
            chrom = cds[0]
            if strand is not None:
                assert cds[6] == strand
            strand = cds[6]
            blocks2.append([int(cds[3]) - 1, int(cds[4])]) 
        if len(blocks1) == 0:
            assert len(blocks2) > 0
            blocks1 = blocks2
        blocks1 = list(sorted(blocks1, key=lambda item: item[0]))
        blocks2 = list(sorted(blocks2, key=lambda item: item[0]))
        start = blocks1[0][0]
        end = blocks1[-1][1]
        sizes = [block[1] - block[0] for block in blocks1]
        offsets = [block[0] - start for block in blocks1]
        tstart = start
        tend = start
        if len(blocks2) > 0:
            tstart, tend = blocks2[0][0], blocks2[-1][1]
        vs = [chrom, start, end, tid, ".", strand, tstart, tend, "0,0,0", len(blocks1), ",".join(map(str, sizes)) + ",", ",".join(map(str, offsets)) + ","]
        line = "\t".join(map(str, vs))
        print(line)


if __name__ == "__main__":
    main()
