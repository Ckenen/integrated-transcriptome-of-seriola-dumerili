#!/usr/bin/env python
import sys
import pandas as pd


def main():
    path1, path2, path3 = sys.argv[1:]

    gids = dict()
    biotypes = dict()
    gnames = dict()
    with open(path1, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            columns = line.strip("\n").split("\t")
            attributes = dict()
            for item in columns[8].split(";"):
                item = item.strip()
                if item != "":
                    i = item.find(" ")
                    if i == -1:
                        print(item)
                    assert i != -1
                    key = item[:i]
                    value = item[i + 1:]
                    value = value.strip("\"")
                    attributes[key] = value
            if columns[2] == "transcript":
                tid = attributes.get("transcript_id", "")
                gid = attributes.get("gene_id", "")
                biotype = attributes.get("gene_biotype", "")
                gname = attributes.get("gene_name", "")
                gids[tid] = gid
                biotypes[tid] = biotype
                gnames[tid] = gname

    headers = ["TranscriptID", "GeneID", "GeneName",
               "BioType", "Chrom", "Start", "End", "Strand", "Length"]
    rows = []
    with open(path2) as f:
        for line in f:
            columns = line.strip("\n").split("\t")
            tid = columns[3]
            gid = gids.get(tid, "")
            gname = gnames.get(tid, "")
            biotype = biotypes.get(tid, "")
            chrom = columns[0]
            start = columns[1]
            end = columns[2]
            strand = columns[5]
            length = sum([int(x) for x in columns[10].split(",")[:-1]])
            vs = [tid, gid, gname, biotype, chrom, start, end, strand, length]
            rows.append(vs)
    anno = pd.DataFrame(rows, columns=headers)

    tids = set()
    for gid, dat in anno.groupby(by="GeneID"):
        dat = dat.sort_values(by="Length", ascending=False)
        tid = dat["TranscriptID"].values[0]
        tids.add(tid)
    anno["Longest"] = [tid in tids for tid in anno["TranscriptID"]]

    anno.to_csv(path3, sep="\t", index=False)


if __name__ == "__main__":
    main()
