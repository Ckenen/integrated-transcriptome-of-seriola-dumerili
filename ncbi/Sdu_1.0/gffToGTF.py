#!/usr/bin/env python
import sys
from collections import defaultdict


def main():
    path1, path2 = sys.argv[1:]

    records = []
    with open(path1) as f:
        for line in f:
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            record = line.split("\t")
            attris = dict()
            for item in record[-1].split(";"):
                item = item.strip()
                n = item.find("=")
                key = item[:n]
                val = item[n + 1:]
                attris[key] = val
            record[-1] = attris
            records.append(record)

    mapper = dict()
    for record in records:
        mapper[record[-1]["ID"]] = record

    array1 = defaultdict(list)
    array2 = defaultdict(list)
    for record in records:
        if record[2] == "exon":
            array1[record[-1]["Parent"]].append(record)
        elif record[2] == "CDS":
            array2[record[-1]["Parent"]].append(record)

    tids1 = list(array1.keys())
    tids2 = list(array2.keys() - array1.keys())

    with open(path2, "w+") as fw:
        for tid in tids1:
            list1 = array1[tid]
            list2 = array2[tid]
            transcript = mapper[tid]
            gene = None
            if "Parent" in transcript[-1].keys():
                gene_id = transcript[-1]["Parent"]
                gene = mapper[gene_id]
            else:
                gene_id = tid

            values = transcript[:-1]
            values[2] = "transcript"
            attri1 = "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, tid)
            s1 = attri1
            attri2 = transcript[-1].copy()
            if gene:
                for k, v in gene[-1].items():
                    attri2[k] = v
            tmp = []
            for k in attri2.keys():
                v = attri2[k]
                if k == "ID" or k == "Parent" or k == "transcript_id" or k == "gene_id":
                    continue
                if k == "gene":
                    k = "gene_name"
                v = v.replace("%2C", ",")
                v = v.replace("%25", "%")
                v = v.replace("%3B", ",")
                tmp.append("%s \"%s\"" % (k, v))
            if len(tmp) > 0:
                s1 = s1 + " " + "; ".join(tmp) + ";"
            values.append(s1)
            fw.write("\t".join(values) + "\n")
            for exon in list1:
                assert exon[2] == "exon"
                values = exon[:-1]
                values.append(
                    "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, tid))
                fw.write("\t".join(values) + "\n")
            for cds in list2:
                assert cds[2] == "CDS"
                values = cds[:-1]
                values.append(
                    "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, tid))
                fw.write("\t".join(values) + "\n")

        for tid in tids2:
            list1 = array1[tid]
            assert len(list1) == 0
            list2 = array2[tid]
            gene = mapper[tid]
            values = gene[:-1]
            values[2] = "transcript"
            s1 = "gene_id \"%s\"; transcript_id \"%s\";" % (
                gene[-1]["ID"], tid)
            tmp = []
            for k in gene[-1].keys():
                v = gene[-1][k]
                if k == "ID" or k == "Parent" or k == "transcript_id" or k == "gene_id":
                    continue
                if k == "gene":
                    k = "gene_name"
                v = v.replace("%2C", ",")
                v = v.replace("%25", "%")
                v = v.replace("%3B", ",")
                tmp.append("%s \"%s\"" % (k, v))
            s1 = s1 + " " + "; ".join(tmp) + ";"
            values.append(s1)
            fw.write("\t".join(values) + "\n")
            for cds in list2:
                assert cds[2] == "CDS"
                values = cds[:-1]
                values.append(
                    "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, tid))
                fw.write("\t".join(values) + "\n")
                values = cds[:-1]
                values[2] = "exon"
                values.append(
                    "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, tid))
                fw.write("\t".join(values) + "\n")


if __name__ == "__main__":
    main()
