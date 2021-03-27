#!/usr/bin/env python
import sys
from collections import defaultdict
from pyBioInfo.IO.File import GffFile

# deprecated
assert False

def main():
    path1, path2 = sys.argv[1:]
    with GffFile(path1) as f:
        records = [x for x in f]

    id_to_record_dict = dict()
    parent_to_childs_dict = defaultdict(list)
    for record in records:
        id_to_record_dict[record.attributes["ID"]] = record
        if "Parent" in record.attributes:
            parent_to_childs_dict[record.attributes["Parent"]].append(record)

    transcript_id_set = set()
    for record in records:
        if record.feature == "exon" or record.feature == "CDS":
            transcript_id_set.add(record.attributes["Parent"])

    with open(path2, "w+") as fw:
        gene_id_set = set()
        for transcript_id in transcript_id_set:
            transcript_record = id_to_record_dict[transcript_id]

            if "gene" in transcript_record.attributes:
                gene_name = transcript_record.attributes["gene"]
            else:
                gene_name = transcript_record.attributes["product"]
                # print(gene_name)
            if "Parent" in transcript_record.attributes:
                gene_record = id_to_record_dict[transcript_record.attributes["Parent"]]
                gene_id = gene_record.attributes["ID"]
            else:
                gene_id = transcript_id
            gene_id_set.add(gene_id)

            attr = "gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"" % (
                gene_id, transcript_id, gene_name)
            vs = [
                transcript_record.chrom, transcript_record.source, "transcript",
                transcript_record.start + 1, transcript_record.end, transcript_record.score,
                transcript_record.strand, ".", attr
            ]
            line = "\t".join(map(str, vs))
            fw.write(line + "\n")

            childs = parent_to_childs_dict[transcript_id]
            assert len(childs) > 0
            items1 = list(filter(lambda x: x.feature == "exon", childs))
            items2 = list(filter(lambda x: x.feature == "CDS", childs))
            if len(items1) == 0:
                items1 = items2
            for item in items1:
                attr = "gene_id \"%s\"; transcript_id \"%s\";" % (
                    gene_id, transcript_id)
                vs = [item.chrom, item.source, "exon", item.start + 1, item.end,
                      item.score, item.strand, item.frame if item.frame else ".", attr]
                line = "\t".join(map(str, vs))
                fw.write(line + "\n")
            for item in items2:
                attr = "gene_id \"%s\"; transcript_id \"%s\";" % (
                    gene_id, transcript_id)
                vs = [item.chrom, item.source, "CDS", item.start + 1, item.end,
                      item.score, item.strand, item.frame if item.frame else ".", attr]
                line = "\t".join(map(str, vs))
                fw.write(line + "\n")

        for gene_id in gene_id_set:
            gene_record = id_to_record_dict[gene_id]

            if "gene" in gene_record.attributes:
                gene_name = gene_record.attributes["gene"]
            else:
                gene_name = gene_record.attributes["product"]

            if "gene_biotype" in gene_record.attributes:
                biotype = gene_record.attributes["gene_biotype"]
            else:
                biotype = gene_record.feature

            items = [
                "gene_id \"%s\"" % gene_id,
                "gene_name \"%s\"" % gene_name,
                "gene_biotype \"%s\"" % biotype
            ]
            attri = "; ".join(items) + ";"
            vs = [
                gene_record.chrom, gene_record.source, "gene", gene_record.start + 1, gene_record.end,
                gene_record.score, gene_record.strand, ".", attri
            ]
            fw.write("\t".join(map(str, vs)) + "\n")

if __name__ == "__main__":
    main()
