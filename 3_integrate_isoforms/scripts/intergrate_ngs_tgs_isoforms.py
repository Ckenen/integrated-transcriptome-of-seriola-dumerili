#!/usr/bin/env python
import sys
import gzip
import pandas as pd


def main():
    ngs_gtf, tgs_gtf, ngs_fpkm, tgs_fpkm, sqanti_class, out_gtf = sys.argv[1:]

    fpkms_ngs = pd.read_csv(ngs_fpkm, sep="\t", index_col=0)
    fpkms_tgs = pd.read_csv(tgs_fpkm, sep="\t", index_col=0)
    sqanti = pd.read_csv(sqanti_class, sep="\t", index_col=0)

    NGS_FPKM_CUTOFF = 0.5
    ngs = fpkms_ngs[fpkms_ngs.max(axis=1) >= NGS_FPKM_CUTOFF].index.values
    print("Loaded %d NGS-based isoforms, %d isoforms with fpkm >= %.2f" %
          (len(fpkms_ngs), len(ngs), NGS_FPKM_CUTOFF))

    intergenic = sqanti[sqanti["structural_category"]
                        == "intergenic"].index.values
    print("%s NGS-bases isoforms were intergenic to TGS-based isoforms" %
          len(intergenic))

    ngs = set(ngs) & set(intergenic)
    print("Finally retain %s NGS-based isoforms" % len(ngs))

    TGS_FPKM_CUTOFF = 0
    tgs = fpkms_tgs[fpkms_tgs.max(axis=1) >= TGS_FPKM_CUTOFF].index.values
    print("Loaded %s TGS-based isoforms, %d isoforms with fpkm >= %.2f" %
          (len(fpkms_tgs), len(tgs), TGS_FPKM_CUTOFF))

    with open(out_gtf, "w+") as fw:
        with gzip.open(ngs_gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                row = line.strip("\n").split("\t")
                attris = dict()
                for item in row[-1].split(";"):
                    item = item.strip()
                    if item == "":
                        continue
                    k, v = item.split(" ")
                    attris[k.strip()] = v.strip()[1:-1]
                tid = attris["transcript_id"]
                if tid in ngs:
                    fw.write(line)

        with gzip.open(tgs_gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fw.write(line)


if __name__ == "__main__":
    main()
