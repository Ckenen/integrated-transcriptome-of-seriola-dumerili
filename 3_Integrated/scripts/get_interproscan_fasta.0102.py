#!/usr/bin/env python
import sys
import re
import pandas as pd
from Bio import SeqIO


diamond_out, orf_fasta, out_fasta = sys.argv[1:]

dat = pd.read_csv(diamond_out, sep="\t")
isoforms = []
for iden in dat["ID"]:
    r = re.search("^lcl\|ORF[0-9]+\_", iden)
    assert r is not None
    i1, i2 = r.span()
    tid = iden[i2:].split(":")[0]
    isoforms.append(tid)
dat["Isoform"] = isoforms
dat = dat[dat["Evalue"] <= 1e-6]
rows = []
for tid, tmp in dat.groupby("Isoform"):
    tmp = tmp.sort_values(by="Evalue")
    rows.append(tmp.iloc[0])
dat = pd.DataFrame(rows)
idens = set(dat["ID"])


with open(out_fasta, "w+") as fw:
    for record in SeqIO.parse(orf_fasta, "fasta"):
        if record.id in idens:
            SeqIO.write(record, fw, "fasta")
