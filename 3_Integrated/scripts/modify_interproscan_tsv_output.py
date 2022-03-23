#!/usr/bin/env python
import sys

columns = [
    "Accession", "MD5", "Length", "Analysis", 
    "Signature_Accession", "Signature_Description", 
    "Match_Start", "Match_End", "Evalue", "Status", "Date", 
    "InterPro_Accession", "InterPro_Description", "GO", "Pathways"]
print("\t".join(columns))

with open(sys.argv[1]) as f:
    for line in f:
        row = line.strip("\n").split("\t")
        while len(row) < 15:
            row.append("-")
        print("\t".join(row))