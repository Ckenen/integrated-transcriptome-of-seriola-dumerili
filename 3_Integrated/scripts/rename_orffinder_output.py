#!/usr/bin/env python
import sys
import re
import pandas as pd
from Bio import SeqIO
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, FastaFile

gtf, fasta, infile, outfile1, outfile2 = sys.argv[1:]

with GtfFile(gtf) as f:
    transcripts = list(GtfTranscriptBuilder(f))
with FastaFile(fasta) as fasta:
    for transcript in transcripts:
        transcript.seq = fasta.fetch(obj=transcript)
transcript_dict = {t.name: t for t in transcripts}
rows = []
with open(outfile1, "w+") as fw:
    for record in SeqIO.parse(infile, "fasta"):
        name = record.name
        ret = re.match("^lcl\|ORF[0-9]+_", name)
        name, start, end = name[ret.end():].split(":")
        start = int(start) # 包含start codon和stop codon
        end = int(end) + 1
        transcript = transcript_dict[name]
        start_codon = transcript.seq[start:start + 3].upper()
        stop_codon = transcript.seq[end - 3:end].upper()
        iden = "%s:%d:%d" % (name, start, end)
        row = [iden, name, len(transcript), start, end, end - start, start_codon, stop_codon]
        fw.write(">%s\n%s\n" % (iden, str(record.seq)))
        rows.append(row)
dat = pd.DataFrame(rows)
dat.columns = ["ID", "Transcript", "Length", "ORF_Start", "ORF_End", "ORF_Length", "Start_Codon", "Stop_Codon"]
dat.to_csv(outfile2, sep="\t", index=False)