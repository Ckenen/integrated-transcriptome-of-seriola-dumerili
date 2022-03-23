#!/usr/bin/env python
import sys
import re
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, FastaFile

gtf, fasta, infile, outfile = sys.argv[1:]
with GtfFile(gtf) as f:
    transcripts = list(GtfTranscriptBuilder(f))
with FastaFile(fasta) as fasta:
    for transcript in transcripts:
        transcript.seq = fasta.fetch(obj=transcript)
transcript_dict = {t.name: t for t in transcripts}

dat = pd.read_csv(infile, sep="\t")
index_list = []
transcript_name_list = []
rna_start_list = []
rna_end_list = []
start_codon_list = []
stop_codon_list = []
for v, start, end in dat[["ID", "ORF_start", "ORF_end"]].values:
    ret = re.search("_ORF_[0-9]+$", v)
    name = v[:ret.start()]
    start = start - 1
    index = "%s:%d:%d" % (name, start, end)
    start_codon = transcript_dict[name].seq[start:start + 3].upper()
    stop_codon = transcript_dict[name].seq[end - 3:end].upper()
    index_list.append(index)
    transcript_name_list.append(name)
    rna_start_list.append(start)
    rna_end_list.append(end)
    start_codon_list.append(start_codon)
    stop_codon_list.append(stop_codon)
dat.index = index_list
dat.index.name = "Index"
dat["Transcript"] = transcript_name_list
dat["ORF_Start"] = rna_start_list
dat["ORF_End"] = rna_end_list
dat["Start_Codon"] = start_codon_list
dat["Stop_Codon"] = stop_codon_list
dat = dat[["Transcript", "mRNA", "ORF_Start", "ORF_End", "ORF", "Start_Codon", "Stop_Codon", "Fickett", "Hexamer", "Coding_prob"]]
dat.columns = ["Transcript", "Length", "ORF_Start", "ORF_End", "ORF_Length", "Start_Codon", "Stop_Codon", "Fickett", "Hexamer", "Coding_Prob"]
dat.index.name = "ID"
dat.to_csv(outfile, sep="\t")