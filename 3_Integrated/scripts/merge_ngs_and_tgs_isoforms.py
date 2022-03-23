#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder
from pyBioInfo.Utils import ShiftLoader

with GtfFile(sys.argv[1]) as f:
    transcript_list_ngs = list(filter(lambda item: item.chrom != "NC_016870.1", GtfTranscriptBuilder(f)))

with GtfFile(sys.argv[2]) as f:
    transcript_list_tgs = list(filter(lambda item: item.chrom != "NC_016870.1", GtfTranscriptBuilder(f)))

loader = ShiftLoader(transcript_list_tgs)
array = []
for transcript in transcript_list_ngs:
    flag = True
    for obj in loader.fetch(obj=transcript):
        if obj.strand == transcript.strand:
            flag = False
            break
    if flag:
        array.append(transcript)
transcripts = list(sorted(array + transcript_list_tgs))

records = []
for transcript in transcripts:
    for values in transcript.records.values():
        records.extend(values)
for record in list(sorted(records)):
    print(record.format())