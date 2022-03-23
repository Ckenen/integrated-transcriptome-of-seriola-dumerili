#!/usr/bin/env python
import os
import sys
import pandas as pd
from pyBioInfo.IO.File import GtfFile

path_list = sys.argv[1:-1]
outfile = sys.argv[-1]

array = []
for path in path_list:
    name = os.path.basename(path)
    transcript_id_list = []
    fpkm_list = []
    with GtfFile(path + "/transcripts.gtf") as f:
        for record in f:
            if record.feature != "transcript":
                continue
            transcript_id_list.append(record.attributes["transcript_id"])
            fpkm_list.append(record.attributes["FPKM"])
    dat = pd.DataFrame(index=transcript_id_list, data={name: fpkm_list})
    dat.index.name = "TranscriptID"
    array.append(dat)
dat = pd.concat(array, axis=1)
dat.index.name = "TranscriptID"
dat.to_csv(outfile, sep="\t")
