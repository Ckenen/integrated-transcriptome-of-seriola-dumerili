#!/usr/bin/env python
import sys
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import GtfFile, GtfGeneBuilder


def main():
    # input.gtf orf.bed lncrna.bed
    path1, path2, path3 = sys.argv[1:]

    with GtfFile(path1) as f:
        records = [x for x in f]
    genes = list(GtfGeneBuilder(records))
    array1 = []
    array2 = []
    for gene in genes:
        if gene.chrom == "NC_016870.1":
            continue
        temp1 = []
        temp2 = []
        for transcript in gene.transcripts:
            biotype = transcript.records["transcript"][0].attributes["gene_biotype"]
            if biotype == "protein_coding" and transcript.thick:
                temp1.append(transcript)
            elif biotype == "lncRNA":
                temp2.append(transcript)
        temp1 = list(sorted(temp1, key=lambda item: len(item)))
        temp2 = list(sorted(temp2, key=lambda item: len(item)))
        if len(temp1) > 0:
            assert len(temp2) == 0
            transcript = temp1[-1]
            thick_start, thick_end = transcript.thick
            obj = transcript.clip(thick_start, thick_end)
            array1.append(obj)
        elif len(temp2) > 0:
            assert len(temp1) == 0
            transcript = temp2[-1]
            array2.append(transcript)
    with open(path2, "w+") as fw:
        for transcript in sorted(array1):
            fw.write(transcript.format("BED") + "\n")
    with open(path3, "w+") as fw:
        for transcript in sorted(array2):
            fw.write(transcript.format("BED") + "\n")


if __name__ == "__main__":
    main()
