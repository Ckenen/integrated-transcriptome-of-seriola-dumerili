#!/usr/bin/env python
import sys
import pysam
from pyBioInfo.IO.File import FastaFile


def main():
    # input.bam genome.fasta output.bam
    path1, path2, path3 = sys.argv[1:]

    print("Chrom\tStart\tEnd\tName\tStrand\tSequence\tPass")
    
    with pysam.AlignmentFile(path1) as f, FastaFile(path2) as fasta, pysam.AlignmentFile(path3, "wb", template=f) as fw:
        for segment in f:
            if segment.is_unmapped:
                continue
            chrom = segment.reference_name
            if chrom == "NC_016870.1":
                continue
            start = segment.reference_start
            end = segment.reference_end
            name = segment.query_name
            strand = "-" if segment.is_reverse else "+"
            flag = True
            if strand == "+":
                start1 = segment.reference_end
            else:
                start1 = segment.reference_start - 20
            end1 = start1 + 20
            
            seq = "X" * 20
            try:
                seq = fasta.fetch(chrom=chrom, start=start1,
                                  end=end1, strand=strand)
                seq = seq.upper()
                if seq.count("A") >= len(seq) * 0.75:
                    flag = False
                if len(seq) < 20:
                    seq = seq + "X" * (20 - len(seq))
            except ValueError:
                pass
            row = [chrom, start, end, name, strand, seq, flag]
            print("\t".join(map(str, row)))
            if flag:
                fw.write(segment)


if __name__ == "__main__":
    main()
