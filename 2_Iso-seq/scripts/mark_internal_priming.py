#!/usr/bin/env python
import sys
import pysam

MAPPER = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "-": "-"
}


def rev_com(seq):
    tmp = []
    i = len(seq) - 1
    while i >= 0:
        tmp.append(MAPPER[seq[i]])
        i = i - 1
    return "".join(tmp)


def is_internal_primer(seq):
    if seq.count("A") / len(seq) >= 0.8:
        return True
    for i in range(0, 5):
        if seq[i:i + 10].count("A") >= 8:
            return True
    return False


def main():
    infile1, infile2, outfile = sys.argv[1:]

    bam = pysam.AlignmentFile(infile1)
    fasta = pysam.FastaFile(infile2)
    out = pysam.AlignmentFile(outfile, "wb", bam)

    WIDTH = 20

    for segment in bam:
        if segment.is_unmapped:
            continue
        chrom = segment.reference_name
        chrom_length = bam.get_reference_length(chrom)
        strand = "-" if segment.is_reverse else "+"
        if strand == "+":
            p1 = segment.reference_end - WIDTH
        else:
            p1 = segment.reference_start - WIDTH
        p2 = p1 + 2 * WIDTH

        seq = fasta.fetch(chrom, max(0, p1), min(p2, chrom_length)).upper()
        if p1 < 0:
            seq = "-" * abs(p1) + seq
        if p2 > chrom_length:
            seq = seq + "-" * (p2 - chrom_length)
        if strand == "-":
            seq = rev_com(seq)

        priming = is_internal_primer(seq[WIDTH:])
        segment.set_tag("XW", WIDTH)
        segment.set_tag("XU", seq[:WIDTH])
        segment.set_tag("XD", seq[WIDTH:])
        segment.set_tag("XP", priming)
        out.write(segment)

    out.close()
    bam.close()
    fasta.close()


if __name__ == "__main__":
    main()
