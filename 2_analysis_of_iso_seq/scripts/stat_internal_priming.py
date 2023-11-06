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

MAX_THREE_END_CLIP = 5


def reverse_complement(seq):
    tmp = []
    i = len(seq) - 1
    while i >= 0:
        tmp.append(MAPPER[seq[i]])
        i = i - 1
    return "".join(tmp)


def is_internal_priming(seq):
    if seq.count("A") / len(seq) >= 0.8:
        return True
    for i in range(0, 5):
        if seq[i:i + 10].count("A") >= 8:
            return True
    return False


def main():
    f_bam, f_fa, f_out = sys.argv[1:]

    WIDTH = 20
    
    with pysam.AlignmentFile(f_bam) as bam, \
        pysam.FastaFile(f_fa) as fasta, \
        pysam.AlignmentFile(f_out, "wb", bam) as out:

        for segment in bam:
            if segment.is_unmapped:
                continue
            chrom = segment.reference_name
            chrom_length = bam.get_reference_length(chrom)
            strand = "-" if segment.is_reverse else "+"
            if strand == "+":
                start = segment.reference_end - WIDTH
            else:
                start = segment.reference_start - WIDTH
            end = start + 2 * WIDTH

            seq = fasta.fetch(chrom, max(0, start), min(end, chrom_length)).upper()
            if start < 0:
                seq = "-" * abs(start) + seq
            if end > chrom_length:
                seq = seq + "-" * (end - chrom_length)
            if strand == "-":
                seq = reverse_complement(seq)

            priming = is_internal_priming(seq[WIDTH:])
            segment.set_tag("XW", WIDTH)
            segment.set_tag("XU", seq[:WIDTH])
            segment.set_tag("XD", seq[WIDTH:])
            segment.set_tag("XP", priming)
            out.write(segment)


if __name__ == "__main__":
    main()
