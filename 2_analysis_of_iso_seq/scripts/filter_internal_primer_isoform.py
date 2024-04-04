#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import BedFile, GtfFile, GtfTranscriptBuilder, FastaFileRandom


def is_internal_primer(seq): 
    if seq.count("A") / len(seq) > 0.75:
        return True
    for i in range(0, 5):
        if seq[i:i + 10].count("A") >= 8:
            return True
    return False


def main():

    f_gtf_bed, f_fasta = sys.argv[1:]
    
    if "gtf" in f_gtf_bed:
        f_gtf = f_gtf_bed

        with GtfFile(f_gtf) as f:
            records = [x for x in f]
        transcripts = list(GtfTranscriptBuilder(records))
        transcripts = list(filter(lambda item: len(item) >= 200, transcripts))
        transcripts = list(filter(lambda item: item.chrom != "NC_016870.1", transcripts))

        with FastaFileRandom(f_fasta) as fasta:
            for transcript in transcripts:
                if transcript.strand == "+":
                    start = transcript.end
                else:
                    start = transcript.start - 20
                end = start + 20
                seq = fasta.fetch(chrom=transcript.chrom, start=start, end=end, strand=transcript.strand)
                transcript.polya_downstream_sequence = seq.upper()

        for transcript in transcripts:
            transcript.is_internal_primer = is_internal_primer(transcript.polya_downstream_sequence)

        transcripts = list(filter(lambda item: not item.is_internal_primer, transcripts))

        records = []
        for transcript in transcripts:
            for key, values in transcript.records.items():
                for item in values:
                    records.append(item)
        records = list(sorted(records))

        for record in records:
            print(record.format())
            
    elif "bed" in f_gtf_bed:
        f_bed = f_gtf_bed
        
        transcripts = list(BedFile(f_bed))
        
        with FastaFileRandom(f_fasta) as fasta:
            for transcript in transcripts:
                if transcript.strand == "+":
                    start = transcript.end
                else:
                    start = transcript.start - 20
                end = start + 20
                seq = fasta.fetch(chrom=transcript.chrom, start=start, end=end, strand=transcript.strand)
                transcript.polya_downstream_sequence = seq.upper()

        for transcript in transcripts:
            transcript.is_internal_primer = is_internal_primer(transcript.polya_downstream_sequence)

        transcripts = list(filter(lambda item: not item.is_internal_primer, transcripts))
        
        for t in transcripts:
            print(t.format())


if __name__ == "__main__":
    main()