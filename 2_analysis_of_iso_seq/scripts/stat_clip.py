#!/usr/bin/env python
import sys
import optparse
from collections import defaultdict
import pysam


def get_clip_count(segment):
    count1 = 0 # 5'end
    count2 = 0 # 3'end
    cigars = segment.cigartuples
    flag, count = cigars[0]
    if flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
        count1 = count
    flag, count = cigars[-1]
    if flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
        count2 = count
    if segment.is_reverse:
        count1, count2 = count2, count1
    return count1, count2

        
def stat_clip():
    parser = optparse.OptionParser(usage="%prog [options] input.bam")
    parser.add_option("-o", "--outbam", dest="outbam", metavar="PATH", 
                      help="Output BAM. [%default]")
    parser.add_option("-s", "--summary", dest="summary", metavar="PATH", 
                      help="[%default]")
    parser.add_option("-c", "--max-clip-3", dest="max_clip", type="int", default=1000000000, metavar="INT", 
                      help="[%default]")
    options, args = parser.parse_args()
    
    assert len(args) == 1
    bamfile = args[0]
    outbam = options.outbam
    max_clip = options.max_clip
    f_summary = options.summary
    
    n1 = 0
    n2 = 0
    counter = defaultdict(int)
    with pysam.AlignmentFile(bamfile) as f:
        fw = None
        if outbam:
            fw = pysam.AlignmentFile(outbam, "wb", f)
        for segment in f:
            n1 += 1
            c1, c2 = get_clip_count(segment)
            counter[(c1, c2)] += 1
            # if fw and c1 <= max_clip and c2 <= max_clip:
            if fw and c2 <= max_clip:
                fw.write(segment)
                n2 += 1
    r = n2 / n1 if n1 > 0 else 0
    sys.stderr.write("Input\tOutput\tRatio\n")
    sys.stderr.write("%d\t%d\t%f\n" % (n1, n2, r))
            
    if f_summary:
        with open(f_summary, "w+") as fw:
            fw.write("HeadClip\tTailClip\tCount\tRatio\n")
            t = sum(counter.values())
            for (c1, c2), v in sorted(counter.items()):
                r = v / t
                fw.write("\t".join(map(str, [c1, c2, v, r])) + "\n")


if __name__ == "__main__":
    stat_clip()
    