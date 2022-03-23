#!/usr/bin/env python
import os
import sys
import multiprocessing
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder

def run(options):
    cmd = "cpat.py -r %s -g %s -d %s -x %s -o %s > %s 2>&1" % (options["fasta"], 
                                                     options["bed"], 
                                                     options["rda"], 
                                                     options["hexamer"], 
                                                     options["prefix"],
                                                     options["log"])
    print(cmd)
    assert os.system(cmd) == 0

gtf = sys.argv[1]
fasta = sys.argv[2]
rda = sys.argv[3]
hexamer = sys.argv[4]
threads = int(sys.argv[5])
outdir = sys.argv[6]

os.mkdir(outdir)

with GtfFile(gtf) as f:
    records = [x for x in f]
transcripts = list(GtfTranscriptBuilder(records))

# 生成任务
n = 0
array = []
capacity = 500
for i in range(0, len(transcripts), capacity):
    j = min(i + capacity, len(transcripts))
    if j <= i:
        break
    options = {
        "batch": n,
        "fasta": "../common/ncbi/serDum.ncbi.fasta",
        "bed": outdir + "/batch.%d.bed" % n,
        "rda": rda,
        "hexamer": hexamer,
        "prefix": outdir + "/batch.%d" % n,
        "log": outdir + "/batch.%d.log" % n,
        "best": outdir + "/batch.%d.ORF_prob.best.tsv" % n,
        "tsv": outdir +  "/batch.%d.ORF_prob.tsv" % n,
        "no_orf": outdir + "/batch.%d.no_ORF.txt" % n,
        "seqs": outdir + "/batch.%d.ORF_seqs.fa" % n,
        "r": outdir + "/batch.%d.r" % n
    }
    with open(options["bed"], "w+") as fw:
        for transcript in transcripts[i:j]:
            fw.write(transcript.format("BED") + "\n")
    array.append(options)
    n += 1

# 执行任务
pool = multiprocessing.Pool(threads)
for options in array:    
    pool.apply_async(run, (options, ))
pool.close()
pool.join()

# 收集结果
fw1 = open(outdir + "/all.no_ORF.txt", "w+")
fw2 = open(outdir + "/all.ORF_prob.best.tsv", "w+")
fw3 = open(outdir + "/all.ORF_prob.tsv", "w+")
fw4 = open(outdir + "/all.ORF_seqs.fa", "w+")
for n, options in enumerate(array):
    with open(options["no_orf"]) as f:
        for line in f:
            fw1.write(line)
    with open(options["best"]) as f:
        for i, line in enumerate(f):
            if i == 0:
                if n == 0:
                    fw2.write(line)
            else:
                fw2.write(line)
    with open(options["tsv"]) as f:
        for i, line in enumerate(f):
            if i == 0:
                if n == 0:
                    fw3.write(line)
            else:
                fw3.write(line)
    with open(options["seqs"]) as f:
        for line in f:
            fw4.write(line)  
    os.remove(options["no_orf"])
    os.remove(options["bed"])
    os.remove(options["best"])
    os.remove(options["tsv"])
    os.remove(options["seqs"])
    os.remove(options["r"])
    os.remove(options["log"])
fw1.close()
fw2.close()
fw3.close()
fw4.close()
