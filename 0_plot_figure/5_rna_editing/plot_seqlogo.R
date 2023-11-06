#!/usr/bin/env Rscript
library(seqLogo)

args <- commandArgs(TRUE)
infile <- args[1]
prefix <- args[2]

m <- read.csv(infile, row.names=1)
p <- makePWM(m)

pdf(paste0(prefix, ".pdf"), width=8, height=4)
seqLogo(p, ic.scale=F)
dev.off()

pdf(paste0(prefix, ".scaled.pdf"), width=8, height=4)
seqLogo(p, ic.scale=T)
dev.off()