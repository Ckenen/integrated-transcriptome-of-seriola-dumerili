#!/usr/bin/env Rscript
library(DESeq2)

args<-commandArgs(T)
infile <- args[1]
outfile <- args[2]
counts <- read.table(infile, sep = "\t", header = T, row.names = 1)
condition = factor(rep(c('control', 'treat'), each = 2))
coldata <- data.frame(condition = factor(rep(c('control', 'treat'), each = 2), levels = c('control', 'treat')))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'treat', 'control'))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
m <- merge(x = counts, y = res1, by = "row.names")
write.table(m, outfile, row.names = F, sep = '\t', quote = FALSE)


