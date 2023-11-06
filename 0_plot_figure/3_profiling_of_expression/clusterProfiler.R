#!/usr/bin/env Rscript
library("clusterProfiler")
args <- commandArgs(TRUE)
infile <- args[1]
prefix <- args[2]

gene <- read.table(infile, sep = "\t", header = T)$GeneID
term2gene <- read.csv("/data/chenzonggui/gaotishi/3_integrate_isoforms/results/function/goterm/term2gene.tsv", header = F, sep = "\t")
term2name <- read.csv("/data/chenzonggui/gaotishi/3_integrate_isoforms/results/function/goterm/term2name.tsv", header=F, sep = "\t")

x <- enricher(gene, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)

write.csv(x, paste0(prefix, ".enriched_terms.csv", sep = ""))

pdf(paste0(prefix, ".barplot.pdf", sep = ""), width = 6, height = 8)
barplot(x)
dev.off()

pdf(paste0(prefix, ".dotplot.pdf", sep = ""), width = 6, height = 8)
dotplot(x)
dev.off()