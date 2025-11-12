#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

gds_prefix <- args[1]  # path to file with name but without extension

library(SNPRelate)

gds <- snpgdsOpen(paste0(gds_prefix, ".gds"))

snpset <- snpgdsLDpruning(gds,
                          method="corr",
                          slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1),
                          verbose=FALSE)

pruned <- unlist(snpset,
                 use.names=FALSE)

write.table(pruned,
            file = paste0(gds_prefix, "_prune.in"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
snpgdsClose(gds)
