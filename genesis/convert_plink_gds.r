#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

plink_prefix <- args[1]

library(SNPRelate)

snpgdsBED2GDS(bed.fn = paste0(plink_prefix, ".bed"), 
              bim.fn = paste0(plink_prefix, ".bim"), 
              fam.fn = paste0(plink_prefix, ".fam"), 
              out.gdsfn = paste0(plink_prefix, ".gds"))
