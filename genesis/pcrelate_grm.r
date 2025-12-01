#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

plink_prefix <- args[1]  # path to file with name but without extension
pruned_list <- args[2]
num_pcs <- as.numeric(args[3])
small_sample <- args[4]  # should be "yes" if < 5000, otherwise "no" 

library(SNPRelate)
library(GENESIS)
library(GWASTools)

geno <- GdsGenotypeReader(paste0(plink_prefix, ".gds"))
genoData <- GenotypeData(geno)
pruned <- read.table(pruned_list)[,1]
mypcair <- readRDS(paste0(plink_prefix, "_pcair.rds"))
small_sample_flag <- small_sample == "yes"

# PC-Relate
genoData <- GenotypeBlockIterator(genoData, snpInclude=pruned)

print(paste0("Running PCRelate with ", num_pcs, " PCs."))
mypcrelate <- pcrelate(genoData,
                       pcs = mypcair$vectors[,1:num_pcs], 
                       training.set = mypcair$unrels,
                       BPPARAM = BiocParallel::MulticoreParam(4),
                       ibd.probs = TRUE,
                       sample.block.size = 5000, # default
                       small.samp.correct = small_sample_flag,
                       maf.thresh = 0.01) # default
close(genoData)
saveRDS(mypcrelate, paste0(plink_prefix, "_pcrelate.rds"))

# Export GRM
grm <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2))
saveRDS(grm, file = paste0(plink_prefix, "_GRM.rds"))

# Check if kinship matrix has negative eigenvalues, if so, correct to
# the nearest positive definite matrix
e_vals <- eigen(grm, symmetric = T, only.values = T)
e_vals_neg <- sum(e_vals$values < 0)

if (e_vals_neg > 0) {
    print(paste0("There are ", e_vals_neg, " negative eigenvalues:"))
    print(e_vals$values[e_vals$values < 0])
    print("These will be corrected with Matrix::nearPD, but it should")
    print("be confirmed that it seems reasonble to do that.")

    grm_pd <- Matrix::nearPD(grm)
    saveRDS(grm_pd, file = paste0(plink_prefix, "_GRM_pd.rds"))
} else {
   print("There are no negative eigenvalues, GRM is positive definite.")
}