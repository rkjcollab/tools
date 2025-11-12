#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

plink_prefix_pruned <- args[1]  # path to file with name but without extension
plink_prefix <- args[2]  # path to file with name but without extension
pruned_list <- args[3]

library(SNPRelate)
library(GENESIS)
library(GWASTools)

# KING kinship - used pruned input
snpgdsBED2GDS(bed.fn = paste0(plink_prefix_pruned, ".bed"), 
              bim.fn = paste0(plink_prefix_pruned, ".bim"), 
              fam.fn = paste0(plink_prefix_pruned, ".fam"), 
              out.gdsfn = paste0(plink_prefix_pruned, ".gds"))

gds <- snpgdsOpen(paste0(plink_prefix_pruned, ".gds"))

krobust <- snpgdsIBDKING(gdsobj = gds,
                         type = "KING-robust")
snpgdsClose(gds)

KINGmat <- kingToMatrix(krobust)

saveRDS(KINGmat, paste0(plink_prefix, "_kingrobust_kinship_matrix.rds"))

# PC-AiR - use pruned input
    #TODO: why did Shane's code read in unpruned but then use snp.include?
    # For now, staying consistent with that here
kinobj <- divobj <- KINGmat
geno <- GdsGenotypeReader(paste0(plink_prefix, ".gds"))
genoData <- GenotypeData(geno)
pruned <- read.table(pruned_list)[,1]

mypcair <- pcair(genoData,
                 kinobj = kinobj,
                 divobj = divobj,
                 snp.include = pruned,
                 num.cores = 4,
                 kin.thresh = 2^(-11/2),
                 div.thresh = 2^(-11/2))
saveRDS(mypcair, paste0(plink_prefix, "_pcair.rds"))

# Make scree plot so can select number PCs to be used in PC-Relate
var_explained <- mypcair$values / sum(mypcair$values)
png(paste0(plink_prefix, "_pcair_scree_plot.png"), width = 400, height = 300)
plot(
  var_explained[1:10],
  type = "b",
  pch = 19,
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Elbow Plot"
)
dev.off()
