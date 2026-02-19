#!/usr/bin/env Rscript

# TO NOTE: this script uses genetic PCs and GMM method to estimate ancestral
# populations in reference and study inputs. Although arguments are described
# in context with Admixture inputs (since that step precedes this one), no
# Admixture outputs are used in this script to determine ancestral populations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plinkFile))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(pcaMethods))
suppressPackageStartupMessages(library(plotly))

parser <- ArgumentParser(
  description = paste0(
    "Rscript 2-3_get_eur.R run by run_admixture.Rmd. Note that the pipeline ",
    "is specific to using 1000G with K=5 as reference."))

# Reference arguments
parser$add_argument(
  "--ref", help=paste0(
    "Prefix of PLINK file input to Admixture for reference (required)"),
    required=TRUE)
parser$add_argument(
  "--demo-ref", help="Demographics file for reference (required)", required=TRUE)

# Study arguments
parser$add_argument(
  "--study", help=paste0(
    "Prefix of PLINK file input to Admixture and .Q file from Admixture for ",
    "study (required)"), required=TRUE)

args <- parser$parse_args()

# Setup ------------------------------------------------------------------------

# Set output prefix
out_path <- dirname(args$study)

# Load study files
q_study <- read.table(paste0(args$study, ".5.Q"))
fam_study <- read.table(paste0(args$study, ".fam"))

# Load reference files
q_ref <- read.table(paste0(args$ref, ".5.Q"))
demo_ref <- read.table(args$demo_ref, header = T)
fam_ref <- read.table(paste0(args$ref, ".fam"))
bed_ref <- readBED(paste0(args$ref, ".bed"))

# Load study files
q_study <- read.table(paste0(args$study, ".5.Q"))
fam_study <- read.table(paste0(args$study, ".fam"))
bed_study <- readBED(paste0(args$study, ".bed"))

# Update reference demo file
demo_ref <- demo_ref %>% 
  rename(IID = sample) %>% 
  select("IID","super_pop")

# Use clustering to define EUR -------------------------------------------------

# First, do PCA on reference
pca_ref <- pca(
  bed_ref,
  nPcs = 10,
  method = "svd",
  center = TRUE,
  scale = "none")

# Project study PCs based on reference
pca_study <- predict(pca_ref,
                     newdata = bed_study)

# Using 4 PCs based on this reference (this method is based off of it)
# https://linkinghub.elsevier.com/retrieve/pii/S1877050920324686
set.seed(123)
X <- scores(pca_ref)[,1:4]
gmm <- Mclust(data = X, G = 5)  # based on 1000G super-populations
# summary(gmm)

# Get summarized results of each reference ancestry group in each cluster
result_indv <- data.frame(
  classification = gmm$classification,
  pop = demo_ref$super_pop[match(names(gmm$classification), demo_ref$IID)]) %>%
  rownames_to_column("IID")
result <- as.data.frame(
  table(result_indv %>% dplyr::select(-IID)))

# Use maximum ancestry group in each cluster to assign population
map <- result %>%
  group_by(classification) %>%
  slice_max(Freq, with_ties = FALSE) %>%   # keep the classification with max counts
  ungroup() %>%
  dplyr::select(-Freq)

# Now,use study PCs for cluster/ancestry assignment
X <- pca_study$scores[,1:4]

study_ancestry <- predict(gmm, 
                          newdata = X)
study_ancestry$IID <- rownames(bed_study)

study_ancestry$pop <- map$pop[
  match(study_ancestry$classification, map$classification)]
# table(study_ancestry$pop)

# Add IDs
fam_study <- fam_study %>% 
  rename(FID = V1, IID = V2) %>% 
  select(FID, IID)
fam_study$pop <- study_ancestry$pop[match(fam_study$IID, study_ancestry$IID)]

# Plot population clusters -----------------------------------------------------

pcplot <- data.frame(
  PC1 = pca_study$scores[,1],
  PC2 = pca_study$scores[,2],
  PC3 = pca_study$scores[,3],
  group = study_ancestry$pop[match(rownames(pca_study$scores), study_ancestry$IID)])

# 2D plot
pops <- sort(unique(pcplot$group))
pop_colors <- RColorBrewer::brewer.pal(n = length(pops), name = "Set1")
names(pop_colors) <- pops
pcplot$color <- pop_colors[pcplot$group]

png(paste0(out_path, "/study_proj_ref_pca_2d.png"), width = 1200, height = 900)
plot(pcplot$PC1,
     pcplot$PC2,
     col = pcplot$color,
     xlab="PC1",
     ylab="PC2",
     main="Study projected onto reference PCA")
legend("topright", legend=names(pop_colors), col=pop_colors, pch=19)
dev.off()

# 3D plot
plot_study_anc_pc <- plot_ly(
  data = pcplot,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~group,
  colors = pop_colors,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    title = "Study projected onto reference PCA")

# Write out --------------------------------------------------------------------

# Update variables for clarity after write out
result_indv_exp <- result_indv %>%
  dplyr::rename(pop_1000g = pop) %>%
  dplyr::mutate(pop_gmm = map$pop[match(classification, map$classification)]) %>%
  dplyr::select(-classification)
fam_study_exp <- fam_study %>%
  dplyr::rename(pop_gmm = pop)

# Add admixture components
colnames(q_study) <- paste0("Ancestry", 1:ncol(q_study))
colnames(q_ref) <- paste0("Ancestry", 1:ncol(q_ref))
q_ref$IID <- fam_ref$V2
result_indv_exp <- result_indv_exp %>%
  dplyr::left_join(., q_ref, by = c("IID"))

q_study$IID <- fam_study$IID  # admixture preserves order
fam_study_exp <- fam_study_exp %>%
  dplyr::left_join(., q_study, by = c("IID"))

# Write out 1000G super-population labels and assigned clusters
write_tsv(
  result_indv_exp,
  paste0(out_path, "/1000G_all_pop.txt"))

# Write out assigned pop for all study individuals
write_tsv(
  fam_study_exp,
  paste0(out_path, "/study_all_pop.txt"))

# Write out list with EUR in study
fam_study_eur <- fam_study %>%
  dplyr::filter(pop == "EUR") %>%
  dplyr::select(FID, IID)
write_tsv(
  fam_study_eur,
  paste0(out_path, "/study_eur_list.txt"),
  col_names = F)

# Write out plots not written out above
htmlwidgets::saveWidget(
  plot_study_anc_pc, paste0(out_path, "/study_proj_ref_pca_3d.html"))
