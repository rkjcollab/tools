#!/usr/bin/env Rscript

library(argparse)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(plinkFile))
suppressPackageStartupMessages(library(pcaMethods))
suppressPackageStartupMessages(library(plotly))

parser <- ArgumentParser(
  description = paste0(
    "Rscript cluster_eur.R run by run_admixture.Rmd. Note that the pipeline ",
    "is specific to using 1000G with K=5 as reference."))

# Reference arguments
parser$add_argument(
  "--ref", help=paste0(
    "Prefix of PLINK file input to Admixture and .Q file from Admixture for ",
    "reference (required)"), required=TRUE)
parser$add_argument(
  "--demo-ref", help="Demographics file for reference (required)", required=TRUE)

# Study arguments
parser$add_argument(
  "--study", help=paste0(
    "Prefix of PLINK file input to Admixture and .Q file from Admixture for ",
    "study (required)"), required=TRUE)

args <- parser$parse_args()

# Setup ------------------------------------------------------------------------

# Load reference files
q_ref <- read.table(paste0(args$ref, ".5.Q"))
demo_ref <- read.table(args$demo_ref, header = T)
fam_ref <- read.table(paste0(args$ref, ".fam"))
bed_ref <- readBED(paste0(args$ref, ".bed"))

# Load study files
q_study <- read.table(paste0(args$study, ".5.Q"))
fam_study <- read.table(paste0(args$study, ".fam"))
bed_study <- readBED(paste0(args$study, ".bed"))

# # TEMP DELETE
# q_ref <- read.table("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/1000G_ref.5.Q")
# demo_ref <- read.table("Immunogenetics_T1D/raw/genetic_maps/1000genomes-phase3/integrated_call_samples_v3.20130502.ALL.panel", header = T)
# fam_ref <- read.table("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/1000G_ref.fam")
# bed_ref <- readBED("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/1000G_ref.bed")
# 
# # Load study files
# q_study <- read.table("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/study.5.Q")
# fam_study <- read.table("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/study.fam")
# bed_study <- readBED("DAISY/genetics/daisy_ask_genetics/admixture/ancestry_estimation/study.bed")
# #####

# Update column names
colnames(q_study) <- paste0("Ancestry", 1:ncol(q_study))
colnames(q_ref) <- paste0("Ancestry", 1:ncol(q_ref))

# Update study and reference .fam files
fam_ref <- fam_ref %>% 
  rename(FID = V1,
         IID = V2) %>% 
  left_join(demo_ref, by = c("IID"="sample")) %>% 
  select("IID","super_pop")
q_ref$IID <- fam_ref$IID
q_ref$pop <- fam_ref$super_pop

# Output prefix for plots that have to be written out as made
out_path <- dirname(args$study)

# Plot reference admixture -----------------------------------------------------

# Factor by reference pop, order by decreasing Ancestry1
q_ref <- q_ref %>%
  arrange(pop) %>%
  mutate(pop = factor(pop)) %>%
  arrange(pop, desc(Ancestry1)) %>%
  dplyr::mutate(id = row_number()) %>%
  ungroup()

# Then pivot longer
q_ref_long <- q_ref %>%
  pivot_longer(
    cols = starts_with("Ancestry"),
    names_to = "AncestryComponent",
    values_to = "Proportion")

# Plot stacked barplot with population labels
anc_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:5]
plot_ref_anc_bar <-
  ggplot(q_ref_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = anc_colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population (from Reference)",
       y = "Ancestry (by Admixture)",
       title = "Reference Populations") +
  facet_grid(~ pop, scales = "free_x", space = "free_x")

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
gmm <- Mclust(data = X, G = 5)  # based on 1000G K=5
# summary(gmm)

# Get summarized results of each reference ancestry group in each cluster
result <- as.data.frame(
  table(data.frame(classification = gmm$classification,
                   pop = q_ref$pop[match(names(gmm$classification), q_ref$IID)])))

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

# Plot study results -----------------------------------------------------------

# Add IDs
fam_study <- fam_study %>% 
  rename(FID = V1, IID = V2) %>% 
  select(FID, IID)
q_study$IID <- fam_study$IID  # admixture preserves order
q_study$FID <- fam_study$FID
q_study$pop <- study_ancestry$pop[match(q_study$IID, study_ancestry$IID)]

# Factor by reference pop, order by decreasing Ancestry1
q_study <- q_study %>%
  arrange(pop) %>%
  mutate(pop = factor(pop)) %>%
  arrange(pop, desc(Ancestry1)) %>%
  dplyr::mutate(id = row_number()) %>%
  ungroup()

# Then pivot longer
q_study_long <- q_study %>%
  pivot_longer(
    cols = starts_with("Ancestry"),
    names_to = "AncestryComponent",
    values_to = "Proportion")

# Plot stacked barplot with population labels
plot_study_anc_bar <- ggplot(q_study_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = anc_colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population (by Clustering)",
       y = "Ancestry (by Admixture)",
       title = "Study Populations") +
  facet_grid(~ pop, scales = "free_x", space = "free_x")

plot_study_anc_bar_v2 <- ggplot(q_study_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = anc_colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population (by Clustering)",
       y = "Ancestry (by Admixture)",
       title = "Study Populations") +
  facet_grid(~ pop, scales = "free_x")


# Plot population clusters
pcplot <- data.frame(
  PC1 = pca_study$scores[,1],
  PC2 = pca_study$scores[,2],
  PC3 = pca_study$scores[,3],
  group = study_ancestry$pop[match(rownames(pca_study$scores), study_ancestry$IID)])

# 2D plot
pop_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set1")[1:5]
names(pop_colors) <- unique(q_study$pop)
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

# Write out plots not written out above
ggsave(paste0(out_path, "/ref_admixture_bar.png"),
       plot_ref_anc_bar, units = "in", width = 20, height = 10)
ggsave(paste0(out_path, "/study_admixture_bar.png"),
       plot_study_anc_bar, units = "in", width = 20, height = 10)
ggsave(paste0(out_path, "/study_admixture_bar_even_scale.png"),
       plot_study_anc_bar_v2, units = "in", width = 20, height = 10)
htmlwidgets::saveWidget(
  plot_study_anc_pc, paste0(out_path, "/study_proj_ref_pca_3d.html"))

# Write out all populations results
write_tsv(
  q_study %>%
    dplyr::select(-id) %>%
    dplyr::relocate(FID, IID),
  paste0(out_path, "/study_pops.txt"))

# Write out list with > 80% EUR in DAISY
q_study_eur <- q_study %>%
  dplyr::filter(pop == "EUR") %>%
  dplyr::select(FID, IID)
write_tsv(
  q_study_eur,
  paste0(out_path, "/study_eur_list.txt"),
  col_names = F)
