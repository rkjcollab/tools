#!/usr/bin/env Rscript

# TO NOTE: this script uses ancestral population groups from GMM in step 2-3
# and Admixture components for visualization of ancestral population separation.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

parser <- ArgumentParser(
  description = paste0(
    "Rscript 2-4_plot_pops.R run by run_admixture.Rmd. Note that the pipeline ",
    "is specific to using 1000G with K=5 as reference."))

# Study arguments
parser$add_argument(
  "--study", help=paste0(
    "Prefix of PLINK file input to Admixture and .Q file from Admixture for ",
    "study (required)"), required=TRUE)

args <- parser$parse_args()

# Setup ------------------------------------------------------------------------

# Set output prefix
out_path <- dirname(args$study)

# Load GMM results from 2-3
gmm_study <- read.table(paste0(out_path, "/study_all_pop.txt"), header = T)
gmm_ref <- read.table(paste0(out_path, "/1000G_all_pop.txt"), header = T)

# Plot reference admixture by 1000G pops -------------------------------------

# This section shows the Admixture ancestry components for each reference
# super-population group. No GMM results included.

# Factor by reference pop, order by decreasing Ancestry1
df_admix <- gmm_ref %>%
  arrange(pop_1000g) %>%
  mutate(pop_1000g = factor(pop_1000g)) %>%
  arrange(pop_1000g, desc(Ancestry1)) %>%
  dplyr::mutate(id = row_number()) %>%
  ungroup()

# Then pivot longer
df_admix_long <- df_admix %>%
  pivot_longer(
    cols = starts_with("Ancestry"),
    names_to = "AncestryComponent",
    values_to = "Proportion")

# Plot stacked barplot with population labels
anc_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:5]
plot_ref_admix_bar <-
  ggplot(df_admix_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(
    values = anc_colors,
    labels = function(x) sub("^Ancestry", "", x)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population Label (1000G)",
       y = "Proportion Ancestry",
       fill = "Ancestry\nComponent\n(Admixture)") +
  facet_grid(~ pop_1000g, scales = "free_x", space = "free_x", switch = "x")

# Plot reference admixture by GMM groups ---------------------------------------

# This section shows the Admixture ancestry components in the reference
# for each GMM group.

# Factor by reference pop, order by decreasing Ancestry1
df_gmm <- gmm_ref %>%
  arrange(pop_gmm) %>%
  mutate(pop_gmm = factor(pop_gmm)) %>%
  arrange(pop_gmm, desc(Ancestry1)) %>%
  dplyr::mutate(id = row_number()) %>%
  ungroup()

# Then pivot longer
df_gmm_long <- df_gmm  %>%
  pivot_longer(
    cols = starts_with("Ancestry"),
    names_to = "AncestryComponent",
    values_to = "Proportion")

# Plot stacked barplot with population labels
anc_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:5]
plot_ref_gmm_bar <-
  ggplot(df_gmm_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(
    values = anc_colors,
    labels = function(x) sub("^Ancestry", "", x)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population Label (GMM)",
       y = "Proportion Ancestry",
       fill = "Ancestry\nComponent\n(Admixture)") +
  facet_grid(~ pop_gmm, scales = "free_x", space = "free_x", switch = "x")


# Plot study results by GMM pops -----------------------------------------------

# This section includes both admixture-estimated components and GMM selected
# ancestry groups for the study.

# Factor by reference pop, order by decreasing Ancestry1
df_gmm_study <- gmm_study %>%
  arrange(pop_gmm) %>%
  mutate(pop_gmm = factor(pop_gmm)) %>%
  arrange(pop_gmm, desc(Ancestry1)) %>%
  dplyr::mutate(id = row_number()) %>%
  ungroup()

# Then pivot longer
df_gmm_study_long <- df_gmm_study %>%
  pivot_longer(
    cols = starts_with("Ancestry"),
    names_to = "AncestryComponent",
    values_to = "Proportion")

# Plot stacked barplot with population labels
plot_study_gmm_bar <- ggplot(
  df_gmm_study_long, aes(x = id, y = Proportion, fill = AncestryComponent)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(
    values = anc_colors,
    labels = function(x) sub("^Ancestry", "", x)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = 10)) +
  labs(x = "Population Label (GMM)",
       y = "Proportion Ancestry",
       fill = "Ancestry\nComponent\n(Admixture)") +
  facet_grid(~ pop_gmm, scales = "free_x", space = "free_x", switch = "x")

plot_study_gmm_bar_equal <-
  plot_study_gmm_bar + facet_grid(~ pop_gmm, scales = "free_x", switch = "x")

# Write out --------------------------------------------------------------------

# Write out plots not written out above
ggsave(paste0(out_path, "/ref_admixture_bar.png"),
       plot_ref_admix_bar, units = "in", width = 7, height = 4)
ggsave(paste0(out_path, "/ref_gmm_bar.png"),
       plot_ref_gmm_bar, units = "in", width = 7, height = 4)

ggsave(paste0(out_path, "/study_gmm_bar.png"),
       plot_study_gmm_bar, units = "in", width = 7, height = 4)
ggsave(paste0(out_path, "/study_gmm_bar_even_scale.png"),
       plot_study_gmm_bar_equal, units = "in", width = 7, height = 4)
