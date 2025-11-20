#!/usr/bin/env Rscript

library(argparse)
suppressPackageStartupMessages(library(tidyverse))

parser <- ArgumentParser(
  description = paste0(
    "Rscript dedup_miss.R run by 1-1_get_comm_variants.sh. ",
    "Adapted from imputation repo."))

parser$add_argument(
  "-v", "--vmiss", help="PLINK2 .vmiss file path (required)", required=TRUE)
parser$add_argument(
  "-p", "--pvar", help="PLINK2 .pvar file path (required)", required=TRUE)

args <- parser$parse_args()

# Duplicate variants -----------------------------------------------------------

vmiss <- read.delim(args$vmiss)
pvar <- read.delim(args$pvar)

merge <- inner_join(
  pvar,
  vmiss,
  by = c("X.CHROM", "ID")) %>%
  dplyr::mutate(chrpos = paste0(X.CHROM, ":", POS, ":", REF, ":", ALT))

merge_dup_list <- merge %>%
  dplyr::filter(duplicated(chrpos))
merge_dups <- merge %>%
  dplyr::filter(chrpos %in% merge_dup_list$chrpos)

merge_dups_filt <- merge_dups %>%
  dplyr::group_by(chrpos) %>%
  dplyr::filter(
    any(F_MISS > 0),
    length(unique(F_MISS)) > 1)

merge_dups_keep <- merge_dups_filt %>%
  slice_min(order_by = F_MISS, n = 1) %>%
  ungroup()

merge_dups_remove <- merge_dups_filt %>%
  dplyr::filter(!ID %in% merge_dups_keep$ID)
merge_dups_remove_ct <- n_distinct(merge_dups_remove$ID)

# Write out --------------------------------------------------------------------

out_path <- paste0(gsub("\\.vmiss", "", args$vmiss), "_rm.txt")
write.table(
  merge_dups_remove$ID,
  file = out_path,
  sep="\t",
  quote=F,
  row.names=F,
  col.names=F)
