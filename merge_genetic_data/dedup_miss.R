#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(
  description = "Choose duplicate variants to remove by missingness.")

parser$add_argument(
  "-v", "--vmiss", help="PLINK2 .vmiss file path (required)", required=TRUE)
parser$add_argument(
  "-p", "--pvar", help="PLINK2 .pvar file path (required)", required=TRUE)

args <- parser$parse_args()

# Duplicate variants -----------------------------------------------------------

vmiss <- read.delim(args$vmiss)
pvar <- read.delim(args$pvar)

pvar$row_id <- seq_len(nrow(pvar))
merge <- merge(pvar, vmiss, by = c("X.CHROM", "ID"), sort = FALSE)
merge <- merge[order(merge$row_id), ]
merge$row_id <- NULL

merge$chrpos <- paste0(merge$X.CHROM, ":", merge$POS, ":", merge$REF, ":", merge$ALT)
merge_dup_list <- merge$chrpos[duplicated(merge$chrpos)]
merge_dups <- merge[merge$chrpos %in% merge_dup_list, ]

if (nrow(merge_dups) > 0) {
  merge_dups_filt_list <-
    lapply(split(merge_dups, merge_dups$chrpos), function(df) {
      if (any(df$F_MISS > 0) && length(unique(df$F_MISS)) > 1) {
        df
      }
    })
  merge_dups_filt_list <- Filter(Negate(is.null), merge_dups_filt_list)

  if (length(merge_dups_filt_list) > 0) {
    merge_dups_filt <- do.call(rbind, merge_dups_filt_list)

    merge_dups_keep <- do.call(
      rbind,
      lapply(split(merge_dups_filt, merge_dups_filt$chrpos), function(df) {
        df[which.min(df$F_MISS), , drop = FALSE]
      })
    )

    merge_dups_remove <-
      merge_dups_filt[!merge_dups_filt$ID %in% merge_dups_keep$ID, ]
    merge_dups_remove_ct <- length(unique(merge_dups_remove$ID))
    merge_dups_remove_exp <- merge_dups_remove$ID
  } else {
    merge_dups_remove_ct <- 0
    merge_dups_remove_exp <- character()
  }
} else {
  merge_dups_remove_ct <- 0
  merge_dups_remove_exp <- character()
}

# Write out --------------------------------------------------------------------

out_path <- paste0(gsub("\\.vmiss", "", args$vmiss), "_rm.txt")
write.table(
  merge_dups_remove_exp,
  file = out_path,
  sep="\t",
  quote=F,
  row.names=F,
  col.names=F)
