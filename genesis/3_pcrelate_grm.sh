#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 3 of GENESIS.

while getopts s:l:n: opt; do
   case "${opt}" in
      s) STUDY_PATH=${OPTARG};;  # path to study QCed PLINK1.9 prefix
      l) PRUNED_LIST=${OPTARG};;  # path to file with list of pruned SNPs
      n) NUM_PCS=${OPTARG};;  # number PCs chosen from scree plot
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run PCRelate and make GRM, set to number of PCs identified by PCAiR
   #TODO: need to add small sample arg to turn on/off
Rscript ${code_dir}/pcrelate_grm.r \
    ${STUDY_PATH} ${PRUNED_LIST} ${NUM_PCS}