#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 2 of GENESIS.

# Purpose: “to get initial kinship estimates using KING-robust, which is robust
# to discrete population structure but not ancestry admixture. KING-robust will
# be able to identify close relatives (e.g. 1st and 2nd degree) reliably, but
# may identify spurious pairs or miss more distant pairs of relatives in the
# presence of admixture.”

# KING recommends not LD Pruning for their Tool but GENESIS has a different purpose.
# KING looks at all pairwise kinship but GENESIS wants to tease out local and
# global ancestry so that the kinship matrix is independent of global ancestry.
# So, we will LD prune before running KING.

#TODO: remove passing pruned file and pruned SNP list!

while getopts s:p:l: opt; do
   case "${opt}" in
      s) STUDY_PATH=${OPTARG};;  # path to study QCed PLINK1.9 prefix
      p) STUDY_PATH_PRUNED=${OPTARG};;  # path to pruned PLINK1.9 prefix
      l) PRUNED_LIST=${OPTARG};;  # path to file with list of pruned SNPs
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run KING and PC-AiR
   #TODO: need to add small sample arg to turn on/off
Rscript ${code_dir}/king_pcair.r ${STUDY_PATH_PRUNED} ${STUDY_PATH} ${PRUNED_LIST}
