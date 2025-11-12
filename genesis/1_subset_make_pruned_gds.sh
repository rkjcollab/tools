#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1 of GENESIS.

# First we subset to the EUR or NHW group that we determined from the ADMIXTURE analysis
# SDS adding additional list subset option here
   #TODO: revisit and find better way to do this!

while getopts s:e:i:o: opt; do
   case "${opt}" in
      s) STUDY_PATH=${OPTARG};;  # path to study QCed PLINK1.9 prefix
      e) EUR_ID_LIST=${OPTARG};;  # path to EUR ID list
      i) OTHER_ID_LIST=${OPTARG};;  # optional path to other ID filtering list
      o) OUT=${OPTARG};;  # path to out dir
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
mkdir -p ${OUT}
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

echo "Filtering to EUR only."
plink --bfile ${STUDY_PATH} \
      --keep ${EUR_ID_LIST} \
      --make-bed \
      --nonfounders \
      --keep-allele-order \
      --out ${TMP}/tmp_study_nhw

if [ -n "${OTHER_ID_LIST}" ]; then
   echo "Also filtering to other ID list only."
   plink --bfile ${TMP}/tmp_study_nhw \
      --keep ${OTHER_ID_LIST} \
      --make-bed \
      --nonfounders \
      --keep-allele-order \
      --out ${OUT}/study_nhw
else
   plink --bfile ${TMP}/tmp_study_nhw \
      --make-bed \
      --nonfounders \
      --keep-allele-order \
      --out ${OUT}/study_nhw
fi

# Convert PLINK1.9 to GDS
Rscript ${code_dir}/convert_plink_gds.r ${OUT}/study_nhw

# LD prune GDS
Rscript ${code_dir}/ld_prune.r ${OUT}/study_nhw

plink \
   --bfile ${OUT}/study_nhw \
   --extract ${OUT}/study_nhw_prune.in \
   --make-bed \
   --nonfounders \
   --keep-allele-order \
   --out ${OUT}/study_nhw_pruned