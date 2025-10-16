#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1-3 of Admixture for ancestry estimation.

while getopts o: opt; do
   case "${opt}" in
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/ancestry_estimation
        # new out dir Immunogenetics_T1D/genetics/teddy_r01/admixture, make subfolder ancestry_estimation
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
mkdir -p ${OUT}/ancestry_estimation
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# LD prune
plink \
    --bfile ${OUT}/study_1000G/study_1000G_merged \
    --keep-allele-order \
    --indep-pairwise 50 5 0.2 \
    --out ${OUT}/ancestry_estimation/study_1000G_ld_prune

plink \
    --bfile ${OUT}/study_1000G/study_1000G_merged \
    --extract ${OUT}/ancestry_estimation/study_1000G_ld_prune.prune.in \
    --maf 0.1 \
    --make-bed \
    --out ${OUT}/ancestry_estimation/study_1000G_pruned_maf0.1

# Make 1000G ref file
plink \
    --bfile ${OUT}/ancestry_estimation/study_1000G_pruned_maf0.1 \
    --keep ${OUT}/study_1000G/1000G_samples.txt \
    --make-bed \
    --nonfounders \
    --keep-allele-order \
    --out ${OUT}/ancestry_estimation/1000G_ref
      
# Make study file
plink \
    --bfile ${OUT}/ancestry_estimation/study_1000G_pruned_maf0.1 \
    --keep ${OUT}/study_1000G/study_samples.txt \
    --make-bed \
    --nonfounders \
    --keep-allele-order \
    --out ${OUT}/ancestry_estimation/study

#TODO: may need to check completely missing SNPs here or above!
