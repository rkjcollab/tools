
#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 2-2 of Admixture for ancestry estimation.

while getopts o:k: opt; do
   case "${opt}" in
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/ancestry_estimation
      k) K=${OPTARG};;
        # new out dir Immunogenetics_T1D/genetics/teddy_r01/admixture, make subfolder ancestry_estimation
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Set wd so admixture output goes to right spot
cd ${OUT}/ancestry_estimation

# # Run unsupervised ADMIXTURE with K selected from 2-1
# admixture ${OUT}/ancestry_estimation/1000G_ref.bed ${K}

# Use learned allele frequencies as (fixed) input to next step
cp ${OUT}/ancestry_estimation/1000G_ref.${K}.P ${OUT}/ancestry_estimation/study.${K}.P.in

# Run projection ADMIXTURE with K=5
admixture -P ${OUT}/ancestry_estimation/study.bed ${K}
