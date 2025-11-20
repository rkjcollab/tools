
#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 2-2 of Admixture for ancestry estimation.

while getopts o:k: opt; do
   case "${opt}" in
      o) OUT=${OPTARG};;  # path to out dir
      k) K=${OPTARG};;  # K value from reference population (script 2-1)
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Set wd so admixture output goes to right spot
cd ${OUT}/ancestry_estimation

# Run unsupervised ADMIXTURE with K selected from 2-1
# TO NOTE: parallel hardcoded as 6 to work for most systems
admixture -j6 ${OUT}/ancestry_estimation/1000G_ref.bed ${K}

# Use learned allele frequencies as (fixed) input to next step
cp ${OUT}/ancestry_estimation/1000G_ref.${K}.P ${OUT}/ancestry_estimation/study.${K}.P.in

# Run projection ADMIXTURE with K from 2-1
# TO NOTE: parallel hardcoded as 6 to work for most systems
admixture -j6 -P ${OUT}/ancestry_estimation/study.bed ${K}
