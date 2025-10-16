#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 2-1 of Admixture for ancestry estimation.
# NOTE: this script can be skipped as it only needs to be run once and
# then its results can be used across projects/studies.

while getopts o: opt; do
   case "${opt}" in
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/ancestry_estimation
        # new out dir Immunogenetics_T1D/genetics/teddy_r01/admixture, make subfolder ancestry_estimation
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Run unsupervised ADMIXTURE with 5 fold cv
for K in 1 2 3 4 5; \
   do admixture --cv=5 ${OUT}/ancestry_estimation/1000G_ref.bed $K | tee log${K}.out; \
done

grep -h CV ${OUT}/log*.out > ${OUT}/CV_errors.txt


