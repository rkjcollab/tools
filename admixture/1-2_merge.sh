#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1-2 of Admixture for ancestry estimation.

# Now let's merge (and flip variants if necessary and remove missing vars).

#TODO: current version follows SER's pipeline. Could also modify to more closely
# follow Michelle Daya's pipeline, modified by SDS:
    # https://github.com/sdslack/merge_genetic_data/blob/master/merge_data_sets.sh

while getopts o: opt; do
   case "${opt}" in
      o) OUT=${OPTARG};;  # path to out dir
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
mkdir -p ${OUT}/study_1000G
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# Write out list of flipped SNPs, merge-mode 6 = no merge, report all mismatching calls
plink \
   --bfile ${OUT}/study_1000G/study.common \
   --keep-allele-order \
   --bmerge ${OUT}/study_1000G/1000G_all_chrs \
   --merge-mode 6 \
   --out ${OUT}/study_1000G/study_mergecheck
      
# If any, flip SNPs
if [ -e "${OUT}/study_1000G/study_mergecheck.missnp" ]
then
    plink \
        --bfile ${OUT}/study_1000G/study.common \
        --keep-allele-order \
        --flip ${OUT}/study_1000G/study_mergecheck.missnp \
        --make-bed \
        --out ${OUT}/study_1000G/study.common.flipped
        
    # Recheck with flipped SNPs, merge-mode 6 = no merge, report all mismatching calls
    plink \
        --bfile ${OUT}/study_1000G/study.common.flipped \
        --keep-allele-order \
        --bmerge ${OUT}/study_1000G/1000G_all_chrs \
        --merge-mode 6 \
        --out ${OUT}/study_1000G/study_mergecheck2

    # If any remaining allele mismatches, exclude them
    if [ -e "${OUT}/study_1000G/study_mergecheck2.missnp" ]
    then
        plink \
            --bfile ${OUT}/study_1000G/study.common.flipped \
            --keep-allele-order \
            --exclude ${OUT}/study_1000G/study_mergecheck2.missnp \
            --make-bed \
            --out ${OUT}/study_1000G/study.common.cleaned
    else
        plink \
            --bfile ${OUT}/study_1000G/study.common.flipped \
            --keep-allele-order \
            --make-bed \
            --out ${OUT}/study_1000G/study.common.cleaned
    fi
else
    plink \
        --bfile ${OUT}/study_1000G/study.common \
        --keep-allele-order \
        --make-bed \
        --out ${OUT}/study_1000G/study.common.cleaned
fi

# Merge
awk '{print $2}' ${OUT}/study_1000G/study.common.cleaned.bim > ${TMP}/common_snp_list.txt
plink \
   --bfile ${OUT}/study_1000G/study.common.cleaned \
   --keep-allele-order \
   --bmerge ${OUT}/study_1000G/1000G_all_chrs \
   --make-bed \
   --extract ${TMP}/common_snp_list.txt \
   --out ${OUT}/study_1000G/study_1000G_merged

# Get sample IDs for study and 1000G
awk '{print $1, $2}' ${OUT}/study_1000G/study.common.cleaned.fam > ${OUT}/study_1000G/study_samples.txt
awk '{print $1, $2}' ${OUT}/study_1000G/1000G_all_chrs.fam > ${OUT}/study_1000G/1000G_samples.txt
