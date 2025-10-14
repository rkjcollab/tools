#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1-1 of Admixture for ancestry estimation.

# Now we find variants that are common between TEDDY exome and 1000 Genomes & Merge.

set -e
set -u

while getopts r:s:o: opt; do
   case "${opt}" in
      r) REF_PATH=${OPTARG};;  # path to chr1 ref after 0-2, ${OUT}/qc/1000G.chr${chr}.qc.pruned
      s) STUDY_PATH=${OPTARG};;  # path to study raw PLINK1.9 prefix
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/1000genomesPCA
        # new out dir Immunogenetics_T1D/genetics/teddy_r01/admixture, make subfolder study_1000g
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Find common variants between TEDDY and 1000 Genomes & Merge, matching
# by "chr:pos" IDs
# RAW="/Users/ridouxs/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/raw"

# Get/make paths
mkdir -p ${OUT}/study_1000G
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# # Write out snplists for each reference chromosome
# for chr in $(seq 1 22)
# do
#     ref_chr_path="${REF_PATH/chr1/chr$chr}"
#     ref_chr_path_chrpos="${REF_PATH/chr1/chr$chr}.chrpos"
#     # First, write out version with chr:pos:ref:alt varIDs
#     plink2 \
#         --bfile "$ref_chr_path" \
#         --set-all-var-ids @:# --new-id-max-allele-len 999 \
#         --make-bed \
#         --out "$ref_chr_path_chrpos"
# 	plink \
# 		--bfile "$ref_chr_path_chrpos" \
#         --keep-allele-order \
# 		--write-snplist \
# 		--out ${TMP}/1000G.chr${chr}
# done

# # Merge reference snplists
# cat ${TMP}/1000G.chr*.snplist | sort | uniq > ${OUT}/study_1000G/1000G.snplist

# Write out snplist from study
#TODO: may have to split this into two commands!




#TODO: restart here - need different strategy... maybe merge by chr:pos:ref:alt,
# then for chr:pos for those that don't merge, then ???




plink2 \
	--bfile ${STUDY_PATH} \
    --set-all-var-ids @:# --new-id-max-allele-len 999 \
	--write-snplist  \
	--out ${OUT}/study_1000G/study

# Find common SNPs between reference and study
sort ${OUT}/study_1000G/study.snplist > ${TMP}/study_snps.sorted.txt
sort ${OUT}/study_1000G/1000G.snplist > ${TMP}/1000G_snps.sorted.txt
comm -12 ${TMP}/study_snps.sorted.txt ${TMP}/1000G_snps.sorted.txt > ${OUT}/study_1000G/common_snps.txt

### Extract only common reference SNPs
# Extract only common SNPs from each chromosome
for chr in $(seq 1 22)
do
    ref_chr_path_chrpos="${REF_PATH/chr1/chr$chr}.chrpos"
	plink \
		--bfile "$ref_chr_path_chrpos" \
        --keep-allele-order \
		--extract ${OUT}/study_1000G/common_snps.txt \
		--make-bed \
		--out ${TMP}/chr${chr}.common.tmp
done

# First merge attempt
ls ${TMP}/chr*.common.tmp.bed | sed 's/.bed$//' > ${TMP}/for_merge_tmp.list
plink \
	--merge-list ${TMP}/for_merge_tmp.list \
    --keep-allele-order \
	--make-bed \
	--write-snplist \
	--out ${TMP}/1000G_all_chrs_tmp

# If any exclude SNPs that failed merge
if [ -e "${TMP}/1000G_all_chrs_tmp-merge.missnp" ]
then
    for chr in $(seq 1 22)
    do
        plink \
            --bfile ${TMP}/chr${chr}.common.tmp \
            --keep-allele-order \
            --exclude ${TMP}/1000G_all_chrs_tmp-merge.missnp \
            --make-bed \
            --out ${TMP}/chr${chr}.common
    done

    # Second merge attempt
    ls ${TMP}/chr*.common.bed | sed 's/.bed$//' > ${TMP}/for_merge.list
    plink \
        --merge-list ${TMP}/for_merge.list \
        --keep-allele-order \
        --make-bed \
        --write-snplist \
        --out ${OUT}/study_1000G/1000G_all_chrs
else
    plink \
        --bfile ${TMP}/1000G_all_chrs_tmp \
        --keep-allele-order \
        --make-bed \
        --out ${OUT}/study_1000G/1000G_all_chrs
fi

### Extract only common study SNPs
# Prepare TEDDY data with only common snps
plink \
	--bfile ${STUDY_PATH} \
    --keep-allele-order \
	--extract ${OUT}/study_1000G/1000G_all_chrs.snplist \
	--make-bed \
	--out ${OUT}/study_1000G/study.common
	