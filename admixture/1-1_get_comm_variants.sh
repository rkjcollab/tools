#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1-1 of Admixture for ancestry estimation.

# Now we find variants that are common between TEDDY exome and 1000 Genomes & Merge.

set -e
set -u

while getopts r:o: opt; do
   case "${opt}" in
      r) REF_PATH=${OPTARG};;  # path to chr1 ref after 0-2, ${OUT}/qc/1000G.chr${chr}.qc.pruned.chrpos
      s) STUDY_PATH=${OPTARG};;  # path to study raw PLINK1.9 prefix
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/1000genomesPCA
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Find common variants between TEDDY and 1000 Genomes & Merge
#TODO: restart here and decide whether to match by RSID or chr:pos
# RAW="/Users/ridouxs/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/raw"


# Get/make paths
mkdir -p ${OUT}/study_1000G
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# Write out snplists for each chromosome
for chr in $(seq 1 22)
do
    ref_chr_path="${REF_PATH/chr1/chr$chr}"
	plink \
		--bfile "$ref_chr_path" \
		--write-snplist \
		--out ${TMP}/1000G.chr${chr}
done

# merge snplists
cat ${TMP}/1000G.chr*.snplist | sort | uniq > ${OUT}/study_1000G/1000G.snplist

# Write out snplist from study
#TODO: may have to split this into two commands!
plink2 \
	--bfile ${STUDY_PATH} \
    --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 999 \
	--write-snplist \
	--out ${OUT}/study_1000G/study

# Find common SNPs
sort ${OUT}/study_1000G/study.snplist > ${TMP}/study_snps.sorted.txt
sort ${OUT}/study_1000G/1000G.snplist > ${TMP}/1000G_snps.sorted.txt
comm -12 ${TMP}/study_snps.sorted.txt ${TMP}/1000G_snps.sorted.txt > ${OUT}/study_1000G/common_snps.txt

# Extract only common SNPs from each chromosome
for chr in $(seq 1 22)
do
	plink \
		--bfile ${OUT}/qc/1000G.chr${chr}.qc.pruned \
		--extract ${OUT}/teddy_1000G/common_snps.txt \
		--make-bed \
		--out ${TMP}/chr${chr}.common
done
# do again but exclude merge conflicted SNPs (4 of them)
for chr in $(seq 1 22)
do
	plink \
		--bfile ${OUT}/qc/1000G.chr${chr}.qc.pruned \
		--extract ${OUT}/teddy_1000G/common_snps.txt \
		--exclude ${OUT}/teddy_1000G/1000G_all_chrs-merge.missnp \
		--make-bed \
		--out ${TMP}/chr${chr}.common
done

# Merge chromosome data
cd ${TMP}
ls chr*.common.bed | sed 's/.bed$//' > ForMerge.list


plink \
	--merge-list ${TMP}/ForMerge.list \
	--make-bed \
	--write-snplist \
	--out ${OUT}/teddy_1000G/1000G_all_chrs

# Prepare TEDDY data with only common snps
plink \
	--bfile ${RAW}/teddy_r01/2025-03-05/OmicsDatasets/t1dexome_masked \
	--extract ${OUT}/teddy_1000G/1000G_all_chrs.snplist \
	--make-bed \
	--out ${OUT}/teddy_1000G/TEDDY.common
	