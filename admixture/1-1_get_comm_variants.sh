#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1-1 of Admixture for ancestry estimation.

# Now we find variants that are common between TEDDY exome and 1000 Genomes & Merge.

while getopts r:s:c:o: opt; do
   case "${opt}" in
      r) REF_PATH=${OPTARG};;  # path to chr1 ref after 0-2, ${OUT}/qc/1000G.chr${chr}.qc.pruned
      s) STUDY_PATH=${OPTARG};;  # path to study raw PLINK1.9 prefix
      c) CODE_DIR=${OPTARG};;  # path to tools/admixture code
      o) OUT=${OPTARG};;  # path to out dir
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Find common variants between TEDDY and 1000 Genomes & Merge, matching
# by "chr:pos" IDs

# Get/make paths
rm -rf ${OUT}/study_1000G  # remove previous results if already present
mkdir -p ${OUT}/study_1000G
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# Write out snplists for each reference chromosome
for chr in $(seq 1 22)
do
    ref_chr_path="${REF_PATH/chr1/chr$chr}"
    ref_chr_path_chrpos="${REF_PATH/chr1/chr$chr}.chrpos"
    # First, write out version with chr:pos:ref:alt varIDs
    plink2 \
        --bfile "$ref_chr_path" \
        --set-all-var-ids @:# --new-id-max-allele-len 999 \
        --make-bed \
        --out "$ref_chr_path_chrpos"
	plink \
		--bfile "$ref_chr_path_chrpos" \
        --keep-allele-order \
		--write-snplist \
		--out ${TMP}/1000G.chr${chr}
done

# Merge reference snplists
cat ${TMP}/1000G.chr*.snplist | sort | uniq > ${OUT}/study_1000G/1000G.snplist

# Remove multiallelic, keep only autosomes in analysis
plink2 --bfile ${STUDY_PATH} \
    --max-alleles 2 \
    --autosome \
    --make-bed --out ${TMP}/tmp_study_autosome

# First, if duplicate IDs have different missingness, remove the SNP with
# more missingness.
plink2 --bfile ${TMP}/tmp_study_autosome \
   --missing --freq \
   --make-pgen \
   --out ${TMP}/tmp_study_dedup
Rscript ${CODE_DIR}/dedup_miss.R \
   -v ${TMP}/tmp_study_dedup.vmiss \
   -p ${TMP}/tmp_study_dedup.pvar
plink2 --bfile ${TMP}/tmp_study_autosome \
   --make-bed \
   --exclude ${TMP}/tmp_study_dedup_rm.txt \
   --out ${OUT}/study_1000G/study.nodupl

# Then, update all IDs to chr:pos:ref:alt and remove the remaining duplicates
# by arbitrarily keeping first
plink2 \
	--bfile ${OUT}/study_1000G/study.nodupl \
    --set-all-var-ids @:# --new-id-max-allele-len 999 \
    --rm-dup force-first \
    --write-snplist \
    --make-bed \
	--out ${OUT}/study_1000G/study.nodupl.chrpos.snps

# Find common SNPs between reference and study
sort ${OUT}/study_1000G/study.nodupl.chrpos.snps.snplist > ${TMP}/study_snps.sorted.txt
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

# If any, exclude SNPs that failed merge
#TODO: in SDS testing, these seem to be multi-allelic and in reference dataset.
# Do we want to revisit this decision?
# TO NOTE: as is, can't run script to fail on error because would get stuck here
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
# Prepare study data with only common snps
plink \
	--bfile ${OUT}/study_1000G/study.nodupl.chrpos.snps \
    --keep-allele-order \
	--extract ${OUT}/study_1000G/1000G_all_chrs.snplist \
	--make-bed \
	--out ${OUT}/study_1000G/study.common
	