#!/bin/bash

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 0-2 of Admixture for ancestry estimation.
# NOTE: this script is step 0 as it only needs to be run once and then its results
# can be used across projects/studies.

# We will split multiallelics and then subset to just biallelic SNPs and remove
# duplicates (exact meaning only removes duplicates in same position with same
# alt allele so there may be snps at the same position with differing alt alleles).
# Now we convert to plink binary from vcf and maf filter at 0.1, LD prune at 5 50 0.2,
# and restrict to good regions in the strict file set provided by 1000 Genomes.

set -e
set -u

while getopts r:m:o: opt; do
   case "${opt}" in
      r) REF_PATH=${OPTARG};;  # path to chr1 ref
      m) MASK_PATH=${OPTARG};;  # path to 1000G mask
      o) OUT=${OPTARG};;  # path to out dir, Immunogenetics_T1D/genetics/teddy_r01/1000genomesPCA
        # new out dir Immunogenetics_T1D/genetics/teddy_r01/admixture, make subfolders normalize_1000g, qc_1000g
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Get/make paths
mkdir -p ${OUT}/normalize
mkdir -p ${OUT}/qc
TMP=$(mktemp -d -p "$OUT")
trap 'rm -rf "$TMP"' EXIT

# Normalize VCF
echo "Step1: Normalize VCF"
echo "------------------------------"

# split multiallelics then subset to just biallelic SNPs then rm duplicates
for chr in $(seq 1 22)
do
  echo "Processing chr${chr}..."
  ref_chr_path="${REF_PATH/chr1/chr$chr}"

  bcftools norm \
      --multiallelics -any "$ref_chr_path" | \
      bcftools view \
          --types snps \
          --min-alleles 2 \
          --max-alleles 2 | \
          bcftools norm \
              --remove-duplicates \
              --output-type z \
              --output ${OUT}/normalize/1000G.chr${chr}.bcfnorm.snps.vcf.gz
  bcftools index ${OUT}/normalize/1000G.chr${chr}.bcfnorm.snps.vcf.gz
done

# Now we convert to plink binary from vcf and maf filter at 0.1, LD prune at 5 50 0.2,
# and restrict to good regions in the strict file set provided by 1000 Genomes.
echo "Step2: Convert VCFs to PLINK Binary"
echo "-----------------------------------"

for chr in $(seq 1 22)
do
	echo "Processing chr${chr}..."

	# Step 1: Convert and filter
  plink \
    --vcf ${OUT}/normalize/1000G.chr${chr}.bcfnorm.snps.vcf.gz \
    --keep-allele-order \
    --double-id \
    --make-bed \
    --snps-only \
    --maf 0.1 \
    --out ${TMP}/1000G.chr${chr}.qc

  # Step 2: LD prune & extract to well recognized snps
  plink \
    --bfile ${TMP}/1000G.chr${chr}.qc \
    --keep-allele-order \
    --make-set "$MASK_PATH" \
    --write-set \
    --out ${TMP}/strict_set
  plink \
    --bfile ${TMP}/1000G.chr${chr}.qc \
    --keep-allele-order \
    --extract ${TMP}/strict_set.set \
    --indep-pairwise 50 5 0.2 \
    --out ${TMP}/1000G.chr${chr}.prune

  # Step 3: Extract pruned SNPs
  plink \
    --bfile ${TMP}/1000G.chr${chr}.qc \
    --keep-allele-order \
    --extract ${TMP}/1000G.chr${chr}.prune.prune.in \
    --make-bed \
    --out ${OUT}/qc/1000G.chr${chr}.qc.pruned
done
