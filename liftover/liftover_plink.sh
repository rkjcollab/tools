#!/bin/bash

# Set arguments
if [ "$#" -lt  "3" ]
then
   echo "Usage: ${0##*/} <input_file> <out_dir> <target_build: hg19|hg38>"
   exit
fi

# Liftover between hg19 and hg38 with for PLINK1.9 files.

# Script requires PLINK2 and CrossMap v0.7.0 or higher. Can install based on
# directions here: https://crossmap.readthedocs.io/en/latest/#installation
# Chain files and helper code files should be in same directory as this script.
# Currently needs to be run from liftover code directory.

# STILL TODO:
# make conda environment and recipe for this?
# update to handle PLINK2 input as well

input_file=$1  # path to and name PLINK1.9 .bim file
out_dir=$2
target_build=$3

if [ "$target_build" = "hg19" ]
then
   chain_file="hg38ToHg19.over.chain"
elif [ "$target_build" = "hg38" ]
then
   chain_file="hg19ToHg38.over.chain"
else
   echo "Error: target_build must be hg19 or hg38"
   exit 1
fi

# Get input_file without path and without extension
input_file_no_path=$(basename "$input_file")
input_file_name=${input_file_no_path%.*}
input_file_no_ext=${input_file%.*}

# Create bed file to crossover from hg38 to hg19
cat $input_file | cut -f1 > ${out_dir}/tmp_c1.txt
cat $input_file | cut -f4 > ${out_dir}/tmp_c2.txt
cat $input_file | cut -f2 > ${out_dir}/tmp_c3.txt
paste ${out_dir}/tmp_c1.txt \
    ${out_dir}/tmp_c2.txt \
    ${out_dir}/tmp_c2.txt \
    ${out_dir}/tmp_c3.txt \
    > ${out_dir}/tmp_in.bed

# Do crossover
CrossMap bed $chain_file \
   ${out_dir}/tmp_in.bed  \
   ${out_dir}/tmp_out.bed

# Extract only those SNPs that were successfully cross-overed
cut -f4 ${out_dir}/tmp_out.bed > ${out_dir}/tmp_snp_keep.txt
plink2 --bfile $input_file_no_ext \
    --extract ${out_dir}/tmp_snp_keep.txt \
    --make-bed --out ${out_dir}/tmp_gwas

# Update bim file positions
Rscript --vanilla update_pos.R \
    ${out_dir}/tmp_out.bed ${out_dir}/tmp_gwas.bim

# Set all varids to chr:pos:ref:alt and sort
plink2 --bfile ${out_dir}/tmp_gwas \
    --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 500 \
    --sort-vars \
    --make-pgen --out ${out_dir}/tmp_gwas_id_sort

# Write out final 1.9 format file
plink2 --pfile ${out_dir}/tmp_gwas_id_sort \
    --make-bed --out ${out_dir}/${input_file_name}_${target_build}

# Report SNP counts
orig_snp_nr=`wc -l ${input_file}`
crossover_snp_nr=`wc -l ${out_dir}/${input_file_name}_${target_build}.bim`

echo "Original SNP nr: $orig_snp_nr"
echo "Crossovered SNP nr: $crossover_snp_nr"
echo "Original SNP nr: $orig_snp_nr" > ${out_dir}/liftover_plink_log.txt
echo "Crossovered SNP nr: $crossover_snp_nr" >> ${out_dir}/liftover_plink_log.txt

# Cleanup
rm ${out_dir}/tmp_*
