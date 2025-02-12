#!/bin/bash

# Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <input_format> <input_file> <out_dir>"
   exit
fi

# Liftover from hg38 to hg19 with for PLINK1.9 files.

# Script requires PLINK2 and CrossMap v0.7.0 or higher. Can install based on
# directions here: https://crossmap.readthedocs.io/en/latest/#installation
# Chain files and helper code files should be in same directory as this script.

#TODO: change .py before sharing!

# STILL TODO:
# make conda environment and recipe for this?
# implement liftover hg19 to hg38
# update to handle PLINK2 input as well

input_file=$1  # path to and name PLINK1.9 .bim file
out_dir=$2

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
CrossMap.py bed hg38ToHg19.over.chain \
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
    --make-bed --out ${out_dir}/${input_file_name}_hg19

# Report SNP counts
orig_snp_nr=`wc -l ${input_file}`
crossover_snp_nr=`wc -l ${out_dir}/${input_file_name}_hg19.bim`

echo "Original SNP nr: $orig_snp_nr"
echo "Crossovered SNP nr: $crossover_snp_nr"
echo "Original SNP nr: $orig_snp_nr" > ${out_dir}/liftover_plink_log.txt
echo "Crossovered SNP nr: $crossover_snp_nr" >> ${out_dir}/liftover_plink_log.txt

# Cleanup
rm ${out_dir}/tmp_*
