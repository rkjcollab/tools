#!/bin/bash

# Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <input_file> <chr_col_num> <pos_col_num>"
   echo "       <id_col_num> <out_dir>"
   exit
fi

# Liftover from hg38 to hg19 with tab/space (txt, tsv) delimited file as input.
# Returns a file in simple bed format: columns chrom pos pos.

# Script requires CrossMap v0.7.0 or higher. Can install based on
# directions here: https://crossmap.readthedocs.io/en/latest/#installation
# Chain files and helper code files should be in same directory as this script.

# STILL TODO:
# make conda environment and recipe for this?
# implement liftover hg19 to hg38

input_file=$1  # path to and name of file with extension
chr_col_num=$2
pos_col_num=$3
id_col_num=$4  # to allow for mapping back
out_dir=$5

# Get input_file without path and without extension
input_file_no_path=$(basename "$input_file")
input_file_name=${input_file_no_path%.*}

# TO NOTE: If need to create id column from chr:pos:ref:alt, can use code like
# this and make modified input file before running rest of script:
# cat $input_file | cut -f$chr_col_num > ${out_dir}/tmp_c1.txt  # -f# number of column with chromosome
# cat $input_file | cut -f$pos_col_num > ${out_dir}/tmp_c2.txt  # -f# number of column with position
# cat $input_file | cut -f3 > ${out_dir}/tmp_c3.txt  # -f# number of column with ref allele
# cat $input_file | cut -f4 > ${out_dir}/tmp_c4.txt  # -f# number of column with alt allele
# awk '{print $1 ":" $2 ":" $3 ":" $4}' $input_file > ${out_dir}/tmp_id.txt  # chr:pos ID column
# paste $input_file ${out_dir}/tmp_id.txt > ${out_dir}/${input_file_name}_id.txt

# Extract chr and pos columns to simple bed format
cat $input_file | cut -f$chr_col_num > ${out_dir}/tmp_c1.txt
cat $input_file | cut -f$pos_col_num > ${out_dir}/tmp_c2.txt
cat $input_file | cut -f$id_col_num > ${out_dir}/tmp_c3.txt

paste ${out_dir}/tmp_c1.txt \
    ${out_dir}/tmp_c2.txt \
    ${out_dir}/tmp_c2.txt \
    ${out_dir}/tmp_c3.txt \
    > ${out_dir}/tmp_bed.txt

# Check if input file had column names
first_c2=$(awk 'NR == 1 {print $1}' ${out_dir}/tmp_c2.txt)
if [[ ! "$first_c2" =~ ^[0-9]+$ ]]
then
    echo "Input file had column names. Removing for liftover."
    tail -n +2 ${out_dir}/tmp_bed.txt > ${out_dir}/${input_file_name}_bed.bed
else
    cp ${out_dir}/tmp_bed.txt ${out_dir}/${input_file_name}_bed.bed
fi

# Do crossover
CrossMap bed hg38ToHg19.over.chain \
   ${out_dir}/${input_file_name}_bed.bed  \
   ${out_dir}/${input_file_name}_bed_hg19.bed

# Report SNP counts
orig_snp_nr=`wc -l ${input_file}`
crossover_snp_nr=`wc -l ${out_dir}/${input_file_name}_bed_hg19.bed`

echo "Original SNP nr: $orig_snp_nr"
echo "Crossovered SNP nr: $crossover_snp_nr"
echo "Original SNP nr: $orig_snp_nr" > ${out_dir}/liftover_bed_log.txt
echo "Crossovered SNP nr: $crossover_snp_nr" >> ${out_dir}/liftover_bed_log.txt

# Cleanup
rm ${out_dir}/tmp_*
