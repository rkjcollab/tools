#!/usr/bin/env bash

set -e
set -u

usage() {
   echo "Usage: ${0##*/} -p <study1_plink_path> -b <study1_build>"
   echo "       -s <study2_plink_path> -t <study2_build> -h <merge_build>"
   echo "       -n <min_maf> -k <min_geno> -s <min_mind> -o <out_dir>"
   echo "Merges two PLINK1.9 datasets after preparing them with"
   echo "prep_for_merge.sh. Each input is lifted over only when its"
   echo "current build differs from merge_build. Then all SNPs are"
   echo "merged, flipped as needed, filtered for allele mismatches,"
   echo "and filtered by missingness and MAF."
}

code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while getopts a:b:c:d:o:m:g:s:t: opt; do
   case "$opt" in
      a) study1_plink_path=${OPTARG%.[bB][eE][dD]} ;;
      b) study1_build=$OPTARG ;;
      c) study2_plink_path=${OPTARG%.[bB][eE][dD]} ;;
      d) study2_build=$OPTARG ;;
      o) out_dir=$OPTARG ;;
      m) min_maf=$OPTARG ;;
      g) min_geno=$OPTARG ;;
      s) min_mind=$OPTARG ;;
      t) merge_build=$OPTARG ;;
      *) usage; exit 1 ;;
   esac
done
shift $((OPTIND - 1))

if [ -z "${study1_plink_path:-}" ] || [ -z "${study1_build:-}" ] || \
   [ -z "${study2_plink_path:-}" ] || [ -z "${study2_build:-}" ] || \
   [ -z "${out_dir:-}" ] || [ -z "${min_maf:-}" ] || \
   [ -z "${min_geno:-}" ] || [ -z "${min_mind:-}" ] || [ -z "${merge_build:-}" ] || [ "$#" -ne 0 ]
then
   usage
   exit 1
fi

mkdir -p "$out_dir"

study1_plink_prefix=$(basename "$study1_plink_path")
study2_plink_prefix=$(basename "$study2_plink_path")
study1_prep_dir="${out_dir}/study1_prep"
study2_prep_dir="${out_dir}/study2_prep"
study1_prepped_prefix="${study1_plink_prefix}_${merge_build}_clean_for_merge"
study2_prepped_prefix="${study2_plink_prefix}_${merge_build}_clean_for_merge"

mkdir -p "$study1_prep_dir" "$study2_prep_dir"

# Prepare each dataset on the same genome build and variant ID scheme.
bash "$code_dir/prep_for_merge.sh" \
   -p "$study1_plink_path" -b "$study1_build" -t "$merge_build" \
   -o "$study1_prep_dir"
bash "$code_dir/prep_for_merge.sh" \
   -p "$study2_plink_path" -b "$study2_build" -t "$merge_build" \
   -o "$study2_prep_dir"

#First merge attempt
plink --bfile "${study1_prep_dir}/${study1_prepped_prefix}" \
   --bmerge "${study2_prep_dir}/${study2_prepped_prefix}" \
   --keep-allele-order \
   --make-bed --out "${out_dir}/tmp_${study1_plink_prefix}_${study2_plink_prefix}_merged1"

#Flip SNPs where necessary
if [ -e "${out_dir}/tmp_${study1_plink_prefix}_${study2_plink_prefix}_merged1-merge.missnp" ]
then
   plink --bfile "${study2_prep_dir}/${study2_prepped_prefix}" \
      --flip "${out_dir}/tmp_${study1_plink_prefix}_${study2_plink_prefix}_merged1-merge.missnp" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study2_prepped_prefix}_flipped"
   plink --bfile "${study1_prep_dir}/${study1_prepped_prefix}" \
      --bmerge "${out_dir}/tmp_${study2_prepped_prefix}_flipped" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2"
else
   plink --bfile "${out_dir}/tmp_${study1_plink_prefix}_${study2_plink_prefix}_merged1" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2"
fi

#Exclude allele mismatch SNPs
if [ -e "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2-merge.missnp" ]
then
   plink --bfile "${study1_prep_dir}/${study1_prepped_prefix}" \
      --exclude "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2-merge.missnp" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study1_prepped_prefix}_no_mismatches"
   plink --bfile "${out_dir}/tmp_${study2_prepped_prefix}_flipped" \
      --exclude "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2-merge.missnp" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study2_prepped_prefix}_no_mismatches"
   plink --bfile "${out_dir}/tmp_${study1_prepped_prefix}_no_mismatches" \
      --bmerge "${out_dir}/tmp_${study2_prepped_prefix}_no_mismatches" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged3"
else
   plink --bfile "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged2" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged3"
fi

#Remove SNPs with large missingness and apply MAF filter if non-zero
if [ "$min_maf" != "0" ]
then
   maf_filter="--maf $min_maf"
else
   maf_filter=""
fi
plink --bfile "${out_dir}/tmp_${study1_prepped_prefix}_${study2_prepped_prefix}_merged3" \
   --geno "$min_geno" \
   --mind "$min_mind" \
   $maf_filter \
   --keep-allele-order \
   --make-bed --out "${out_dir}/${study1_plink_prefix}_${study2_plink_prefix}_merged_clean"

# Clean up
rm "${out_dir}"/tmp_*