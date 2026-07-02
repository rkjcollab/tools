#!/usr/bin/env bash

set -e
set -u

usage() {
   echo "Usage: ${0##*/} -p <plink_path> -b <from_build (hg19 or hg38)>"
   echo "      -t <to_build (hg19 or hg38)> -o <out_dir>"
   echo "Cleans PLINK1.9 data using PLINK2 in preparation of merge:"
   echo "removes SNPs with duplicate positions, lifts over from current"
   echo "build to desired build, updates varIDs to chr:pos:ref:alt, and"
   echo "removes SNPs with duplicate IDs. If no liftover is needed, set"
   echo "from_build and to_build inputs to the same value."
}

code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while getopts p:o:n:b:t:s:k:h:l: opt; do
   case "$opt" in
      p) plink_path=${OPTARG%.[bB][eE][dD]} ;;
      o) out_dir=$OPTARG ;;
      b) from_build=$OPTARG ;;
      t) to_build=$OPTARG ;;
      *) usage; exit 1 ;;
   esac
done
shift $((OPTIND - 1))

if [ -z "${plink_path:-}" ] || [ -z "${out_dir:-}" ] || \
   [ -z "${from_build:-}" ] || [ -z "${to_build:-}" ] || [ "$#" -ne 0 ]
then
   usage
   exit 1
fi

plink_prefix=$(basename "$plink_path")

mkdir -p "$out_dir"

# Deduplicate and set all variant IDs to chr:pos:ref:alt
# First, if duplicate IDs have different missingness, remove the SNP with
# more missingness.
plink2 --bfile "$plink_path" \
   --missing --freq \
   --nonfounders \
   --make-pgen \
   --out "${out_dir}/tmp_dedup"
Rscript --vanilla "${code_dir}/dedup_miss.R" \
   -v "${out_dir}/tmp_dedup.vmiss" \
   -p "${out_dir}/tmp_dedup.pvar"
plink2 --bfile "$plink_path" \
   --make-bed \
   --exclude "${out_dir}/tmp_dedup_rm.txt" \
   --out "${out_dir}/tmp_dedup_rm"

# Then, update all IDs to chr:pos:ref:alt and remove the remaining duplicates
# by arbitrarily keeping first
plink2 --bfile "${out_dir}/tmp_dedup_rm" \
   --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000 \
   --rm-dup force-first \
   --make-bed \
   --out "${out_dir}/tmp_no_dupl"

# Create bed file to crossover from hg19 to hg38 
cut -f1 "${out_dir}/tmp_no_dupl.bim" | sed 's/^/chr/' > "${out_dir}/tmp_c1.txt"
cut -f4 "${out_dir}/tmp_no_dupl.bim" > "${out_dir}/tmp_c2.txt"
cut -f4 "${out_dir}/tmp_no_dupl.bim" > "${out_dir}/tmp_c3.txt"
cut -f2 "${out_dir}/tmp_no_dupl.bim" > "${out_dir}/tmp_c4.txt"
paste  "${out_dir}/tmp_c1.txt" \
       "${out_dir}/tmp_c2.txt" \
       "${out_dir}/tmp_c3.txt" \
       "${out_dir}/tmp_c4.txt" \
       >  "${out_dir}/tmp_in.bed"

# Check if liftover requested and do liftover
if [ "$from_build" == "$to_build" ]; then
   echo "Not lifting over, input is already desired ${to_build}"

elif [[ "$to_build" =~ .*38.* ]]; then
   CrossMap.py bed "${code_dir}/hg19ToHg38.over.chain" \
      "${out_dir}/tmp_in.bed"  \
      "${out_dir}/tmp_out.bed"

elif [[ "$to_build" =~ .*19.* ]]; then
   CrossMap.py bed "${code_dir}/hg38ToHg19.over.chain" \
      "${out_dir}/tmp_in.bed"  \
      "${out_dir}/tmp_out.bed"
fi

# Follow up liftover
if [ "$from_build" == "$to_build" ]; then
   plink2 --bfile "${out_dir}/tmp_no_dupl" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_gwas"

else
   # Extract only those SNPs that were successfully cross-overed
   cut -f4 "${out_dir}/tmp_out.bed" > "${out_dir}/tmp_snp_keep.txt"
   plink2 --bfile "${out_dir}/tmp_no_dupl" \
      --extract "${out_dir}/tmp_snp_keep.txt" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/tmp_gwas"

   # Update bim file positions
   Rscript --vanilla "${code_dir}/update_pos.R" \
      --crossmap-bed "${out_dir}/tmp_out.bed" \
      --success-bim "${out_dir}/tmp_gwas.bim"
fi

# Update SNP IDs, remove ID dups, sort output (need to make .pgen)
plink2 --bfile "${out_dir}/tmp_gwas" \
   --set-all-var-ids @:#:\$r:\$a --rm-dup --sort-vars \
   --make-pgen --out "${out_dir}/tmp_output"

# Make final file
plink2 --pfile "${out_dir}/tmp_output" \
   --keep-allele-order \
   --make-bed --out "${out_dir}/${plink_prefix}_${to_build}_clean_for_merge"

# Clean up
rm "${out_dir}"/tmp_*
