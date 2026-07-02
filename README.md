# merge_genetic_data

Scripts for preparing and merging two PLINK datasets.

## Requirements

- PLINK 1.9
- PLINK 2
- CrossMap (only needed when lifting between genome builds)
- R
- R package `argparse`

The conda enviroment included the the RKJcollab inmputation repo has all the dependencies.

## How the scripts fit together

`prep_for_merge.sh` is the per-dataset preparation step. It optionally lifts over from `hg19` to `hg38` or from `hg38` to `hg19`, standardizes variant IDs to `chr:pos:ref:alt`, removes variants with duplicate chr:pos:ref:alt IDs, and writes a PLINK bed/bim/fam set ready for merging.

`merge_data_sets.sh` calls `prep_for_merge.sh` for each input then merges the prepared outputs.

## Prepare one dataset for merging

```bash
./prep_for_merge.sh -p <plink_prefix> -b <from_build> -t <to_build> -o <out_dir>

```

`from_build` and `to_build` should be `hg19` or `hg38`. If the builds match, no liftOver is performed. The script always uses the chain files in this repository. If doing liftOver, this script calls helper `update_pos.R`.

Output:

```text
<out_dir>/<plink_prefix>_<to_build>_clean_for_merge.{bed,bim,fam}

```

## Merge two datasets

```bash
./merge_data_sets.sh -p <study1.bed> -b <study1_build> -s <study2.bed> -t <study2_build> -h <merge_build> -n <min_maf> -k <min_geno> -o <out_dir>

```

This script calls `prep_for_merge.sh` for each input. If an input build already matches `merge_build`, no liftOver is performed for that dataset. After preparation, the script flips SNPs when needed, excludes remaining allele mismatches, and applies the final missingness and MAF filters.

Output:

```text
<out_dir>/<study1>_<study2>_merged_clean.{bed,bim,fam}

```

## Helper scripts

- `update_pos.R`: updates BIM positions after CrossMap liftOver.
