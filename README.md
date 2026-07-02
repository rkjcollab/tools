<<<<<<< HEAD
# tools

*TO NOTE: these tools are under active development.*

Small tools/pipelines commonly used by the RKJcollab.

### **liftover**

Pipeline to run liftover using CrossMap - in development. Currently can do hg38 to hg19 for PLINK1.9 or bed input.

Examples:

``` bash
bash liftover_plink.sh daisyask.bim /Users/slacksa/support
```

``` bash
bash liftover_bed.sh T1DGC_saige_EUR_id.txt 1 2 14 /Users/slacksa/support
```

### **admixture**

Initial version developed by Shane E. Ridoux. Currently implemented to identify primarily European individuals.
Default is to use 1000 Genomes (1000G) Phase 3 Version 5a as the reference population.

After restricting to common variants, 10 genetic principal components (PCs) are estimated in 1000G, and study is
projected into the 1000G PC space. Then, a gaussian mixture model with G=5 clusters for the five 1000G ancestral
populations is used on the first four 1000G PCs (Budiarto A, et al. 2021 Procedia Computer Science), with each
cluster labeled by its predominant 1000G super-population. Study individuals that were predicted to map to the
1000G European population using the first four PCs were included in European list. To evaluate population
separation, the pipeline also uses unsupervised Admixture with five-fold cross validation in 1000G to choose K=5
ancestral populations and estimate population structure, before using Admixture to project the study onto 1000G.


*TODOs:*

+ convert methods above into documentation and use instructions
+ revisit and collapse steps 1-2 & 1-3
+ make more files TMP
+ confirm MAF 0.1 at time LD prune merged input
=======
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
- `dedup_miss.R`: chooses duplicates to remove by missingness.
>>>>>>> merge_genetic_data/master
