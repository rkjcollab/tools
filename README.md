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
