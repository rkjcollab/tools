# tools

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

Initial version developed by Shane E. Ridoux.

*TODOs:*

+ add documentation/use instructions here
+ revisit and collapse 1-2 & 1-3
+ make more files TMP only
+ confirm MAF 0.1 at time LD prune merged input
+ remove study-only Admixture step & shorten final script
