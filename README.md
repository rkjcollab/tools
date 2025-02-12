# liftover
Pipeline to run liftover using CrossMap - in development. Currently can do hg38 to hg19 for PLINK1.9 or bed input.

Examples:

``` bash
bash liftover_plink.sh daisyask.bim /Users/slacksa/support
```

``` bash
bash liftover_bed.sh T1DGC_saige_EUR_id.txt 1 2 14 /Users/slacksa/support
```