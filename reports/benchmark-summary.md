Benchmarking MHcut with different parameters
============================================

MHcut was run using different parameters on 100K randomly selected variants from the dbSNP/ClinVar analysis. We want to see the effect of the following parameters on the numbers we get:

-   Testing guides unicity: Jellyfish vs BLAST.
-   Micro-homology definition: maximum 1 mismatch vs 2 mismatches.
-   Flanks used: best flank vs both flanks.

Jellyfish vs BLAST
------------------

As a quick annotation, MHcut checks if protospacers are unique in the genome. The two ways implemented uses BLAST or JellyFish to test how often a sequence is present in the genome. Here we compare the results for each method on the same set of input variants.

``` r
## For data manipulation
library(data.table) # good for big data
library(dplyr) # To manipulate/transform data.frames
library(knitr)
```

``` r
var.blast = fread('../scripts-dbSNP-ClinVar/dbsnp-clinvar-deletion-100KforBenchmark-1CMM-blast-bestFl-variants.tsv',
                  select=c(20:22,25,29:30))
var.default = fread('../scripts-dbSNP-ClinVar/dbsnp-clinvar-deletion-100KforBenchmark-1CMM-jellyfish-bestFl-variants.tsv',
                  select=c(1:4,20:22,25,29:30))
```

We want to know how many variants have:

1.  MH length &gt;= 3, `mhL>=3`.
2.  an available PAM, `pamMot>0`.
3.  with a unique protospacer, `pamUniq>0`.

The first two numbers should be exactly the same because related to steps upstream of the protospacer mapping.

``` r
blast.sum = tibble(method='BLAST',
                   nb.mh3 = nrow(var.blast[mhL>=3]),
                   nb.mh3.pam = nrow(var.blast[mhL>=3 & pamMot>0]),
                   nb.mh3.pam.uniq = nrow(var.blast[mhL>=3 & pamMot>0 & pamUniq>0]))
jf.sum = tibble(method='JellyFish',
                   nb.mh3 = nrow(var.default[mhL>=3]),
                   nb.mh3.pam = nrow(var.default[mhL>=3 & pamMot>0]),
                   nb.mh3.pam.uniq = nrow(var.default[mhL>=3 & pamMot>0 & pamUniq>0]))

rbind(blast.sum, jf.sum) %>% kable(format.args=list(big.mark=','))
```

| method    |  nb.mh3|  nb.mh3.pam|  nb.mh3.pam.uniq|
|:----------|-------:|-----------:|----------------:|
| BLAST     |  30,150|       5,420|            3,213|
| JellyFish |  30,150|       5,420|            3,267|

Very close numbers of variants with unique guides.

1 mismatch vs 2 mismatches
==========================

By default we allow for at most one mismatch we extending the micro-homology. What if we allow for 2 mismatches? Do we get a lot more micro-homologies?

``` r
var.2mm = fread('../scripts-dbSNP-ClinVar/dbsnp-clinvar-deletion-100KforBenchmark-2CMM-jellyfish-bestFl-variants.tsv',
                  select=c(20:22,25,29:30))

mm1.sum = tibble(method='1 mismatch',
                   nb.mh3 = nrow(var.default[mhL>=3]),
                   nb.mh3.pam = nrow(var.default[mhL>=3 & pamMot>0]),
                   nb.mh3.pam.uniq = nrow(var.default[mhL>=3 & pamMot>0 & pamUniq>0]))
mm2.sum = tibble(method='2 mismatches',
                   nb.mh3 = nrow(var.2mm[mhL>=3]),
                   nb.mh3.pam = nrow(var.2mm[mhL>=3 & pamMot>0]),
                   nb.mh3.pam.uniq = nrow(var.2mm[mhL>=3 & pamMot>0 & pamUniq>0]))

rbind(mm1.sum, mm2.sum) %>% kable(format.args=list(big.mark=','))
```

| method       |  nb.mh3|  nb.mh3.pam|  nb.mh3.pam.uniq|
|:-------------|-------:|-----------:|----------------:|
| 1 mismatch   |  30,150|       5,420|            3,267|
| 2 mismatches |  31,229|       5,873|            3,586|

Only around 1.036 times more micro-homologies found if we allow for 2 consecutive mismatches.

Best flank vs both flanks
-------------------------

When looking for micro-homolgy, MHcut looks at both possible flanks configurations (inner-outer and outer-inner). By default it then picks the ones with highest homology before continuing. However we might be missing some candidates on the other flanks. Maybe the flanks with best homology don't have a valid PAM that we could use but the other one, although with weaker homology, do.

We ran MHcut in both modes: best flank and 2 flanks. Now we want to see how many candidates we gained.

``` r
var.2fls = fread('../scripts-dbSNP-ClinVar/dbsnp-clinvar-deletion-100KforBenchmark-1CMM-jellyfish-2fls-variants.tsv',
                  select=c(1:4, 20:23,27,31:32))
var.2fls$variant = paste(var.2fls$chr, var.2fls$start, var.2fls$end, var.2fls$RS)

fl1.sum = tibble(method='Best flank',
                   nb.mh3 = nrow(var.default[mhL>=3]),
                   nb.mh3.pam = nrow(var.default[mhL>=3 & pamMot>0]),
                   nb.mh3.pam.uniq = nrow(var.default[mhL>=3 & pamMot>0 & pamUniq>0]))
fl2.sum = tibble(method='Both flanks',
                   nb.mh3 = length(unique(var.2fls[mhL>=3]$variant)),
                   nb.mh3.pam = length(unique(var.2fls[mhL>=3 & pamMot>0]$variant)),
                   nb.mh3.pam.uniq = length(unique(var.2fls[mhL>=3 & pamMot>0 & pamUniq>0]$variant)))

rbind(fl1.sum, fl2.sum) %>% kable(format.args=list(big.mark=','))
```

| method      |  nb.mh3|  nb.mh3.pam|  nb.mh3.pam.uniq|
|:------------|-------:|-----------:|----------------:|
| Best flank  |  30,150|       5,420|            3,267|
| Both flanks |  30,151|       5,421|            3,267|

We gain only 1 variant by testing both flanking configurations instead of using the one with highest MH score.
