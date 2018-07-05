MHcut modes: 1 flank vs 2 flanks
================================

The first version of MHcut looks for homology in both possible flanks (inner-outer and outer-inner) and picks the ones with highest homology before continuing. However we might be missing candidates on the other flanks. Maybe the flanks with best homology don't have a valid PAM that we could use but the others, although with weaker homology, do.

I ran MHcut in both modes: 1 flank (i.e. best flank only) and 2 flanks. Now I want to see how many candidates we gained.

First how many MH were found in each run:

| mode     |    MH|
|:---------|-----:|
| 1 flank  |  1332|
| 2 flanks |  1336|

Wow, only 4 variants had MH in both flanks...that's disappointing.

Using the PAM status based on `pamMot`, `pamUniq` and `guidesNoOT` columns, I will count:

-   the number of variants with at least one PAM motif.
-   the number of variants with at least one PAM motif and unique guide
-   the number of variants with at least one PAM motif, a unique guide and no off-target MH.

| mode     |  pam|  pam.uniq|  pam.uniq.noOt|
|:---------|----:|---------:|--------------:|
| 1 flank  |  487|       359|            109|
| 2 flanks |  488|       360|            109|

In the end we gained nothing, not surprising considering the number of additional variants that the "2 flanks" mode gives.

Sanity check: the "1 flank" run is a subset of the "2 flanks" run
-----------------------------------------------------------------

| mode     |  pam|  pam.uniq|  pam.uniq.noOt|
|:---------|----:|---------:|--------------:|
| 1 flank  |  487|       359|            109|
| 2 flanks |  488|       360|            109|

What !?

Hum, could it be that some variants from the "2 flanks" mode are not included in the other run?

| Chromosome |      Start|       Stop| pam.1 | pam.uniq.1 | pam.uniq.noOt.1 | pam.2 | pam.uniq.2 | pam.uniq.noOt.2 |
|:-----------|----------:|----------:|:------|:-----------|:----------------|:------|:-----------|:----------------|
| chr2       |   47402369|   47422413| NA    | NA         | NA              | FALSE | FALSE      | FALSE           |
| chrX       |  108687509|  108687535| NA    | NA         | NA              | TRUE  | TRUE       | FALSE           |

Ah, these two variants were not included in the "1 flank" mode results.

| Chromosome |      Start|       Stop| GeneSymbol |      dbSNP| dbVar | reviewStatus                   |  nbSubmitters|   varL|  flank|  flankScore|  mhL|  mh1L|   hom|  nbMM|  mhDist| MHseq1         | MHseq2         |  pamMot|  bestPamHet|  pamUniq|  guidesNoOT|  guidesMinOT| size.class |
|:-----------|----------:|----------:|:-----------|----------:|:------|:-------------------------------|-------------:|------:|------:|-----------:|----:|-----:|-----:|-----:|-------:|:---------------|:---------------|-------:|-----------:|--------:|-----------:|------------:|:-----------|
| chr2       |   47402369|   47422413| MSH2       |         -1| -     | reviewed by expert panel       |             1|  20045|      1|          22|   14|     9|  0.93|     1|   20031| aaccctccgactcc | aaccctccggctcc |       0|          NA|        0|          NA|           NA| \>30       |
| chrX       |  108687509|  108687535| COL4A5     |  104886412| -     | no assertion criteria provided |             1|     27|      2|           8|    4|     4|  1.00|     0|      23| GTCC           | GTCC           |       6|          18|        5|           0|            5| 11-30      |

One is a large deletion, the other is more surprising. I checked the cartoon for the second one and it's because the best pair of flanks were filtered out right away, most likely by the `m1L>=3` default filter. We might want to first check for these criteria and then chose the pair of flanks. Or change the scoring system that we use to choose a pair of flanks.
